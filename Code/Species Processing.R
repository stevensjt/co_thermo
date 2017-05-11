#Process species attribute data
#Jens Stevens stevensjt@gmail.com

library(tidyverse)
library(readxl)
library(Hmisc)

####1. Read data####
sp.d <- #Species list data
  read_excel("./Data/Raw/tc understory species frequency pre post hayman 1997-2012.xlsx",
             sheet = "data")
sp.d <- #Percent cover data. Warnings ok (have to do with NA's)
  read_excel("Data/Raw/tc understory species cover hayman 1997-2012 for jens.xlsx",
             sheet = "data",
             col_types = c(rep("text",2),"numeric",rep("text",7),rep("numeric",7)))
#sp.d[1,c(11:17)] <- round(as.numeric(sp.d[,c(11:17)]),2) #Deprecated
sp.attr <- 
  read_excel("./Data/Raw/RavenAxelrodOrigins_JTS.xlsx")

####2. Process data####

#2a: Get origins for species with exact name matches in the RnA data
sp.d$SciName <- gsub('\\s+','_',sp.d$SciName ) #Replace space with underscore for consistency
index <- pmatch(sp.d$SciName,sp.attr$mrt) #Identify position matches in sp.d
sp.d[,"origin"] <- sp.attr[index,"RnA"]

#2b: Get origins for species with genera having one and only one origin in the RnA data, but the species itself is not represented
genera <- gsub("_.*$","",sp.d$SciName)
genera[!is.na(index)] <- NA #Exclude species you've already classified.
for(g in 1:length(genera)){
  n <-
    genera[g]%>%
    grep(x=sp.attr$mrt)
  if(!NA%in%n & length(unique(sp.attr$RnA[n] ) )==1 ){
    sp.d[g,"origin"] <- 
      sp.attr$RnA[n] %>%
      unique() %>%
      paste()
  }
}
length(which(!is.na(sp.d$origin)))/length(sp.d$origin) #Have origin data for 86% of species at this point, if you run this on species list data
missing.sp <- sp.d %>%
  group_by(SciName) %>%
  summarise(origin = origin[1])
#write_csv(missing.sp,"./Data/Derived/HaymanOrigins_JTS.csv")
#sp.d[,6:12][is.na(sp.d[,6:12])] <- 0 #only relevant for Species list data
#CHECKME Need to reclassify Unk as NA rather than Non-NTM


####3. Exploratory analysis on Species list data####
sp.d[,"fire_change"] <- (sp.d[,"2003"] - sp.d[,"1997"])
sp.d[,"origin_binary"] <- sp.d[,"origin"]
sp.d[,"origin_binary"][sp.d[,"origin_binary"]!="NTm"] <- "Non-NTm"
sp.d=sp.d[-which(is.na(sp.d$origin)),] #Remove unknown origins
ggplot(sp.d,aes(x=origin_binary,y=fire_change))+
  geom_point()+
  stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 0.2)

#START HERE
sp.d.long=gather(sp.d,year,cover,-Family,-SciName,-FinalCode,-FunctionalGrp,-NativeStatus,-origin,-fire_change,-origin_binary)
sp.d.long$year=as.integer(sp.d.long$year)
sp.d.grp <- sp.d.long %>% 
  group_by(origin_binary,year) %>% 
  summarise(cover = mean(cover),sd=sd(cover))

ggplot(sp.d.long,aes(x=year,y=sqrt(cover),col=origin_binary))+
  geom_line(aes(group=SciName))+
  geom_line(data=sp.d.grp,alpha = .8, size = 3)


####3b. Exploratory analysis on percent cover data

#3b.1 Filter data
sp.d[,"origin_binary"] <- sp.d[,"origin"]
sp.d[,"origin_binary"][sp.d[,"origin_binary"]!="NTm"] <- "Non-NTm"
sp.d$FireSeverity <- factor(sp.d$FireSeverity,levels=c("High","Moderate","Low"))

pc.d.long <- #Create long-form data frame; cover data is from year columns (last argument)
  gather(sp.d,year,cover, which(!is.na(as.integer(names(sp.d)))))
pc.d.long$year=as.integer(pc.d.long$year)
pc.d.long <- pc.d.long[-which(is.na(pc.d.long$origin_binary)),] #Remove unclassified species

pc.d.long <-
  group_by(pc.d.long,SciName,year,Plot) %>%
  mutate(presence = ifelse(is.na(cover),0,1)) %>%
  group_by(SciName,Plot) %>%
  mutate(colonize = ifelse(presence==1, 
                           ifelse(lag(presence)==0,1,0),NA),
         extinct = ifelse(presence==0,
                          ifelse(lag(presence)==1,1,0),NA),
         abs.colonize = ifelse(presence==1&year!=1997, #Absolute colonizer, 1 if present but absent in 1997
                               ifelse(presence[1]==0,1,0),NA),
         abs.extinct = ifelse(presence==0&year!=1997, #Absolute extinction, 1 if absent but present in 1997
                               ifelse(presence[1]==1,1,0),NA),
         all.abs.colonize = ifelse(presence[1]==0 & #Absolute colonizer, 1 if absent in 1997, colonized in current year,
                                     #present in 2012, current year is not 2012, and only colonized once.
                                     colonize==1 & presence[length(presence)]==1 & year!=2012 &
                                     sum(colonize,na.rm=T)==1, 1, 0),
         all.abs.extinct = ifelse(presence[1]==1 & #Absolute extinction, 1 if present in 1997, extinct in current year,
                                     #absent in 2012, current year is not 2012, and only extinct once.
                                     extinct==1 & presence[length(presence)]==0 & year!=2012 &
                                     sum(extinct,na.rm=T)==1, 1, 0)
         )
  
#START HERE
pc.d.grp <- pc.d.long %>% 
  group_by(Plot, FireSeverity, FireMortality, TopoClass, Direction, origin_binary,year) %>% 
  summarise(tot_cover = sum(cover, na.rm=T),
            sd = sd(cover,na.rm=T),
            richness = sum(presence), 
            colonizations = sum(lag(colonize==1),na.rm=T),
            extinctions = sum(lag(extinct==1),na.rm=T),
            abs_colonizations = sum(all.abs.colonize,na.rm=T),
            abs_extinctions = sum(all.abs.extinct,na.rm=T))



pc.d.ratio <- pc.d.grp %>%
  group_by(Plot, FireSeverity, FireMortality, TopoClass, Direction, year) %>% 
  summarise(Prop.NTM = richness[2]/sum(richness),
            richness.tot = sum(richness),
            richness.NTM = richness[2],
            richness.Non_NTM = richness[1])



#Analysis of NTM data
#png(file = paste0("./Figures/EDA/Fig1_PropNTM_Time_",Sys.Date(),".png"),width=8,height=12,units="in",res=200)

ggplot(pc.d.ratio[pc.d.ratio$TopoClass!="Riparian",],
       aes(x=year,y=Prop.NTM,col=FireSeverity))+
  geom_point()+
  geom_smooth(method="lm")
#dev.off()
summary(lm(Prop.NTM~year*FireSeverity, data=pc.d.ratio[pc.d.ratio$TopoClass!="Riparian",]))
summary(lm(Prop.NTM~year, data=pc.d.ratio[pc.d.ratio$FireSeverity=="High",]))
summary(lm(Prop.NTM~year, data=pc.d.ratio[pc.d.ratio$FireSeverity=="Moderate",]))
summary(lm(Prop.NTM~year, data=pc.d.ratio[pc.d.ratio$FireSeverity=="Low",]))

png(file = paste0("./Figures/EDA/Fig2_RelativeChange",Sys.Date(),".png"),width=8,height=12,units="in",res=200)
ggplot(pc.d.grp[pc.d.grp$TopoClass!="Riparian",],aes(col=origin_binary,fill=FireSeverity)) +
  geom_bar(aes(x=year,y=-extinctions),
           position = "dodge", stat = "summary", fun.y = "mean") + 
  geom_bar(aes(x=year,y=colonizations),
           position = "dodge", stat = "summary", fun.y = "mean") +
  stat_summary(aes(year,colonizations),fun.data = mean_se, geom = "errorbar",position="dodge")+ #Needs work
  stat_summary(aes(year,-extinctions),fun.data = mean_se, geom = "errorbar",position="dodge")+ #Needs work
  scale_colour_manual(values=c("gray","black"))+
  geom_hline(aes(yintercept=0))+
  labs(y="species colonizations or extinctions from previous year")+
  xlim(c(2002,2008))
dev.off()

png(file = paste0("./Figures/EDA/Fig3_PermanentChange",Sys.Date(),".png"),width=8,height=12,units="in",res=200)
ggplot(pc.d.grp[pc.d.grp$TopoClass!="Riparian",],aes(col=origin_binary,fill=FireSeverity)) +
  geom_bar(aes(x=year,y=-abs_extinctions),
           position = "dodge", stat = "summary", fun.y = "mean") + 
  geom_bar(aes(x=year,y=abs_colonizations),
           position = "dodge", stat = "summary", fun.y = "mean") +
  stat_summary(aes(year,abs_colonizations),fun.data = mean_se, geom = "errorbar",position="dodge")+ #Needs work
  stat_summary(aes(year,-abs_extinctions),fun.data = mean_se, geom = "errorbar",position="dodge")+ #Needs work
  scale_colour_manual(values=c("gray","black"))+
  geom_hline(aes(yintercept=0))+
  labs(y="permanent species colonizations or extinctions from pre-fire")+
  xlim(c(2002,2008))
dev.off()

rel.m=perm.m=list(); l=1
for(y in 2003:2007){
  rel.m[[l]]<- summary(lm(colonizations/(extinctions+0.5)~FireSeverity+origin_binary,
                          data=pc.d.grp[pc.d.grp$year==y & pc.d.grp$TopoClass!="Riparian",]))
  perm.m[[l]]<- summary(lm(abs_colonizations/(abs_extinctions+0.5)~FireSeverity+origin_binary,
                          data=pc.d.grp[pc.d.grp$year==y & pc.d.grp$TopoClass!="Riparian",]))
  l=l+1
}
perm.m
