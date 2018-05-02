#Process species attribute data
#Jens Stevens stevensjt@gmail.com

library(tidyverse)
library(readxl)
library(Hmisc)
library(gridExtra)
stdErr <- function(x) sqrt(var(x)/length(x))

####1. Read data####
#sp.d <- #Species list data (first pass analysis)
#  read_excel("./Data/Raw/tc understory species frequency pre post hayman 1997-2012.xlsx",
#             sheet = "data")
sp.d <- #Percent cover data. Warnings ok (have to do with NA's)
  read_excel("Data/Raw/tc understory species cover hayman 1997-2012 for jens.xlsx",
             sheet = "data",
             col_types = c(rep("text",2),"numeric",rep("text",7),rep("numeric",7)))
#sp.d[1,c(11:17)] <- round(as.numeric(sp.d[,c(11:17)]),2) #Deprecated
sp.attr <- 
  read_excel("./Data/Raw/RavenAxelrodOrigins_JTS.xlsx")

####2. Assign species origins####

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
rm(sp.attr)



####3. Additional data processing####
#NA in spp cover data means not present in macroplot
#0 in spp cover data means present in macroplot but not subplot.
sp.d[,c(11:17)][is.na(sp.d[,c(11:17)])] <-0 #Set "NA" to 0 (we are focusing on subplots)
absent_spp <- which(rowSums(sp.d[,c(11:17)])==0) #Find speciesXplots where spp is absent
sp.d <- sp.d[-absent_spp,] #Remove 1125 absent speciesXplots

sp.d[grep("Triticosecale",sp.d$SciName),"origin"] <- NA #Deal with wierd hybrid species.
sp.d[grep("Unk",sp.d$origin),"origin"] <- NA #Reclassify 94 Unk speciesXplots as NA rather than Non-NTM (14 species)

missing.sp <- sp.d %>%
  group_by(SciName) %>%
  summarise(origin = origin[1])
missing.sp[grep("Triticosecale",missing.sp$SciName),"origin"] <- NA #Have to do this twice for some reason
length(which(is.na(missing.sp$origin)))/nrow(missing.sp) #195 species, have origin data for 149 (76%)
#write_csv(missing.sp,"./Data/Derived/HaymanOrigins_JTS.csv")

sp.d[,"origin_binary"] <- sp.d[,"origin"]
sp.d[,"origin_binary"][sp.d[,"origin_binary"]!="NTm"] <- "Non-NTm"
#sp.d=sp.d[-which(is.na(sp.d$origin)),] #Remove unknown origins

####4. Group into different datasets####


sp.d.long=gather(sp.d,key = year,value = cover, -Plot, -FireSeverity, -FireMortality,
                 -TopoClass, -Direction, -FinalCode, -SciName, -Family, -FunctionalGrp, 
                 -NativeStatus,-origin,-origin_binary)

sp.d.long$year=as.integer(sp.d.long$year)
sp.d.long$presence = ifelse(sp.d.long$cover>0,1,0)
sp.d.long = sp.d.long[-which(is.na(sp.d.long$origin_binary)),] #Remove unknown origins
sp.d.long = sp.d.long[-which(sp.d.long$TopoClass=="Riparian"),] #Remove unknown origins
sp.d.long$origin_binary[sp.d.long$origin_binary=="NTm"] = "Northern"
sp.d.long$origin_binary[sp.d.long$origin_binary=="Non-NTm"] = "Southern"

sp.d.long <- #identify colonizations/extinctions
  group_by(sp.d.long,SciName,Plot) %>%
  mutate(colonize = ifelse(presence==1, 
                           ifelse(lag(presence)==0,1,0),NA),
         extinct = ifelse(presence==0,
                          ifelse(lag(presence)==1,1,0),NA),
         abs.colonize = #Absolute colonizer, 1 if present but absent in 1997
           ifelse(presence==1&year!=1997, 
                               ifelse(presence[1]==0,1,0),NA),
         abs.extinct = #Absolute extinction, 1 if absent but present in 1997
           ifelse(presence==0&year!=1997, 
                              ifelse(presence[1]==1,1,0),NA),
         all.abs.colonize = 
           #Absolute colonizer, 1 if absent in 1997, colonized in current year,
           #present in 2012, current year is not 2012, and only colonized once.
           ifelse(presence[1]==0 & 
                    colonize==1 & presence[length(presence)]==1 & year!=2012 &
                    sum(colonize,na.rm=T)==1, 1, 0),
         all.abs.extinct = 
           #Absolute extinction, 1 if present in 1997, extinct in current year,
           #absent in 2012, current year is not 2012, and only extinct once.
           ifelse(presence[1]==1 &
                    extinct==1 & presence[length(presence)]==0 & year!=2012 &
                    sum(extinct,na.rm=T)==1, 1, 0)
         )
which.extinct <- sp.d.long[which(sp.d.long$all.abs.extinct == 1),]
p.d <- sp.d.long %>% #Data with plots separate
  group_by(Plot, FireSeverity, origin_binary, year) %>% 
  summarise(mean_cover = mean(cover),sd=sd(cover),
            richness = length(which(cover>0)),
            richness2 = sum(presence),
            colonizations = sum(lag(colonize==1),na.rm=T),
            extinctions = sum(lag(extinct==1),na.rm=T),
            abs_colonizations = sum(all.abs.colonize,na.rm=T),
            abs_extinctions = sum(all.abs.extinct,na.rm=T))

p.d.ratio <- p.d %>% #Plot data with ratios
  group_by(Plot, FireSeverity, year) %>% 
  summarise(Prop.NTM = richness[2]/sum(richness),
            richness.tot = sum(richness),
            richness.NTM = richness[2],
            richness.Non_NTM = richness[1])

sev.class.grp <- sp.d.long %>% #Data with all plots grouped into severity class
  group_by(FireSeverity, origin_binary,year) %>% 
  summarise(mean_cover = mean(cover), sd = sd(cover),
            se = stdErr(cover),
            n_plots = length(unique(Plot))  )

sev.class.grp.ratio <- p.d.ratio %>% #Data with all plots grouped into severity class
  #Plus NTM ratio
  group_by(FireSeverity,year) %>% 
  summarise(mean_Prop.NTM = mean(Prop.NTM), sd = sd(Prop.NTM),
            se = stdErr(Prop.NTM),
            n_plots = length(unique(Plot))  
            )


####5. Plots- cover####
F1a_LowSev <- ggplot(
  sev.class.grp[!is.na(sev.class.grp$origin_binary) & 
                  sev.class.grp$FireSeverity == "Low" ,]
  ) +
  geom_point(aes(x=year, y=mean_cover, col = origin_binary)) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("blue","darkred"))+
  geom_errorbar(aes(x=year, ymin = mean_cover - se, ymax = mean_cover + se, 
                    col = origin_binary))+
  geom_line(aes(x=year, y=mean_cover, col = origin_binary))+
  labs(title = "Low severity", y = "mean % cover", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none")

F1b_ModSev <- ggplot(
  sev.class.grp[!is.na(sev.class.grp$origin_binary) & 
                  sev.class.grp$FireSeverity == "Moderate" ,]
) +
  geom_point(aes(x=year, y=mean_cover, col = origin_binary)) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("blue","darkred"))+
  geom_errorbar(aes(x=year, ymin = mean_cover - se, ymax = mean_cover + se, 
                    col = origin_binary))+
  geom_line(aes(x=year, y=mean_cover, col = origin_binary))+
  labs(title = "Moderate severity", y = "mean % cover", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

F1c_HighSev <- ggplot(
  sev.class.grp[!is.na(sev.class.grp$origin_binary) & 
                  sev.class.grp$FireSeverity == "High" ,]
) +
  geom_point(aes(x=year, y=mean_cover, col = origin_binary)) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("blue","darkred"))+
  geom_errorbar(aes(x=year, ymin = mean_cover - se, ymax = mean_cover + se, 
                    col = origin_binary))+
  geom_line(aes(x=year, y=mean_cover, col = origin_binary))+
  labs(title = "High severity", y = "mean % cover", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.78, 0.91))

pdf(file = paste0("./Figures/NAC/Fig1Horiz",Sys.Date(),".pdf"),width=10,height=7)
grid.arrange(F1a_LowSev,F1b_ModSev,F1c_HighSev,nrow=1)
dev.off()


####6. Plots- ratio####
pdf(file = paste0("./Figures/NAC/Fig3_propNTM",Sys.Date(),".pdf"),width=10,height=7)
ggplot(p.d.ratio,
       aes(x=year,y=Prop.NTM,col=FireSeverity))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title = "proportion of flora \nwitn north-temperate affinity", y= "proportion",
     col = "Fire severity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


summary(lm(Prop.NTM~year, data=p.d.ratio[p.d.ratio$FireSeverity=="High",]))
summary(lm(Prop.NTM~year, data=p.d.ratio[p.d.ratio$FireSeverity=="Moderate",]))
summary(lm(Prop.NTM~year, data=p.d.ratio[p.d.ratio$FireSeverity=="Low",]))

pdf(file = paste0("./Figures/NAC/Fig4_meanpropNTM",Sys.Date(),".pdf"),width=10,height=7)
ggplot(sev.class.grp.ratio)+
  geom_point(aes(x=year, y=mean_Prop.NTM, col = FireSeverity)) + 
  geom_errorbar(aes(x=year, ymin = mean_Prop.NTM - se, ymax = mean_Prop.NTM + se, 
                    col = FireSeverity))+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  geom_line(aes(x=year, y=mean_Prop.NTM, col = FireSeverity))+
  labs(title = "mean proportion of flora \nwitn north-temperate affinity", y= "proportion",
     col = "Fire severity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

####7. Plots- colonizations/extinctions####
pdf(file = paste0("./Figures/NAC/Fig5_Colonization_Extinction",Sys.Date(),".pdf"),width=10,height=7)
ggplot(p.d,aes(col=origin_binary,fill=FireSeverity)) +
  geom_bar(aes(x=year,y=-abs_extinctions),
           position = "dodge", stat = "summary", fun.y = "mean") + 
  geom_bar(aes(x=year,y=abs_colonizations),
           position = "dodge", stat = "summary", fun.y = "mean") +
  stat_summary(aes(year,abs_colonizations),fun.data = mean_se, geom = "errorbar",position="dodge")+ #Needs work
  stat_summary(aes(year,-abs_extinctions),fun.data = mean_se, geom = "errorbar",position="dodge")+ #Needs work
  scale_colour_manual(values=c("gray","black"))+
  geom_hline(aes(yintercept=0))+
  labs(title="permanent species colonizations or extinctions \nfrom pre-fire",
       y = "number of species")+
  xlim(c(2002,2008))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()
