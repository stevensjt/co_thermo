

####0. Read libraries####
library(readxl) ##version 1.0.0; for read_excel()
library(tidyverse) ##version 1.2.1; for read_csv (from readr v 1.1.1)
library(gridExtra) #version 2.3; for grid.arrange()
#library(Hmisc) #for mean_se, maybe deprecated?
stdErr <- function(x) sqrt(var(x, na.rm = T)/length(na.exclude(x)))

####1. Read and process data####
##Read component data
sp.d <- #Percent cover data. Warnings ok (have to do with NA's)
  read_excel("Data/Raw/tc understory species cover hayman 1997-2012 for jens.xlsx",
             sheet = "data",
             col_types = c(rep("text",2),"numeric",rep("text",7),rep("numeric",7)))
sp.d$SciName <-  #Replace space w/underscore for consistency w origin data
  gsub('\\s+','_',sp.d$SciName )
sp.attr <- #Species origin data
  read_csv("./Data/Raw/HaymanOrigins.csv")
env.d <- #Read other environmental data
    read_excel("Data/Raw/tc plot characteristics.xlsx", sheet = "Sheet1")


##Remove bad species
sp.d <- sp.d[-grep("Triticosecale",sp.d$SciName),] #Sterile hybrid; N = 30 
sp.attr <- sp.attr[-grep("Triticosecale",sp.attr$SciName),] #Sterile hybrid; N = 2
sp.d <- sp.d[-grep("unable",sp.d$SciName),] #Unable to identify; N = 23
sp.attr <- sp.attr[-grep("unable",sp.attr$SciName),] #Unable to identify; N = 1

##Remove riparian plots
sp.d <- #Remove riparian plots (N=3)
  sp.d[-which(sp.d$TopoClass=="Riparian"),] 
sp.attr <- #Reduce species attributes to only those species in upland plots
  sp.attr[which(sp.attr$SciName%in%sp.d$SciName),]

##Set affinity for non-natives to NA
sp.attr$Origin [grep("yes",sp.attr$Exotic)] <- NA

##Characterize species origin data
length(unique(sp.d$SciName)) 
#N = 188 unique taxa in Hayman plot data
length(unique(sp.attr$SciName)) 
#N = 188 unique taxa, confirmed in origin data
length(grep("_",sp.attr$SciName)) 
#N = 171 identified species, 17 genera not identified to species
length(grep("no",sp.attr$Exotic)) 
#Of 188 taxa, 166 (88%) were native
length(which(!is.na(sp.attr$Origin))) #Of 166 native taxa, 154 (93%) were classified by biogeographic affinity. 12 were not classified.
length(which(sp.attr[which( !is.na(sp.attr$Origin) & is.na(sp.attr$Expert_class)),
              "CA_match"]=="spp")) 
#Of 154 classified species, 68 (44%) had exact species matches with CA
length(which(sp.attr[which( !is.na(sp.attr$Origin) & is.na(sp.attr$Expert_class)),
                     "CA_match"]=="gen")) 
#Of 154 classified species, 64 (42%) had genus matches with CA
length(which(sp.attr[which( !is.na(sp.attr$Origin) ),
                     "Expert_class"]=="yes" ))
#Of 154 classified species, 22 (14%) were classified by us without clear matches to CA

##Add species origin data to percent cover data
index <- #Identify position matches in sp.d
  pmatch(sp.d$SciName,sp.attr$SciName,duplicates.ok = TRUE) 
sp.d[,"origin"] <- sp.attr[index,"Origin"]
sp.d[,"origin_binary"] <- sp.d[,"origin"]
sp.d[,"origin_binary"][sp.d[,"origin_binary"]!="NTm"] <- "Non-NTm"
sp.d[sp.d$SciName == "Arctostaphylos_uva-ursi","origin_binary"] <- 
  "NTm"   #Arctostaphylos uva-ursi set to NTm, per expert opinion.

####2. Group data####
sp.d.long <- #create long dataset
  gather(sp.d,key = year,value = cover, -Plot, -FireSeverity, -FireMortality,
                 -TopoClass, -Direction, -FinalCode, -SciName, -Family, -FunctionalGrp, 
                 -NativeStatus,-origin,-origin_binary)

#Add/re-label/re-classify variables in long dataset
sp.d.long$year <- as.integer(sp.d.long$year) #Convert year to number
sp.d.long$presence <- #presence in macroplot ("0" or positive values in subplot)
  ifelse(!is.na(sp.d.long$cover),1,0) 
sp.d.long$presence_sub <- #presence in subplot (positive values in subplot)
  ifelse(is.na(sp.d.long$cover) | sp.d.long$cover==0,0,1) 
sp.d.long <- #Remove unknown origins (includes exotics); 154 species remain
  sp.d.long[-which(is.na(sp.d.long$origin_binary)),] 
sp.d.long$origin_binary[sp.d.long$origin_binary=="NTm"] <- "Northern"
sp.d.long$origin_binary[sp.d.long$origin_binary=="Non-NTm"] <- "Southern"
sp.d.long$FireSeverity <- factor(sp.d.long$FireSeverity, 
                                 levels = c("Low", "Moderate", "High"))

####3. Identify colonizations/extinctions####
sp.d.long <- 
  group_by(sp.d.long,SciName,Plot) %>% #For a species in a given plot:
  mutate(colonize = #1 = species present in current year, absent previous year
           #0 = species present in current year, and in previous year
           #NA = species absent in current year, or current year is 1997
           ifelse(presence==1, ifelse(lag(presence)==0,1,0),NA),
         extinct = #1 = species absent in current year, present previous year
           #0 = species absent in current year and previous year
           #NA = species present in current year, or current year is 1997
           ifelse(presence==0, ifelse(lag(presence)==1,1,0),NA),
         fire.colonize = #Post-fire colonization, 
           #1 = species present in current year but absent in pre-fire year (1997)
           #0 = species present in current year and present pre-fire (1997)
           #NA = species absent in current year, or current year is 1997
           ifelse(presence==1 & year!=1997, ifelse(presence[1]==0,1,0),NA),
         fire.extinct = #Post-fire extinction, 
           #1 = species absent in current year but present in pre-fire year (1997)
           #0 = species absent in current year and absent in pre-fire year (1997)
           #NA = species present in current year, or current year is 1997
           ifelse(presence==0 & year!=1997, ifelse(presence[1]==1,1,0),NA),
         abs.fire.colonize = #Absolute (long-term) post-fire colonization
           #1 = species absent in 1997, colonized in current year,
           #present in 2012, current year is not 2012, 
           #and only colonized once (no double counting).
           #NA if species eventually becomes an absolute colonizer but hasn't colonized yet so "colonize" = NA
           #0 if other conditions are not met.
           ifelse(presence[1]==0 & colonize== 1 & 
                    presence[length(presence)]==1 & year!=2012 &
                    sum(colonize,na.rm=T)==1, 1, 0),
         abs.fire.extinct = #Absolute (long-term) post-fire extinction
           #1 = species present in 1997, went extinct in current year,
           #absent in 2012, current year is not 2012, 
           #and only went extinct once (no double counting).
           #NA if species eventually becomes absolutely extinct but hasn't gone extinct yet so "extinct" = NA
           #0 if other conditions are not met.
           ifelse(presence[1]==1 & extinct==1 & 
                    presence[length(presence)]==0 & year!=2012 &
                    sum(extinct,na.rm=T)==1, 1, 0)
  )

##3b: Set up table for extinctions/colonizations
which.extinct <- 
  sp.d.long[which(sp.d.long$abs.fire.extinct == 1 & sp.d.long$year==2003),] %>%
  group_by(SciName) %>%
  summarise(hs_plots = length(Plot[FireSeverity == "High"]),
            ms_plots = length(Plot[FireSeverity == "Moderate"]),
            ls_plots = length(Plot[FireSeverity == "Low"]),
            origin = unique(origin_binary)
            )
which.extinct <- which.extinct[with(which.extinct, order(origin,SciName)),]
#clipr::write_clip(which.extinct)
which.colonize <- 
  sp.d.long[which(sp.d.long$abs.fire.colonize == 1 & sp.d.long$year==2003),] %>%
  group_by(SciName) %>%
  summarise(hs_plots = length(Plot[FireSeverity == "High"]),
            ms_plots = length(Plot[FireSeverity == "Moderate"]),
            ls_plots = length(Plot[FireSeverity == "Low"]),
            origin = unique(origin_binary)
  )
which.colonize <- which.colonize[with(which.colonize, order(origin,SciName)),]
#clipr::write_clip(which.colonize)
  

####4. Set up data for plots####
p.d <- sp.d.long %>% #Data with plots separate
  group_by(Plot, FireSeverity, origin_binary, year) %>% 
  summarise(mean_cover = mean(cover, na.rm = T),sd=sd(cover, na.rm = T),
            richness = sum(presence),
            richness_sub = sum(presence_sub),
            colonizations = sum(lag(colonize==1),na.rm=T),
            extinctions = sum(lag(extinct==1),na.rm=T),
            abs_colonizations = sum(abs.fire.colonize,na.rm=T),
            abs_extinctions = sum(abs.fire.extinct,na.rm=T))

p.d.ratio <- p.d %>% #Plot data with ratios
  group_by(Plot, FireSeverity, year) %>% 
  summarise(Prop.NTM = richness[1]/sum(richness), #Northern is first factor level; previously had been subplot level
            richness.tot = sum(richness),
            richness.NTM = richness[1],
            richness.Non_NTM = richness[2])

sev.class.grp.cover <- sp.d.long %>% #Data with all plots grouped into severity class, for cover
  group_by(FireSeverity, origin_binary,year) %>% 
  summarise(mean_cover = mean(cover, na.rm = TRUE), sd = sd(cover, na.rm = T),
            se = stdErr(cover),
            n_plots = length(unique(Plot))  )

#p.d.ratio <- p.d %>% #Plot data with ratios #Deprecated, contained above.
#  group_by(Plot, FireSeverity, year) %>% 
#  summarise(Prop.NTM = richness[1]/sum(richness), #Northern is first factor level; previously had been subplot level
#            richness.tot = sum(richness),
#            richness.NTM = richness[1],
#            richness.Non_NTM = richness[2])

sev.class.grp.richness <- p.d.ratio %>% #Data with all plots grouped into severity class, for richness
  group_by(FireSeverity, year) %>% 
  summarise(#mean_richness_tot = mean(richness.tot), 
            #se_richness_tot = stdErr(richness.tot),
            mean_northern = mean(richness.NTM), 
            se_northern = stdErr(richness.NTM),
            mean_southern = mean(richness.Non_NTM), 
            se_southern = stdErr(richness.Non_NTM),
            n_plots = length(unique(Plot))  ) %>%
  gather(key , value , mean_northern, mean_southern, 
         se_northern, se_southern) %>%
  extract(key, c("question","origin"), "(.*)\\_(.*)") %>%
  spread(question,value)

sev.class.grp.ratio <- p.d.ratio %>% #Data with all plots grouped into severity class
  #Plus NTM ratio
  group_by(FireSeverity,year) %>%  
  summarise(mean_Prop.NTM = mean(Prop.NTM), sd = sd(Prop.NTM),
            se = stdErr(Prop.NTM),
            n_plots = length(unique(Plot))  
  )

####5. Plots- cover and richness####
##Cover is pretty low and not really that interesting.
F1a_LowSev <- ggplot(
  sev.class.grp.cover[!is.na(sev.class.grp.cover$origin_binary) & 
                  sev.class.grp.cover$FireSeverity == "Low" ,]) +
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
  sev.class.grp.cover[!is.na(sev.class.grp.cover$origin_binary) & 
                        sev.class.grp.cover$FireSeverity == "Moderate" ,]) +
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
  sev.class.grp.cover[!is.na(sev.class.grp.cover$origin_binary) & 
                        sev.class.grp.cover$FireSeverity == "High" ,]) +
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

#pdf(file = paste0("./Figures/MS/FigA1_",Sys.Date(),".pdf"),width=10,height=7)
grid.arrange(F1a_LowSev,F1b_ModSev,F1c_HighSev,nrow=1)
#dev.off()

##Richness
F2<-
  ggplot(sev.class.grp.richness) +
  facet_grid(. ~ FireSeverity) +
  geom_point(aes(x=year, y=mean, col = origin)) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("gray","black"))+
  geom_errorbar(aes(x=year, ymin = mean - se, ymax = mean + se, 
                    col = origin))+
  geom_line(aes(x=year, y=mean, col = origin))+
  labs(title = "burn severity", y = "mean plot richness", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.92, 0.48))

#pdf(file = paste0("./Figures/MS/Fig2_richness_",Sys.Date(),".pdf"),width=8,height=5)
F2
#grid.arrange(F2a_LowSev,F2b_ModSev,F2c_HighSev,nrow=1)
#dev.off()

library(lme4)
library(pbkrtest)
GetME_PVals=function(m){
  require(pbkrtest) #Kenward-Rodger approximation
  # get the KR-approximated degrees of freedom
  df.KR <- get_ddf_Lb(m, fixef(m))
  coefs <- data.frame(coef(summary(m)))
  # get p-values from the t-distribution using the t-values and approximated degrees of freedom
  coefs$p.KR <- round(2 * (1 - pt(abs(coefs$t.value), df.KR)),6); #coefs #Everything is significantly different from UB
  coefs$df = df.KR
  return(coefs)
}

m2.1 <- lmer(richness~origin_binary + (1|Plot), data=p.d)
GetME_PVals(m2.1)
m2.2 <- 
  lmer(richness~year + (1|Plot), #Can adjust origin to N/S and Fire Severity to L/M/H
       data=p.d[p.d$origin_binary=="Northern" & p.d$FireSeverity == "High",])
GetME_PVals(m2.2)
m2.2.outlier <- #Effect of removing outlier; only testing this for "Southern" and "High" as per reviewer comment.
  lmer(richness~year + (1|Plot), #Can adjust origin to N/S and Fire Severity to L/M/H
       data=p.d[p.d$origin_binary=="Northern" & p.d$FireSeverity == "High" & p.d$Plot!="tc0523",])
GetME_PVals(m2.2.outlier)
#mixed effects model not converging because runs out of degrees of freedom (because not looking at treatment effect, just temporal trends within given treatment type). Use regular linear models for this as above. Deprecated, this has been fixed and now using mixed-effects models as appropriate

####6. Plots and stats- ratio####
F3a <- 
  ggplot(p.d.ratio, aes(x=year,y=Prop.NTM,col=FireSeverity))+
  geom_point()+
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  geom_smooth(method="lm")+
  scale_color_manual(values = c("forestgreen","darkgoldenrod1","red2"))+
  labs(y= "proportion of flora \nwitn north-temperate affinity",
       col = "Fire severity", tag = "a") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

F3b <- 
  ggplot(sev.class.grp.ratio)+
  geom_point(aes(x=year, y=mean_Prop.NTM, col = FireSeverity)) + 
  geom_errorbar(aes(x=year, ymin = mean_Prop.NTM - se, ymax = mean_Prop.NTM + se, 
                    col = FireSeverity))+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  geom_line(aes(x=year, y=mean_Prop.NTM, col = FireSeverity))+
  scale_color_manual(values = c("forestgreen","darkgoldenrod1","red2"))+
  labs(y= "meanproportion of flora \nwitn north-temperate affinity",
       col = "Fire severity", tag = "b") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = paste0("./Figures/MS/Fig3_propNTM",Sys.Date(),".pdf"),width=5,height=7)
grid.arrange(F3a, F3b ,nrow=2)
dev.off()

#Significant differences in pre-fire ratio due to sampling year:
p.d.ratio$Actual_Year <- 1997
p.d.ratio$Actual_Year[
  which(p.d.ratio$Plot%in%c("tc0312", "tc0523", "tc1222", "tc1421", "tc2426", "tc2515"))
  ] <- 1996
summary(lm(Prop.NTM~Actual_Year, data=p.d.ratio[p.d.ratio$year==1997,]))
plot(Prop.NTM~Actual_Year, data=p.d.ratio[p.d.ratio$year==1997,])
#No significant differences in pre-fire ratio due to sampling year

#Significant trends in ratio over time?
mr.1 <- #Change between "Low", "Moderate" and "High"
  lmer(Prop.NTM ~ year + (1|Plot), data = p.d.ratio[p.d.ratio$FireSeverity=="High",])
GetME_PVals(mr.1)
#Yes for moderate and high

mr.1.outlier <- #Test effect of removing outlier; negligible. 
  lmer(Prop.NTM ~ year + (1|Plot), 
       data = p.d.ratio[p.d.ratio$FireSeverity=="High" & p.d.ratio$Plot!="tc0523",])
GetME_PVals(mr.1.outlier)



####7. Plots- colonizations/extinctions####
#pdf(file = paste0("./Figures/MS/Fig4_Colonization_Extinction",Sys.Date(),".pdf"),width=8,height=5)
ggplot(p.d,aes(col=origin_binary,fill=FireSeverity)) +
  geom_bar(aes(x=year,y=-abs_extinctions),
           position = "dodge", stat = "summary", fun.y = "mean") + 
  geom_bar(aes(x=year,y=abs_colonizations),
           position = "dodge", stat = "summary", fun.y = "mean") +
  stat_summary(aes(year,abs_colonizations),
               fun.data = mean_se, geom = "errorbar",position="dodge")+ 
  stat_summary(aes(year,-abs_extinctions),
               fun.data = mean_se, geom = "errorbar",position="dodge")+ 
  scale_fill_manual(values = c("forestgreen","darkgoldenrod1","red2"))+
  scale_colour_manual(values=c("gray","black"))+
  scale_x_continuous(breaks=c(2002:2007), limits = c(2002, 2007.5),
                     labels = c("2002 \n(fire year)", c(2003:2007)))+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  labs(title="permanent species colonizations or extinctions \nfrom pre-fire",
       y = "# extinctions                           # colonizations")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#dev.off()

pdf(file = paste0("./Figures/MS/Fig4_Colonization_Extinction",Sys.Date(),".pdf"),width=8,height=5)
ggplot(p.d,aes(fill=FireSeverity)) +
  facet_wrap(facets = vars(origin_binary)) +
  geom_bar(aes(x=year,y=-abs_extinctions),
           position = "dodge", stat = "summary", fun.y = "mean") + 
  geom_bar(aes(x=year,y=abs_colonizations),
           position = "dodge", stat = "summary", fun.y = "mean") +
  stat_summary(aes(year,abs_colonizations),
               fun.data = mean_se, geom = "errorbar",position="dodge")+ 
  stat_summary(aes(year,-abs_extinctions),
               fun.data = mean_se, geom = "errorbar",position="dodge")+ 
  scale_fill_manual(values = c("forestgreen","darkgoldenrod1","red2"))+
  scale_colour_manual(values=c("gray","black"))+
  scale_x_continuous(breaks=c(2002:2007), limits = c(2002, 2007.5),
                     labels = c("2002 \n(fire year)", c(2003:2007)))+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  labs(y = "# extinctions                            # colonizations")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9,0.18))
dev.off()



####8. Analyze and plot environmental data####
env.d$FireSeverity <- factor(env.d$FireSeverity, levels = c("low", "moderate", "high"))
env.d.pre <- env.d[env.d$Year == 1997,]
env.d.pre$Aspect <- ifelse(between(env.d.pre$Orientation,135,315),"SW","NE")
env.d.pre$prop.NTM <- c(p.d.ratio[p.d.ratio$year == 1997,"Prop.NTM"])$Prop.NTM
  
##8.1: ComparePlots sampled in 1996 vs 1997:
env.d.pre$Actual_Year <- 1997
env.d.pre$Actual_Year[grep("1996",env.d.pre$VisitDate)] <- 1996
table(env.d.pre[,c("FireSeverity","Actual_Year")])
chisq.test(table(env.d.pre[,c("FireSeverity","Actual_Year")]))
chisq.test(table(env.d.pre[,c("FireSeverity","Actual_Year")])[c(1,2),])
chisq.test(table(env.d.pre[,c("FireSeverity","Actual_Year")])[c(2,3),])
chisq.test(table(env.d.pre[,c("FireSeverity","Actual_Year")])[c(1,3),])

#No significant differences from expected distribution of actual sample years across burn severity classes.

##8.2: Compare environmental variables in 1997 (prefire)
pairwise.t.test(env.d.pre$Elevation,env.d.pre$FireSeverity)
labs_nsd <- rep("a",3) #L-M, M-H, L-H
p_elev <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = Elevation)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  #geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 2500, label = labs_nsd) +
  theme_bw() +
  labs(y = "elevation (m)") +
  theme(axis.title.x = element_blank())

chisq.test(table(env.d.pre[,c("FireSeverity","Aspect")]))
#Pairwise chisq:
chisq.test(table(env.d.pre[,c("FireSeverity","Aspect")])[c(1,2),])
chisq.test(table(env.d.pre[,c("FireSeverity","Aspect")])[c(2,3),])
chisq.test(table(env.d.pre[,c("FireSeverity","Aspect")])[c(1,3),])
p_asp <- 
  ggplot(env.d.pre, aes (x = FireSeverity, fill = Aspect)) +
  geom_bar() +
  lims(y = c(0,12)) +
  annotate(geom = "text", x = c(1,2,3), y = 11.5, label = labs_nsd) +
  theme_bw() +
  scale_fill_manual(values=c("cadetblue3","darkgoldenrod2"))+
  scale_y_continuous(breaks =  c(0,2,4,6,8,10)) +
  theme(legend.position = c(0.816,0.7), 
        legend.background = element_blank(),
        axis.title.x = element_blank())

pairwise.t.test(env.d.pre$Slope,env.d.pre$FireSeverity)
p_slp <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = Slope)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 42, label = labs_nsd) +
  theme_bw()  +
  labs(y = "slope (%)") +
  theme(axis.title.x = element_blank())

pairwise.t.test(env.d.pre$soilcov,env.d.pre$FireSeverity)
labs_lm <- c("a","b","ab") #L-M, M-H, L-H
p_soil <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = soilcov)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 42, label = labs_lm) +
  theme_bw() +
  labs(y = "soil cover (%)") +
  theme(axis.title.x = element_blank())

pairwise.t.test(env.d.pre$littduffcov,env.d.pre$FireSeverity)
labs_lm <- c("a","b","ab") #L-M, M-H, L-H
p_litt <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = littduffcov)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 80, label = labs_lm) +
  theme_bw() +
  labs(y = "litter and duff cover (%)")+
  theme(axis.title.x = element_blank())

pairwise.t.test(env.d.pre$woodcov,env.d.pre$FireSeverity)
p_wood <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = woodcov)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 20, label = labs_nsd) +
  theme_bw() +
  labs(y = "coarse woody debris cover (%)") +
  theme(axis.title.x = element_blank())

hist(env.d.pre$TPHlivetotal)
hist(log(env.d.pre$TPHlivetotal))
pairwise.t.test(log(env.d.pre$TPHlivetotal),env.d.pre$FireSeverity)
p_tph <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = TPHlivetotal)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 2100, label = labs_nsd) +
  theme_bw() +
  labs(y = expression(paste("tree density (trees ", ha^-1, ")")),
       x = "fire severity")

hist(env.d.pre$Balivetotal)
hist(log(env.d.pre$Balivetotal))
pairwise.t.test(log(env.d.pre$Balivetotal),env.d.pre$FireSeverity)
p_ba <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = Balivetotal)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 24, label = labs_nsd) +
  theme_bw() +
  labs(y = expression(paste("basal area (", m^2*ha^-1, ")")),
       x = "fire severity")


#hist(env.d.pre$CanCovLive_FVS)
pairwise.t.test(env.d.pre$CanCovLive_FVS,env.d.pre$FireSeverity)
p_cc <- 
  ggplot(env.d.pre, aes (x = FireSeverity, y = CanCovLive_FVS)) +
  geom_jitter(width = 0.1) +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 60, label = labs_nsd) +
  theme_bw() +
  labs(y = "canopy cover (%)",
       x = "fire severity")

pdf(file = paste0("./Figures/MS/FigA2_",Sys.Date(),".pdf"),width=8,height=10)
grid.arrange(p_elev,p_asp,p_slp,p_soil,p_litt,p_wood,p_tph,p_ba,p_cc, ncol=3)
dev.off()

##8.3: Compare biogeographic ratio in 1997 (prefire) as a function of canopy cover
summary(lm(prop.NTM ~ CanCovLive_FVS, data = env.d.pre))
p_cc_ntm_all <-
  ggplot(env.d.pre, aes (y = prop.NTM, x = CanCovLive_FVS)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  theme_bw() +
  labs(title = "all plots", 
       y = "proportion of flora \nwitn north-temperate affinity", 
       x = "pre-fire canopy cover (%)", tag = "a") +
  theme(plot.title = element_text(hjust = 0.5))
summary(lm(prop.NTM ~ CanCovLive_FVS, data = env.d.pre[env.d.pre$FireSeverity=="low",]))
p_cc_ntm_ls <-
  ggplot(env.d.pre[env.d.pre$FireSeverity=="low",], aes (y = prop.NTM, x = CanCovLive_FVS)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  #geom_text(aes(label = Plot, hjust = 0, vjust = 0)) +
  theme_bw()+
  labs(title = "low severity plots only", 
       y = "proportion of flora \nwitn north-temperate affinity", 
       x = "pre-fire canopy cover (%)", tag = "b") +
  theme(plot.title = element_text(hjust = 0.5))


pdf(file = paste0("./Figures/MS/FigA3_",Sys.Date(),".pdf"),width=9,height=9)
grid.arrange(p_cc_ntm_all,p_cc_ntm_ls,ncol=1)
dev.off()
