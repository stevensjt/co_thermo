##This is the analysis code for:
##Stevens, J. T., Miller, J. E. D., Fornwalt, P. M. 2019. Fire severity and changing composition of forest understory plant communities. Journal of Vegetation Science. DOI: https://doi.org/10.1111/jvs.12796
##Code author: Jens T. Stevens; stevensjt@gmail.com
##Data DOI: https://doi.org/10.2737/RDS-2019-0031

####0. Read libraries and define functions####
library(readxl) ##version 1.0.0; for read_excel()
library(tidyverse) ##version 1.2.1; for read_csv (from readr v 1.1.1)
library(gridExtra) #version 2.3; for grid.arrange()
library(lme4) #version 1.1-19; for lmer()
library(pbkrtest) #version 0.4-7; for get_ddf_Lb()
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
GetME_PVals <- function(m){
  require(pbkrtest) #Kenward-Rodger approximation
  # get the KR-approximated degrees of freedom
  df.KR <- get_ddf_Lb(m, fixef(m))
  coefs <- data.frame(coef(summary(m)))
  # get p-values from the t-distribution using the t-values and approximated degrees of freedom
  coefs$p.KR <- round(2 * (1 - pt(abs(coefs$t.value), df.KR)),6); #coefs 
  coefs$df = df.KR
  return(coefs)
}
stdErr <- function(x) sqrt(var(x, na.rm = T)/length(na.exclude(x)))

####1. Read and process data####
##Read component data
sp.d <- #Species presence data. Warnings ok (have to do with NA's)
  read_csv("Data/PlotUnderstoryPlantComposition.csv")
sp.attr <- #Species origin data
  read_csv("./Data/UnderstoryPlantCharacteristics.csv")
env.d <- #Read other environmental data
    read_csv("Data/PlotEnvironmentalCharacteristics.csv")

##Remove bad species
sp.d <- sp.d[-grep("Triticosecale",sp.d$SciName),] #Sterile hybrid; N = 25 plots 
sp.attr <- sp.attr[-grep("Triticosecale",sp.attr$SciName),] #Sterile hybrid; N = 2 taxa
sp.d <- sp.d[-grep("unable",sp.d$SciName),] #Unable to identify; N = 18 plots
sp.attr <- sp.attr[-grep("unable",sp.attr$SciName),] #Unable to identify; N = 1 taxa

##Set affinity for non-natives to NA
sp.attr$Origin [grep("yes",sp.attr$Exotic)] <- NA

##Characterize species origin data:
length(unique(sp.d$SciName)) 
#N = 188 unique taxa in Hayman plot data
length(unique(sp.attr$SciName)) 
#N = 188 unique taxa, confirmed in origin data
length(grep("_",sp.attr$SciName)) 
#N = 171 identified species, 17 genera not identified to species
length(grep("no",sp.attr$Exotic)) 
#Of 188 taxa, 166 (88%) were native
length(which(!is.na(sp.attr$Origin))) 
#Of 166 native taxa, 154 (93%) were classified by biogeographic affinity. 12 were not classified.
length(which(sp.attr[which( !is.na(sp.attr$Origin) ),"CAMatch"]=="spp" ))-1
#Of 154 classified species, 63 (41%) had a match to a R&A classified species 
#(-1 is because we changed classification of Artostaphylos_uva_ursi per Notes column in data)
length(which(sp.attr[which( !is.na(sp.attr$Origin) ),"CAMatch"]=="gen" ))
#Of 154 classified species, 64 (42%) had a match to a R&A classified genus 
length(which(sp.attr[which( !is.na(sp.attr$Origin) ),"ExpertClass"]=="yes" ))
#Of 154 classified species, 22 (14%) were classified by us because they didn't have clear matches to CA
#Subsequent to initial submission, we realized that 5 additional species were 
#classified by us but not properly coded "ExpertClass=yes" in the original database:
#Arctostaphylos_uva-ursi, Gutierrezia_sarothrae, Physaria_vitulifera,
#Schizachyrium_scoparium, and Townsendia_grandiflora.
#Bringing the total to 27 (17%).

##Add species origin data to species presence data
index <- #Identify position matches in sp.d
  pmatch(sp.d$SciName,sp.attr$SciName,duplicates.ok = TRUE) 
sp.d[,"origin"] <- sp.attr[index,"Origin"]
sp.d[,"origin_binary"] <- sp.d[,"origin"]
sp.d[,"origin_binary"][sp.d[,"origin_binary"]!="NTm"] <- "Non-NTm"

##Clean working directory at the end of the section
rm(index) 

####2. Group data####
sp.d.long <- #create long dataset
  gather(sp.d,key = year,value = presence, -Plot, -FireSeverity, -SciName, -origin, -origin_binary)

#Add/re-label/re-classify variables in long dataset
sp.d.long$year <- as.integer(sp.d.long$year) #Convert year to number
sp.d.long <- #Remove unknown origins (includes exotics); 154 species remain
  sp.d.long[-which(is.na(sp.d.long$origin_binary)),] 
sp.d.long$origin_binary[sp.d.long$origin_binary=="NTm"] <- "cool-mesic"
sp.d.long$origin_binary[sp.d.long$origin_binary=="Non-NTm"] <- "warm-xeric"
sp.d.long$FireSeverity <- factor(sp.d.long$FireSeverity, 
                                 levels = c("Low", "Moderate", "High"),
                                 labels = c("low", "moderate", "high"))

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
           #NA if species eventually becomes an absolute colonizer but hasn't colonized yet
           #0 if other conditions are not met.
           ifelse(presence[1]==0 & colonize== 1 & 
                    presence[length(presence)]==1 & year!=2012 &
                    sum(colonize,na.rm=T)==1, 1, 0),
         abs.fire.extinct = #Absolute (long-term) post-fire extinction
           #1 = species present in 1997, went extinct in current year,
           #absent in 2012, current year is not 2012, 
           #and only went extinct once (no double counting).
           #NA if species eventually becomes absolutely extinct but hasn't gone extinct yet
           #0 if other conditions are not met.
           ifelse(presence[1]==1 & extinct==1 & 
                    presence[length(presence)]==0 & year!=2012 &
                    sum(extinct,na.rm=T)==1, 1, 0)
  )

##3b: Set up table for extinctions/colonizations
which.extinct <- 
  sp.d.long[which(sp.d.long$abs.fire.extinct == 1 & sp.d.long$year==2003),] %>%
  group_by(SciName) %>%
  summarise(hs_plots = length(Plot[FireSeverity == "high"]),
            ms_plots = length(Plot[FireSeverity == "moderate"]),
            ls_plots = length(Plot[FireSeverity == "low"]),
            origin = unique(origin_binary)
            )
which.extinct <- which.extinct[with(which.extinct, order(origin,SciName)),]

which.colonize <- 
  sp.d.long[which(sp.d.long$abs.fire.colonize == 1 & sp.d.long$year==2003),] %>%
  group_by(SciName) %>%
  summarise(hs_plots = length(Plot[FireSeverity == "high"]),
            ms_plots = length(Plot[FireSeverity == "moderate"]),
            ls_plots = length(Plot[FireSeverity == "low"]),
            origin = unique(origin_binary)
  )
which.colonize <- which.colonize[with(which.colonize, order(origin,SciName)),]

####4. Set up data for plots####
p.d <- #Data with plots separate
  #Ignore warning messages re: uninitialised column
  sp.d.long %>% #Data with plots separate
  group_by(Plot, FireSeverity, origin_binary, year) %>% 
  summarise(richness = sum(presence),
            colonizations = sum(lag(colonize==1),na.rm=T),
            extinctions = sum(lag(extinct==1),na.rm=T),
            abs_colonizations = sum(abs.fire.colonize,na.rm=T),
            abs_extinctions = sum(abs.fire.extinct,na.rm=T))
  
p.d.ratio <- p.d %>% #Ignore warning messages re: uninitialised column
  group_by(Plot, FireSeverity, year) %>% 
  summarise(Prop.NTM = richness[1]/sum(richness), #NTM = cool-mesic;
            #cool-mesic is first factor level; previously had been subplot level
            richness.tot = sum(richness),
            richness.NTM = richness[1],
            richness.Non_NTM = richness[2])

sev.class.grp.richness <- p.d.ratio %>% #Data with all plots grouped into severity class, for richness
  group_by(FireSeverity, year) %>% 
  summarise(mean_northern = mean(richness.NTM), 
            se_northern = stdErr(richness.NTM),
            mean_southern = mean(richness.Non_NTM), 
            se_southern = stdErr(richness.Non_NTM),
            n_plots = length(unique(Plot))  ) %>%
  gather(key , value , mean_northern, mean_southern, 
         se_northern, se_southern) %>%
  extract(key, c("question","origin"), "(.*)\\_(.*)") %>%
  spread(question,value)

sev.class.grp.ratio <- p.d.ratio %>% #Data with all plots grouped into severity class
  #Plus cool-mesic (formerly north-temperate; NTM) ratio
  group_by(FireSeverity,year) %>%  
  summarise(mean_Prop.NTM = mean(Prop.NTM), sd = sd(Prop.NTM),
            se = stdErr(Prop.NTM),
            n_plots = length(unique(Plot))  
  )

####5. Plots and stats- richness####
F2<-
  ggplot(sev.class.grp.richness) +
  facet_grid(. ~ FireSeverity) +
  geom_point(aes(x=year, y=mean, 
                 col = factor(origin, labels =c("cool-mesic","warm-xeric")))) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("gray","black"))+
  geom_errorbar(aes(x=year, ymin = mean - se, ymax = mean + se, 
                    col = factor(origin, labels =c("cool-mesic","warm-xeric"))))+
  geom_line(aes(x=year, y=mean, 
                col = factor(origin, labels =c("cool-mesic","warm-xeric"))))+
  labs(title = "fire severity", y = "mean plot richness", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.55, 0.48))

pdf(file = paste0("./Figures/Fig2_richness_",Sys.Date(),".pdf"),width=8,height=5)
F2
dev.off()

m2.1 <- lmer(richness~origin_binary + (1|Plot), data=p.d)
GetME_PVals(m2.1)
m2.2 <- 
  lmer(richness~year + (1|Plot), #Can adjust origin to N/S and Fire Severity to L/M/H
       data=p.d[p.d$origin_binary=="cool-mesic" & p.d$FireSeverity == "high",])
GetME_PVals(m2.2)
m2.2.outlier <- #Effect of removing outlier; only testing this for "warm-xeric" and "high" as per reviewer comment.
  lmer(richness~year + (1|Plot), #Can adjust origin to N/S and Fire Severity to L/M/H
       data=p.d[p.d$origin_binary=="cool-mesic" & p.d$FireSeverity == "high" & p.d$Plot!="tc0523",])
GetME_PVals(m2.2.outlier)

####6. Plots and stats- ratio####
F3a <- 
  ggplot(p.d.ratio, aes(x=year,y=Prop.NTM,col=FireSeverity))+
  geom_point()+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  geom_smooth(method="lm")+
  scale_color_manual(values = c("cyan3","darkgoldenrod1","red2"))+
  labs(y= "proportion of flora \nwitn cool-mesic affinity",
       col = "fire severity", tag = "a") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

F3b <- 
  ggplot(sev.class.grp.ratio)+
  geom_point(aes(x=year, y=mean_Prop.NTM, col = FireSeverity)) + 
  geom_errorbar(aes(x=year, ymin = mean_Prop.NTM - se, ymax = mean_Prop.NTM + se, 
                    col = FireSeverity))+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  geom_line(aes(x=year, y=mean_Prop.NTM, col = FireSeverity))+
  scale_color_manual(values = c("cyan3","darkgoldenrod1","red2"))+
  labs(y= "mean proportion of flora \nwitn cool-mesic affinity",
       col = "fire severity", tag = "b") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = paste0("./Figures/Fig3_propNTM",Sys.Date(),".pdf"),width=5,height=7)
grid.arrange(F3a, F3b ,nrow=2)
dev.off()

#Significant differences in pre-fire ratio due to sampling year:
p.d.ratio$Actual_Year <- 1997
p.d.ratio$Actual_Year[
  which(p.d.ratio$Plot%in%c("tc0312", "tc0523", "tc1222", "tc1421", "tc2426", "tc2515"))
  ] <- 1996
summary(lm(Prop.NTM~Actual_Year, data=p.d.ratio[p.d.ratio$year==1997,]))
#plot(Prop.NTM~Actual_Year, data=p.d.ratio[p.d.ratio$year==1997,])
#No significant differences in pre-fire ratio due to sampling year

#Significant trends in ratio over time?
mr.1 <- #Change between "low", "moderate" and "high"
  lmer(Prop.NTM ~ year + (1|Plot), data = p.d.ratio[p.d.ratio$FireSeverity=="high",])
GetME_PVals(mr.1)
#Yes for moderate and high

mr.1.outlier <- #Test effect of removing outlier; negligible (still significant). 
  lmer(Prop.NTM ~ year + (1|Plot), 
       data = p.d.ratio[p.d.ratio$FireSeverity=="high" & p.d.ratio$Plot!="tc0523",])
GetME_PVals(mr.1.outlier)



####7. Plots- colonizations/extinctions####
pdf(file = paste0("./Figures/Fig4_Colonization_Extinction",Sys.Date(),".pdf"),width=8,height=5)
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
  scale_fill_manual(values = c("cyan3","darkgoldenrod1","red2"))+
  scale_colour_manual(values=c("gray","black"))+
  scale_x_continuous(breaks=c(2002:2007), limits = c(2002, 2007.5)) +
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  labs(fill = "fire severity",
       y = "# extinctions                            # colonizations")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9,0.18))
dev.off()



####8. Analyze and plot environmental data####
env.d$FireSeverity <- factor(env.d$FireSeverity, levels = c("low", "moderate", "high"))
env.d$Actual_Year <- substrRight(env.d$VisitDate, 4)
env.d$Aspect <- ifelse(between(env.d$Orientation,135,315),"SW","NE")
env.d$prop.NTM <- c(p.d.ratio[p.d.ratio$year == 1997,"Prop.NTM"])$Prop.NTM
  
##8.1: ComparePlots sampled in 1996 vs 1997:
table(env.d[,c("FireSeverity","Actual_Year")])
chisq.test(table(env.d[,c("FireSeverity","Actual_Year")]))
chisq.test(table(env.d[,c("FireSeverity","Actual_Year")])[c(1,2),])
chisq.test(table(env.d[,c("FireSeverity","Actual_Year")])[c(2,3),])
chisq.test(table(env.d[,c("FireSeverity","Actual_Year")])[c(1,3),])
#No significant differences from expected distribution of actual sample years across burn severity classes.

##8.2: Compare environmental variables in 1997 (prefire)
pairwise.t.test(env.d$Elevation,env.d$FireSeverity)
labs_nsd <- rep("a",3) #L-M, M-H, L-H
p_elev <- 
  ggplot(env.d, aes (x = FireSeverity, y = Elevation)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 2500, label = labs_nsd) +
  theme_bw() +
  labs(y = "elevation (m)") +
  theme(axis.title.x = element_blank())

chisq.test(table(env.d[,c("FireSeverity","Aspect")]))
#Pairwise chisq:
chisq.test(table(env.d[,c("FireSeverity","Aspect")])[c(1,2),])
chisq.test(table(env.d[,c("FireSeverity","Aspect")])[c(2,3),])
chisq.test(table(env.d[,c("FireSeverity","Aspect")])[c(1,3),])
p_asp <- 
  ggplot(env.d, aes (x = FireSeverity, fill = Aspect)) +
  geom_bar() +
  lims(y = c(0,12)) +
  annotate(geom = "text", x = c(1,2,3), y = 11.5, label = labs_nsd) +
  theme_bw() +
  scale_fill_manual(values=c("cadetblue3","darkgoldenrod2"))+
  scale_y_continuous(breaks =  c(0,2,4,6,8,10)) +
  theme(legend.position = c(0.816,0.7), 
        legend.background = element_blank(),
        axis.title.x = element_blank())

pairwise.t.test(env.d$Slope,env.d$FireSeverity)
p_slp <- 
  ggplot(env.d, aes (x = FireSeverity, y = Slope)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 42, label = labs_nsd) +
  theme_bw()  +
  labs(y = "slope (%)") +
  theme(axis.title.x = element_blank())

pairwise.t.test(env.d$SoilCov,env.d$FireSeverity)
labs_lm <- c("a","b","ab") #L-M, M-H, L-H
p_soil <- 
  ggplot(env.d, aes (x = FireSeverity, y = SoilCov)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 42, label = labs_lm) +
  theme_bw() +
  labs(y = "soil cover (%)") +
  theme(axis.title.x = element_blank())

pairwise.t.test(env.d$LittDuffCov,env.d$FireSeverity)
labs_lm <- c("a","b","ab") #L-M, M-H, L-H
p_litt <- 
  ggplot(env.d, aes (x = FireSeverity, y = LittDuffCov)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 80, label = labs_lm) +
  theme_bw() +
  labs(y = "litter and duff cover (%)")+
  theme(axis.title.x = element_blank())

pairwise.t.test(env.d$WoodCov,env.d$FireSeverity)
p_wood <- 
  ggplot(env.d, aes (x = FireSeverity, y = WoodCov)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 20, label = labs_nsd) +
  theme_bw() +
  labs(y = "woody debris cover (%)") +
  theme(axis.title.x = element_blank())

#hist(env.d$TPHLiveTotal)
#hist(log(env.d$TPHLiveTotal))
pairwise.t.test(log(env.d$TPHLiveTotal),env.d$FireSeverity)
p_tph <- 
  ggplot(env.d, aes (x = FireSeverity, y = TPHLiveTotal)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 2100, label = labs_nsd) +
  theme_bw() +
  labs(y = expression(paste("tree density (trees ", ha^-1, ")")),
       x = "fire severity")

#hist(env.d$BALiveTotal)
#hist(log(env.d$BALiveTotal))
pairwise.t.test(log(env.d$BALiveTotal),env.d$FireSeverity)
p_ba <- 
  ggplot(env.d, aes (x = FireSeverity, y = BALiveTotal)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 24, label = labs_nsd) +
  theme_bw() +
  labs(y = expression(paste("basal area (", m^2*ha^-1, ")")),
       x = "fire severity")


pairwise.t.test(env.d$CanCovLiveFVS,env.d$FireSeverity)
p_cc <- 
  ggplot(env.d, aes (x = FireSeverity, y = CanCovLiveFVS)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = mean_se, geom ="errorbar", col = "blue") + 
  stat_summary(fun.y = mean, geom = "point", col = "blue") +
  annotate(geom = "text", x = c(1,2,3), y = 60, label = labs_nsd) +
  theme_bw() +
  labs(y = "canopy cover (%)",
       x = "fire severity")

#pdf(file = paste0("./Figures/FigA2_",Sys.Date(),".pdf"),width=8,height=10)
grid.arrange(p_elev,p_asp,p_slp,p_soil,p_litt,p_wood,p_tph,p_ba,p_cc, ncol=3)
#dev.off()

##8.3: Compare biogeographic ratio in 1997 (prefire) as a function of canopy cover
summary(lm(prop.NTM ~ CanCovLiveFVS, data = env.d))
p_cc_ntm_all <-
  ggplot(env.d, aes (y = prop.NTM, x = CanCovLiveFVS, col = FireSeverity)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  labs(
       y = "proportion of pre-fire flora \nwitn cool-mesic affinity", 
       x = "pre-fire canopy cover (%)", col = "fire severity") +
  scale_color_manual(values = c("cyan3","darkgoldenrod1","red2")) +
  theme(plot.title = element_text(hjust = 0.5))
summary(lm(prop.NTM ~ CanCovLiveFVS, data = env.d[env.d$FireSeverity=="low",]))

p_cc_ntm_ls <-
  ggplot(env.d[env.d$FireSeverity=="low",], aes (y = prop.NTM, x = CanCovLiveFVS)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_bw()+
  labs(title = "low severity plots only", 
       y = "proportion of pre-fire flora \nwitn cool-mesic affinity", 
       x = "pre-fire canopy cover (%)", tag = "b") +
  theme(plot.title = element_text(hjust = 0.5))


pdf(file = paste0("./Figures/Fig5_",Sys.Date(),".pdf"),width=9,height=4.5)
p_cc_ntm_all
dev.off()
