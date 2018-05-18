

####0. Read libraries####
library(readxl) ##version 1.0.0; for read_excel()
library(tidyverse) ##version 1.2.1; for read_csv (from readr v 1.1.1)
library(gridExtra) #version 2.3; for grid.arrange()
stdErr <- function(x) sqrt(var(x, na.rm = T)/length(na.exclude(x)))

####1. Read and process data####
##Read component data
sp.d <- #Percent cover data. Warnings ok (have to do with NA's)
  read_excel("Data/Raw/tc understory species cover hayman 1997-2012 for jens.xlsx",
             sheet = "data",
             col_types = c(rep("text",2),"numeric",rep("text",7),rep("numeric",7)))
sp.attr <- #Species origin data
  read_csv("./Data/Raw/HaymanOrigins.csv")

##Remove bad species
sp.d <- sp.d[-grep("Triticosecale",sp.d$SciName),] #Sterile hybrid
sp.attr <- sp.attr[-grep("Triticosecale",sp.attr$SciName),] #Sterile hybrid
sp.d <- sp.d[-grep("unable",sp.d$SciName),] #Unable to identify
sp.attr <- sp.attr[-grep("unable",sp.attr$SciName),] #Unable to identify

##Set affinity for non-natives to NA
sp.attr$Origin [grep("yes",sp.attr$Exotic)] <- NA

##Characterize species origin data
length(unique(sp.d$SciName)) 
#N = 263 unique taxa in Hayman plot data
length(unique(sp.attr$SciName)) 
#N = 263 unique taxa, confirmed in origin data
length(grep("_",sp.attr$SciName)) 
#N = 238 identified species, 25 genera not identified to species
length(grep("no",sp.attr$Exotic)) 
#Of 263 taxa, 231 (87%) were native
length(which(!is.na(sp.attr$Origin))) #Of 231 native taxa, 218 (94%) were classified by biogeographic affinity.
length(which(sp.attr[which( !is.na(sp.attr$Origin) & is.na(sp.attr$Expert_class)),
              "CA_match"]=="spp")) 
#Of 216 classified species, 98 had exact species matches with CA
length(which(sp.attr[which( !is.na(sp.attr$Origin) & is.na(sp.attr$Expert_class)),
                     "CA_match"]=="gen")) 
#Of 216 classified species, 92 had genus matches with CA
length(which(sp.attr[which( !is.na(sp.attr$Origin) ),
                     "Expert_class"]=="yes" ))
#Of 216 classified species, 28 were classified by us without clear matches to CA

##Add species origin data to percent cover data
sp.d$SciName <- gsub('\\s+','_',sp.d$SciName ) #Replace space w/underscore for consistency
index <- pmatch(sp.d$SciName,sp.attr$SciName,duplicates.ok = TRUE) #Identify position matches in sp.d
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
sp.d.long$presence_sub <- #presence in subplot
  ifelse(is.na(sp.d.long$cover) | sp.d.long$cover==0,0,1) #(positive values in subplot)
sp.d.long <- #Remove unknown origins (includes exotics); 218 species remain
  sp.d.long[-which(is.na(sp.d.long$origin_binary)),] 
sp.d.long <- #Remove riparian plots; 160 species remain
  sp.d.long[-which(sp.d.long$TopoClass=="Riparian"),] 
sp.d.long$origin_binary[sp.d.long$origin_binary=="NTm"] <- "Northern"
sp.d.long$origin_binary[sp.d.long$origin_binary=="Non-NTm"] <- "Southern"
sp.d.long$FireSeverity <- factor(sp.d.long$FireSeverity, 
                                 levels = c("High", "Moderate", "Low"))

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

which.extinct <- 
  sp.d.long[which(sp.d.long$abs.fire.extinct == 1 & sp.d.long$year==2003),] %>%
  group_by(SciName,FireSeverity) %>%
  summarise(n_plots = length(Plot),
            origin = unique(origin)
            )
clipr::write_clip(which.extinct)

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

p.d.ratio <- p.d %>% #Plot data with ratios
  group_by(Plot, FireSeverity, year) %>% 
  summarise(Prop.NTM = richness[1]/sum(richness), #Northern is first factor level; previously had been subplot level
            richness.tot = sum(richness),
            richness.NTM = richness[1],
            richness.Non_NTM = richness[2])

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
                  sev.class.grp.cover$FireSeverity == "Low" ,]
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
  sev.class.grp.cover[!is.na(sev.class.grp.cover$origin_binary) & 
                        sev.class.grp.cover$FireSeverity == "Moderate" ,]
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
  sev.class.grp.cover[!is.na(sev.class.grp.cover$origin_binary) & 
                        sev.class.grp.cover$FireSeverity == "High" ,]
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

pdf(file = paste0("./Figures/MS/FigA1_",Sys.Date(),".pdf"),width=10,height=7)
grid.arrange(F1a_LowSev,F1b_ModSev,F1c_HighSev,nrow=1)
dev.off()

##Richness
F2a_LowSev <- ggplot(
  sev.class.grp.richness[sev.class.grp.richness$FireSeverity == "Low" ,] ) +
  geom_point(aes(x=year, y=mean, col = origin)) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("blue","darkred"))+
  geom_errorbar(aes(x=year, ymin = mean - se, ymax = mean + se, 
                    col = origin))+
  geom_line(aes(x=year, y=mean, col = origin))+
  labs(title = "Low severity", y = "mean plot richness", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none")

F2b_ModSev <- ggplot(
  sev.class.grp.richness[sev.class.grp.richness$FireSeverity == "Moderate" ,] ) +
  geom_point(aes(x=year, y=mean, col = origin)) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("blue","darkred"))+
  geom_errorbar(aes(x=year, ymin = mean - se, ymax = mean + se, 
                    col = origin))+
  geom_line(aes(x=year, y=mean, col = origin))+
  labs(title = "Moderate severity", y = "mean plot richness", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none")

F2c_HighSev <- ggplot(
  sev.class.grp.richness[sev.class.grp.richness$FireSeverity == "High" ,] ) +
  geom_point(aes(x=year, y=mean, col = origin)) + 
  geom_vline(aes(xintercept = 2002), lty = 2) +
  scale_color_manual(values = c("blue","darkred"))+
  geom_errorbar(aes(x=year, ymin = mean - se, ymax = mean + se, 
                    col = origin))+
  geom_line(aes(x=year, y=mean, col = origin))+
  labs(title = "High severity", y = "mean plot richness", 
       col = "biogeographic \naffinity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.78, 0.51))

pdf(file = paste0("./Figures/MS/Fig1_richness_",Sys.Date(),".pdf"),width=10,height=7)
grid.arrange(F2a_LowSev,F2b_ModSev,F2c_HighSev,nrow=1)
dev.off()


####6. Plots- ratio####
F3a <- ggplot(p.d.ratio, aes(x=year,y=Prop.NTM,col=FireSeverity))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title = "proportion of flora \nwitn north-temperate affinity", y= "proportion",
       col = "Fire severity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
F3b <- ggplot(sev.class.grp.ratio)+
  geom_point(aes(x=year, y=mean_Prop.NTM, col = FireSeverity)) + 
  geom_errorbar(aes(x=year, ymin = mean_Prop.NTM - se, ymax = mean_Prop.NTM + se, 
                    col = FireSeverity))+
  geom_vline(aes(xintercept = 2002), lty = 2) +
  geom_line(aes(x=year, y=mean_Prop.NTM, col = FireSeverity))+
  labs(title = "mean proportion of flora \nwitn north-temperate affinity", y= "proportion",
       col = "Fire severity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = paste0("./Figures/MS/Fig3_propNTM",Sys.Date(),".pdf"),width=5,height=7)
grid.arrange(F3a, F3b ,nrow=2)
dev.off()

#Need to make these mixed-effects models and control for plot-level random effects.
summary(lm(Prop.NTM~year, data=p.d.ratio[p.d.ratio$FireSeverity=="High",]))
summary(lm(Prop.NTM~year, data=p.d.ratio[p.d.ratio$FireSeverity=="Moderate",]))
summary(lm(Prop.NTM~year, data=p.d.ratio[p.d.ratio$FireSeverity=="Low",]))

#pdf(file = paste0("./Figures/MS/Fig4_meanpropNTM",Sys.Date(),".pdf"),width=10,height=7)
#Deprecated and folded into Fig3
#dev.off()

####7. Plots- colonizations/extinctions####
pdf(file = paste0("./Figures/MS/Fig4_Colonization_Extinction",Sys.Date(),".pdf"),width=10,height=7)
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
