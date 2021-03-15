##Functional groups through time
##This code imports the community matrix, collapses by functional groups and makes some plots
##Code written by Sam Starko
library(tidyverse)
library(Rmisc)
library(cowplot)

#Import functional groups by species
FunGroups<-read_csv("./Data/taxa/Algae_functional_groups.csv")
FunGroups$taxon<-gsub(" ",".",FunGroups$taxon)
# group seagrass with large browns
FunGroups$funct_Sep2020 <- gsub("large_brown","canopy",FunGroups$funct_Sep2020)
FunGroups$funct_Sep2020 <- gsub("seagrass","canopy",FunGroups$funct_Sep2020)

#Import community matrix
comm<-read_csv("R/output/community.csv")
#Check which species names don't match #These should all be animals and other non-algal fields
colnames(comm)[colnames(comm) %in% FunGroups$taxon == "FALSE"] 
#Remove non-algal fields from community matrix
comm2<-comm[,colnames(comm) %in% FunGroups$taxon]
# keep inverts
commuse <- comm

#Creat matrix of functional groups summed
taxon<-data.frame(taxon = colnames(commuse))
taxon.key<-left_join(taxon, FunGroups, by = "taxon")
taxon.key$funct_Sep2020[taxon.key$taxon %in% c('Barnacles','Anemone','Mytilus.sp.')] <- "animal"
taxon.key$funct_Sep2020[taxon.key$taxon %in% c('Bryozoan','Tube.worms','Tunicata/Porifera',
                                               'Pollicipes.polymerus','Pododesmus.sp.',
                                               'Hydroid')] <- "animal"
taxon.key$taxon[is.na(taxon.key$funct_Sep2020)]

# how many taxa (S for species richness) in each group
group_richness <- taxon.key %>% 
  dplyr::group_by(funct_Sep2020) %>% 
  dplyr::summarise(S=length(taxon))
colnames(commuse)<-taxon.key$funct_Sep2020
comm3<-commuse[,which(colnames(commuse) != "NA")]
comm.sum <- t(rowsum(t(comm3), group = colnames(comm3), na.rm = T)) %>% as_tibble()
comm.sum$Year<-comm$Year
comm.sum$Site<-comm$Site
comm.sum$Zone<-comm$Zone

# blade<-summarySE(data = comm.sum, measurevar = "blade", groupvars = c(as.character("Year"), "Site", "Zone"))
# colnames(blade)[5]<-"mean"
# blade$FunGroup<-"blade"
# crust<-summarySE(data = comm.sum, measurevar = "crust", groupvars = c(as.character("Year"), "Site", "Zone"))
# colnames(crust)[5]<-"mean"
# crust$FunGroup<-"crust"
# thin_turf<-summarySE(data = comm.sum, measurevar = "thin_turf", groupvars = c(as.character("Year"), "Site", "Zone"))
# colnames(thin_turf)[5]<-"mean"
# thin_turf$FunGroup<-"thin_turf"
# canopy<-summarySE(data = comm.sum, measurevar = "canopy", groupvars = c(as.character("Year"), "Site", "Zone"))
# colnames(canopy)[5]<-"mean"
# canopy$FunGroup<-"canopy"
# # seagrass<-summarySE(data = comm.sum, measurevar = "seagrass", groupvars = c(as.character("Year"), "Site", "Zone"))
# # colnames(seagrass)[5]<-"mean"
# # seagrass$FunGroup<-"seagrass"
# turf<-summarySE(data = comm.sum, measurevar = "turf", groupvars = c(as.character("Year"), "Site", "Zone"))
# colnames(turf)[5]<-"mean"
# turf$FunGroup<-"turf"

# tidy and calculate summaries
comm.tidy <- comm.sum %>%  
  pivot_longer( cols = animal:turf, names_to = "FunGroup" )
  # pivot_longer( cols = blade:turf, names_to = "FunGroup" )

d <- comm.tidy %>% 
  dplyr::group_by( Year, Site, Zone, FunGroup ) %>% 
  dplyr::summarize( mean=mean(value) )

# d<-rbind(blade, crust, thin_turf, canopy, turf) %>% as_tibble()
d$Zone<-factor(d$Zone, levels = c("LOW", "MID", "HIGH"))
d$zone_num <- as.numeric(d$Zone)
d$Site <- factor(d$Site, levels=c("North Beach","Fifth Beach","Foggy Cove"))
d$site_num <- as.numeric(d$Site)
# combine Site and Zone
# d_sitezone <- d %>% 
#   unite( "sitezone",Site,Zone, sep=" ", remove = FALSE) 
# d <- d_sitezone
# d$sitezone <- factor(d$sitezone, levels=unique(d$sitezone[order(d$Site,d$Zone)]), ordered=TRUE)

d$FunGroup<-factor(d$FunGroup, levels = rev(c("animal","turf", "thin_turf", "blade", "crust", "canopy")))
# d$FunGroup<-factor(d$FunGroup, levels = rev(c("crust", "blade", "thin_turf",  "turf", "canopy")))

# add taxon richness to fun group names
group_richness <- group_richness %>% 
  select( FunGroup=funct_Sep2020, S) %>% 
  filter( !is.na(FunGroup) )
d <- left_join( d, group_richness )
dplot <- d %>% 
  mutate( FunGroup = gsub("_"," ",FunGroup) ) %>% 
  mutate( `Functional Group` = factor(paste0(FunGroup, " (",S,")")) )
dplot$`Functional Group` <-   forcats::fct_relevel( dplot$`Functional Group`, 
                                                    "canopy (13)","blade (15)","crust (13)","thin turf (30)","turf (47)","animal (9)")

# add bare rock data
bareraw <- read_csv( "R/output/bare.csv")
bare2020 <- read_csv( "Data/Excel Files/2020_short_survey/bare_rock.csv" )
bare <- bind_rows( bareraw, bare2020 )
bare$Zone<-factor(bare$Zone, levels = c("LOW", "MID", "HIGH"))
bare$Site <- factor(bare$Site, levels=c("North Beach","Fifth Beach","Foggy Cove"))


dplot2 <- dplot %>% 
  dplyr::group_by( Year, FunGroup, `Functional Group` ) %>% 
  dplyr::summarize( mean=mean(mean) )
fun1 <- ggplot(dplot2, aes(x = Year, y = mean)) +
  geom_bar(aes(fill = `Functional Group`), position="stack", 
           stat="identity", col='black', lwd=0.25)+
  geom_smooth( data=bare, aes(x=Year,y=Abundance, group=1), 
               fill="black",col="yellow" ) +
  theme_classic()+
  scale_fill_manual(values = c("white", "darkred", "red","pink", "darkgrey", "#996633") %>% rev())+  #"darkgreen",
  theme( legend.position = 'top',
         legend.title = element_blank(),
         legend.text = element_text(size=8),
         legend.key.size = unit(0.25, "cm") ) +
  guides( fill =  guide_legend(nrow=2,byrow=T) ) +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         legend.text=element_text(size=7)) +
  ylab("Mean cover (%)")
ggsave( "R/Figs/FunGroups_time.svg", fun1, width=3,height=3.5)

ggplot(dplot, aes(x = Year, y = mean)) +
  geom_bar(aes(fill = `Functional Group`), position="stack", 
           stat="identity", col='black', lwd=0.25)+
  geom_smooth( data=bare, aes(x=Year,y=Abundance, group=1), 
               fill="black",col="yellow" ) +
  theme(plot.title = element_text(size = 5))+
  theme_cowplot()+
  facet_wrap(Site~Zone,dir = "h")+
  # facet_grid(Site~Zone)+
  scale_fill_manual(values = c("white", "darkred", "red","pink", "darkgrey", "#996633") %>% rev())+  #"darkgreen",
  # scale_fill_manual(values = c("pink", "darkgrey", "darkred", "red", "#996633") %>% rev())+  #"darkgreen", 
  scale_x_continuous(guide = guide_axis(n.dodge = 2)) +
  theme( axis.text.x = element_text(size=10) ) +
  # theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  ylab("Mean cover (%)")
ggsave( "R/Figs/FunGroups_time_zone_site.svg",width=8,height=7)


# what is the temporal coeffient of variation for each functional group? Which group varied the most?
dcv <- d %>% 
  dplyr::group_by( FunGroup ) %>% 
  dplyr::summarize( sd=sd(mean), CV = sd(mean)/mean(mean) ) %>% 
  dplyr::group_by( FunGroup ) %>% 
  dplyr::mutate( CVmean=mean(CV) )
dcv$CVrank <- forcats::fct_reorder( dcv$FunGroup, dcv$CVmean, .desc=T ) 

ggplot( data=dcv, aes( x=CVrank, y=CV)) + #facet_grid(Site~Zone) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_point( size= 4 ) 

# year-to-year differences
dsum <- d %>% 
  dplyr::group_by(FunGroup, Site, Zone) %>% 
  dplyr::summarize( change=diff(mean) ) %>%
  mutate( year = gl(7,1,315,labels = c(2013:2019) ))

# point colors based on positive or negative anomaly
dsum$sign <- ifelse(dsum$change>0,"red","blue")
ggplot( data=dsum, aes( x=year, y=change, col=FunGroup)) + facet_grid(Site~Zone) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_hline(yintercept=0) +
  xlab("Year") + ylab("Difference in percent cover from previous year") +
  geom_line(aes(col=FunGroup,group=FunGroup)) + geom_point(size= 3 )  +
  scale_color_manual(values = c("darkred", "red","pink", "darkgrey", "#996633") %>% rev())  #"darkgreen",)
ggsave( "R Code and Analysis/Figs/FunGroups_year-to-year_differences.svg",width=6.25,height=4.75)

library(broom)
dsum %>% 
  nest(data = -FunGroup) %>% 
  mutate(
    test = map(data, ~ lm(change~1, data=.x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(tidied)

dsum$year_ord <- factor(dsum$year,ordered=T)
lm1 <- lm( change~year*FunGroup, data=dsum)
anova(lm1)
library(lsmeans)
lsm1 <- lsmeans(lm1, ~year | FunGroup, cov.reduce=F)
lsmip( lsm1, FunGroup~year )
plot(lsm1, intervals = TRUE, int.adjust = "none", comparisons = TRUE)
lst1 <-  lstrends (lm1, ~ FunGroup, var = "year")
