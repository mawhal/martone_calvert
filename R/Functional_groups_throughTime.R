##Functional groups through time
##This code imports the community matrix, collapses by functional groups and makes some plots
##Code written by Sam Starko and Matt Whalen
library(tidyverse)
# library(Rmisc)
library(cowplot)

#Import functional groups by species
FunGroups<-read_csv("./Data/taxa/Algae_functional_groups.csv")
FunGroups$taxon<-gsub(" ",".",FunGroups$taxon)
FunGroups <- FunGroups %>% distinct()

#Import community matrix
comm <- read_csv("R/output/community_all.csv")
#Check which species names don't match #These should all be animals and other non-algal fields
colnames(comm)[colnames(comm) %in% FunGroups$taxon == "FALSE"] 
# select taxa
comm2 <- comm %>% 
  select(Acrosiphonia:Unknown.red.blade)
commuse <- comm2

#Creat matrix of functional groups summed
taxon<-data.frame(taxon = colnames(commuse))
taxon.key<-left_join(taxon, FunGroups, by = "taxon")
taxon.key$funct_2021[taxon.key$taxon %in% c('Barnacles','Anemone','Mytilus.sp.')] <- "animal"
taxon.key$funct_2021[taxon.key$taxon %in% c('Bryozoan','Tube.worms','Tunicata/Porifera',
                                               'Pollicipes.polymerus','Pododesmus.sp.',
                                               'Hydroid','Styela.sp.')] <- "animal"
taxon.key$taxon[is.na(taxon.key$funct_2021)]
taxon.key$FG <- taxon.key$funct_2021

# how many taxa (S for species richness) in each group
group_richness <- taxon.key %>% 
  dplyr::group_by(FG) %>% 
  dplyr::summarise(S=length(taxon))
colnames(commuse)<-taxon.key$FG
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
d$Site <- factor(d$Site, levels=c("Fifth Beach","North Beach","Foggy Cove"))
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
  select( FunGroup=FG, S) %>% 
  filter( !is.na(FunGroup) )
# manually edit based on revised numbers
group_richness$S <- c(10,12,12,13,28,41)
group_richness
d <- left_join( d, group_richness )
dplot <- d %>% 
  mutate( FunGroup = gsub("_"," ",FunGroup) ) %>% 
  mutate( `Functional Group` = factor(paste0(FunGroup, " (",S,")")) )
dplot$`Functional Group` <-   forcats::fct_relevel( dplot$`Functional Group`, 
                                                    "canopy (12)","blade (12)","crust (13)","thin turf (28)","turf (41)","animal (10)")
dplot$Zone <- factor( dplot$Zone, ordered = T )
dplot$site <- as.character(dplot$Site)

# add bare rock data
bareraw <- read_csv( "R/output/bare.csv")
# bare2020 <- read_csv( "Data/Excel Files/2020_short_survey/bare_rock.csv" )
# bare <- bind_rows( bareraw, bare2020 )
bare <- bareraw
bare$Zone<-factor(bare$Zone, levels = c("LOW", "MID", "HIGH"))
bare$Site <- factor(bare$Site, levels=c("Fifth Beach","North Beach","Foggy Cove"))
bare$site <- as.character(bare$Site)


dplot2 <- dplot %>% 
  dplyr::group_by( Year, FunGroup, `Functional Group` ) %>% 
  dplyr::summarize( mean=mean(mean) )
write_csv(dplot2, "R/output/funtional_groups_annual_mean.csv")

fun1 <- ggplot(dplot2, aes(x = Year, y = mean)) +
  geom_bar(aes(fill = `Functional Group`), position="stack", 
           stat="identity", col='black', lwd=0.25)+
  geom_smooth( data=bareraw, aes(x=Year,y=Abundance, group=1), 
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
  facet_wrap(Site~Zone,dir = "h")+
  geom_bar(aes(fill = `Functional Group`), position="stack", 
           stat="identity", col='black', lwd=0.25)+
  geom_smooth( data=bare, aes(x=Year,y=Abundance, group=1), 
               fill="black",col="yellow" ) +
  theme(plot.title = element_text(size = 5))+
  theme_cowplot()+
  # facet_grid(Site~Zone)+
  scale_fill_manual(values = c("white", "darkred", "red","pink", "darkgrey", "#996633") %>% rev())+  #"darkgreen",
  # scale_fill_manual(values = c("pink", "darkgrey", "darkred", "red", "#996633") %>% rev())+  #"darkgreen", 
  scale_x_continuous(guide = guide_axis(n.dodge = 2)) +
  theme( axis.text.x = element_text(size=10) ) +
  # theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  ylab("Mean cover (%)")
ggsave( "R/Figs/FunGroups_time_zone_site.svg",width=7,height=6.5)

# make a smaller plot showing total producer cover over time for each transect
# first get total producer cover in each quadrat then average within transect
names(comm)
taxon.key$taxon[ taxon.key$FG == "animal"]
comm3 <- comm2 %>% 
  select( !(taxon.key$taxon[ taxon.key$FG == "animal"]) )
quad.prod.total <- rowSums(comm3)  
comm4 <- data.frame( comm[c(2:4,134)], quad.prod.total )
comm5 <- comm4 %>% 
  unite( "Transect", Site, Zone, sep=" ", remove = FALSE) %>% 
  group_by(Year, Site, Zone, Transect ) %>% 
  summarize( value = mean(quad.prod.total) )
comm5$Site <- factor( comm5$Site, levels = c("Foggy Cove","Fifth Beach","North Beach"))
comm5$Zone <- factor( comm5$Zone, levels = c("LOW","MID","HIGH"))
prod_empir_trend_transect <- ggplot( comm5, aes(x = Year, y = (value), group = Transect,
                                                col = Zone)) +
  # geom_path() +
  # geom_smooth( aes(group=1), method = 'lm', se = F, lwd = 1.5, col = 'slateblue', alpha=0.5 ) +
  geom_smooth( method = 'lm', se = F,   alpha=0.5, lwd = 1.5) +
  # scale_y_continuous(trans = 'log', breaks = c(20,50,100,150)) +
  theme_classic() +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         legend.position = "top",
         legend.title = element_blank(),
         legend.text = element_text(size = 7),
         legend.key.size = unit(0.5, "cm"),
         legend.key = element_rect(colour = NA, fill = NA),
         legend.box.margin=margin(-10,-10,-10,-10)) +
  coord_cartesian( ylim = c(0,180) ) +
  scale_color_manual(values = c("black","grey","grey90")) +
  guides( lty = guide_legend(ncol=2), lwd = guide_legend(ncol=3,byrow = FALSE) ) +
  ylab("Mean seaweed cover (%)") + xlab("Year")
prod_empir_trend_transect
ggsave("R/Figs/producer_cover_time_transect.svg", width = 3, height = 3)


library(broom)
initial <- comm5 %>% filter(Year==2012)
total_change <- comm5 %>% 
  nest(data = -Transect) %>% 
  mutate(
    test = map(data, ~ lm((value)~Year, data=.x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(tidied) %>% 
  filter( term == "Year")  %>%  mutate( total_change = estimate*8 )
total_change$Site <- unlist(lapply(strsplit( total_change$Transect, split = " " ),function(z) paste(z[1],z[2])))
total_change$Zone <- unlist(lapply(strsplit( total_change$Transect, split = " " ),function(z) z[3]))

windows(4,4)
plot(x=total_change$total_change,y=rep(1,9), axes=F, xlab="", ylab="", pch=1, cex=1.5)
axis(1,line = -5 )
mtext( "change in seaweed % cover", side=1, line=-3)

library(lme4)
library(lmerTest)
lmm1 <- lmer( log(value) ~ Year*Zone + (1|Transect), data = comm5 )
summary(lmm1)
anova(lmm1)
lm1 <- lm( log(value) ~ Year*Transect, data = comm5 )
summary(lm1)
anova(lm1)

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
  mutate( year = gl(7,1,labels = c(2013:2019)) )

# point colors based on positive or negative anomaly
dsum$sign <- ifelse(dsum$change>0,"red","blue")
ggplot( data=dsum, aes( x=year, y=change, col=FunGroup)) + facet_grid(Site~Zone) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_hline(yintercept=0) +
  xlab("Year") + ylab("Difference in percent cover from previous year") +
  geom_line(aes(col=FunGroup,group=FunGroup)) + geom_point(size= 3 )  +
  scale_color_manual(values = c("darkred", "red","pink", "darkgrey", "#996633","white") %>% rev())  #"darkgreen",)
ggsave( "R/Figs/FunGroups_year-to-year_differences.svg",width=6.25,height=4.75)

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

