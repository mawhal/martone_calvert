# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# created 25 September 2019 -- adpated from Species_TimeSeries

# This script produces figures of the density of algae and invertebrates, saving figures as pdf


# load libraries
library(tidyverse)
library(viridis)

## read data files
# all data
ad <- read.csv( "Data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv")
# all metadata
am <- read.csv("Data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )




# taxa to remove...for now
ad <- ad[ ad$Taxon != "Black spots on Fucus", ]



## Customizations to carry through to figures

# Choose the group
dfilt <- ad %>%
  filter( motile_sessile == "sessile" ) 
  

# merge data and metadata
# get all instances of a particular taxon
# to include all quadrats us full_join, or use left_join for quads with the taxon
d <- full_join( dfilt, am )
d$Abundance[ is.na( d$Abundance) ] <- 0

# remove Meay Channel, and 2011
d <- d %>%
  filter( Site != "Meay Channel", Year != 2011 )


# make abundances numeric
sort(unique(d$Abundance))
d$Abundance <- as.numeric( d$Abundance )
# define leveles for zones
d$Zone <- factor( d$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
d$Site <- factor( d$Site, levels = c("Fifth Beach", "West Beach", "North Beach", "Meay Channel"))
# define Year factor
d$Year <- factor( d$Year, ordered= TRUE )


# because we are focused on comparing two groups rather that the individual taxa, we will sum by group
d <- d %>%
  group_by( UID, Site, Zone, Year, Quadrat, non.alga.flag ) %>%
  summarize( Abundance=sum(Abundance) )

# set West Beach as the first level
d$Site <- relevel(d$Site, ref="West Beach") 
# rename animal to sessile invertebrate
d$non.alga.flag <- as.character(d$non.alga.flag)
d$non.alga.flag[d$non.alga.flag=="Animal"] <- "Sessile\nInvertebrate"

# all sites
# time trends in different tidal zones
windows(5,6)
(ggzone <- ggplot( d, aes(x=as.numeric(as.character(Year)),y=Abundance, col=non.alga.flag)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    # geom_smooth( se=F ) +
    stat_summary( fun.data = "mean_cl_boot", size = 0.5 ) +
    stat_summary( fun.y = "mean", geom="line", size = 0.5 ) +
    # geom_point( alpha=0.4,col='slateblue' ) + ggtitle( taxon ) + 
    xlab("Year") +
    scale_x_continuous(breaks = seq(2010,2018,2) ) ) +
    scale_color_manual( values=c("slateblue4","slateblue1") ) +
    theme_bw() + theme(legend.title=element_blank())
ggsave( "R Code and Analysis/Figs/algae+invert_means.pdf", ggzone, "pdf" )

d$Year2 <- paste0("'",substr(as.character(d$Year),3,4))
d$fakeZone <- factor(d$Zone,labels=c("HIGH","MID","LOW"))
ggplot( d, aes(x=as.numeric(as.character(Year)),y=Abundance, fill=Zone)) + 
    facet_grid(non.alga.flag~Site, scales="free_y") + 
    geom_smooth( aes(group=Zone), col='black',fill="gray25",alpha=0.15, se=T, show.legend = FALSE ) +
  geom_point(aes(group=1,fill=Zone),col='black',size=0.5,alpha=0.25) +
    # stat_summary( fun.y = "mean", geom="line", size = 0.5 ) +
  stat_summary( fun.data = "mean_cl_boot", size = 0.5, shape=21 ) +
    xlab("Year") +
    scale_x_continuous( breaks = seq(2012,2018,2), labels = unique(d$Year2)[seq(1,8,2)] )  +
  # scale_color_manual( values=c("black","gray50","gray99") ) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  ylab("Total cover (%)") + xlab("Year (2000+)") +
  theme_bw() + theme(legend.title=element_blank(),
                     legend.position = c(0.01,0.48),
                     legend.direction = "horizontal",
                     legend.justification = c(0,1),
                     legend.key.size = unit(0.2, "cm"),
                     legend.spacing.x = unit(0.1, "cm"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) #+
  guides(fill=guide_legend(override.aes=list(fill=NA)))
ggsave( "R Code and Analysis/Figs/cover_algae+invert_6panel.svg", width = 6, height=4 )


## each site at a time
ds <- d %>% group_by(Site) %>% group_split()

ylimits1 <- c(0,265)
ylimits2 <- c(-5,120)

a1 <- ggplot( filter(ds[[1]],non.alga.flag=="Algae"), aes(x=as.numeric(as.character(Year)),y=Abundance, fill=Zone)) + 
  geom_smooth( aes(group=Zone), col='black',fill="gray25",alpha=0.15, se=T, show.legend = FALSE ) +
  geom_point(aes(group=1,fill=Zone),col='black',size=0.5,alpha=0.25) +
  stat_summary( fun.data = "mean_cl_boot", size = 0.5, shape=21 ) +
  scale_x_continuous( breaks = seq(2012,2018,2), labels = unique(d$Year2)[seq(1,8,2)] )  +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  ylim(ylimits1)+
  ylab("Algae") + xlab("Year (2000+)") +
  theme_bw() + theme(legend.title=element_blank(), legend.position = 'none', 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) #+
b1 <- ggplot( filter(ds[[1]],non.alga.flag=="Sessile\nInvertebrate"), aes(x=as.numeric(as.character(Year)),y=Abundance, fill=Zone)) + 
  geom_smooth( aes(group=Zone), col='black',fill="gray25",alpha=0.15, se=T, show.legend = FALSE ) +
  geom_point(aes(group=1,fill=Zone),col='black',size=0.5,alpha=0.25) +
  stat_summary( fun.data = "mean_cl_boot", size = 0.5, shape=21 ) +
  scale_x_continuous( breaks = seq(2012,2018,2), labels = unique(d$Year2)[seq(1,8,2)] )  +
  scale_y_continuous( limits = ylimits2, position = "right") +
    scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  ylab("Sessile Invertebrate") + xlab("Year (2000+)") +
  theme_bw() + theme(legend.title=element_blank(),legend.position = 'none', 
                     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
cowplot::plot_grid(a1,b1,ncol=2)
ggsave( "R Code and Analysis/Figs/cover_algae+invert_WB.svg", width = 4, height=1.5 )

a2 <- ggplot( filter(ds[[2]],non.alga.flag=="Algae"), aes(x=as.numeric(as.character(Year)),y=Abundance, fill=Zone)) + 
  geom_smooth( aes(group=Zone), col='black',fill="gray25",alpha=0.15, se=T, show.legend = FALSE ) +
  geom_point(aes(group=1,fill=Zone),col='black',size=0.5,alpha=0.25) +
  stat_summary( fun.data = "mean_cl_boot", size = 0.5, shape=21 ) +
  scale_x_continuous( breaks = seq(2012,2018,2), labels = unique(d$Year2)[seq(1,8,2)] )  +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  ylim(ylimits1)+
  ylab("Algae") + xlab("Year (2000+)") +
  theme_bw() + theme(legend.title=element_blank(), legend.position = 'none', 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) #+
b2 <- ggplot( filter(ds[[2]],non.alga.flag=="Sessile\nInvertebrate"), aes(x=as.numeric(as.character(Year)),y=Abundance, fill=Zone)) + 
  geom_smooth( aes(group=Zone), col='black',fill="gray25",alpha=0.15, se=T, show.legend = FALSE ) +
  geom_point(aes(group=1,fill=Zone),col='black',size=0.5,alpha=0.25) +
  stat_summary( fun.data = "mean_cl_boot", size = 0.5, shape=21 ) +
  scale_x_continuous( breaks = seq(2012,2018,2), labels = unique(d$Year2)[seq(1,8,2)] )  +
  scale_y_continuous( limits = ylimits2, position = "right") +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  ylab("Sessile Invertebrate") + xlab("Year (2000+)") +
  theme_bw() + theme(legend.title=element_blank(),legend.position = 'none', 
                     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
cowplot::plot_grid(a2,b2,ncol=2)
ggsave( "R Code and Analysis/Figs/cover_algae+invert_FB.svg", width = 4, height=1.5 )

a3 <- ggplot( filter(ds[[3]],non.alga.flag=="Algae"), aes(x=as.numeric(as.character(Year)),y=Abundance, fill=Zone)) + 
  geom_smooth( aes(group=Zone), col='black',fill="gray25",alpha=0.15, se=T, show.legend = FALSE ) +
  geom_point(aes(group=1,fill=Zone),col='black',size=0.5,alpha=0.25) +
  stat_summary( fun.data = "mean_cl_boot", size = 0.5, shape=21 ) +
  scale_x_continuous( breaks = seq(2012,2018,2), labels = unique(d$Year2)[seq(1,8,2)] )  +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  ylim(ylimits1)+
  ylab("Algae") + xlab("Year (2000+)") +
  theme_bw() + theme(legend.title=element_blank(), legend.position = 'none', 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) #+
b3 <- ggplot( filter(ds[[3]],non.alga.flag=="Sessile\nInvertebrate"), aes(x=as.numeric(as.character(Year)),y=Abundance, fill=Zone)) + 
  geom_smooth( aes(group=Zone), col='black',fill="gray25",alpha=0.15, se=T, show.legend = FALSE ) +
  geom_point(aes(group=1,fill=Zone),col='black',size=0.5,alpha=0.25) +
  stat_summary( fun.data = "mean_cl_boot", size = 0.5, shape=21 ) +
  scale_x_continuous( breaks = seq(2012,2018,2), labels = unique(d$Year2)[seq(1,8,2)] )  +
  scale_y_continuous( limits = ylimits2, position = "right") +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  ylab("Sessile Invertebrate") + xlab("Year (2000+)") +
  theme_bw() + theme(legend.title=element_blank(),legend.position = 'none', 
                     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
cowplot::plot_grid(a3,b3,ncol=2)
ggsave( "R Code and Analysis/Figs/cover_algae+invert_NB.svg", width = 4, height=1.5 )

# make a stacked version
# calculate means
windows(5,4)
dmean <- d %>%
  group_by( Site, Zone, Year, non.alga.flag ) %>%
  summarize( Abundance = mean(Abundance,na.rm=T) )
(ggstack <- ggplot( dmean, aes(x=as.numeric(as.character(Year)),y=Abundance, fill=non.alga.flag)) + 
    facet_grid(Site~Zone) + 
    # geom_smooth( se=F ) +
    geom_area( ) + 
    xlab("Year") +
    scale_x_continuous(breaks = seq(2010,2018,2) ) ) +
    scale_fill_manual( values=c("slateblue1","slateblue4") ) +
  theme_bw() + theme(legend.title=element_blank()) 


# subset of sites where elevation has been measured
delev <- d[ d$Site != "Meay Channel", ]
windows(10,4)
(ggheight <- ggplot( delev, aes(x=Shore_height_cm,y=Abundance)) + 
    facet_grid(Site~Year, scales = "free_y") + 
    geom_point(alpha=0.2) +  ggtitle( taxon ) + 
    geom_smooth(method="glm", method.args=list(family="quasipoisson"), 
                formula = ceiling(y) ~ poly(x,2), 
                se=FALSE, lwd=0.5) )
# ggsave( paste0("R Code and Analysis/Figs/",taxon,"_elevation_wide.pdf"), ggheight, "pdf" )

windows(4,6)
(ggheight2 <- ggplot( delev, aes(x=Shore_height_cm,y=Abundance,group=Year,col=Year )) + 
    facet_wrap(~Site,ncol=1, scales = "free_y") + 
    geom_point(alpha=0.75) +  ggtitle( taxon ) + 
    geom_smooth(method="glm", method.args=list(family="poisson"), 
                formula = ceiling(y) ~ poly(x,2), 
                se=FALSE, lwd=0.5) +
    # geom_smooth(aes(group=1)) + 
    scale_x_continuous(trans='log10') ) +
  scale_color_viridis_d( direction=-1 )

# ggsave( paste0("R Code and Analysis/Figs/",taxon,"_elevation.pdf"), ggheight2, "pdf" )

