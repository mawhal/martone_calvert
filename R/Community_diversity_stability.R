# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen (adapted in part from Sam Starko's diversity script)
# created 26 Feb 2019


# load libraries
library(tidyverse)
library(vegan)
library(psych)
library(plotrix)


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011", ]
# remove Meay Channel
am <- am[ am$Site != "Meay Channel", ]

# remove taxa that are not coutned towards subtratum cover (i.e. mobile invertebrates)
# make it easier by replacing NA values for substratum
ds <- ad
ds$motile_sessile[ is.na(ds$motile_sessile) ] <- "Substratum"
ds <- ds %>% filter( motile_sessile=="sessile", !is.na(non.alga.flag) )  %>%
  mutate( Taxon=taxon_lumped3)


d <- ds
# # split UID so we can average by transect
splits <- strsplit( as.character(d$UID), " " )
d$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))


# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which wil reduce the size of the dataset a bit
# restrict this to sessile taxa
d.simple <- d %>%
  # filter( motile_sessile=="sessile" ) %>%
  group_by( UID, transect, Taxon, non.alga.flag ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T)) %>% 
  arrange(UID)
# # calculate average abundance by transect
# d.trans <- d.simple %>%
#   group_by( transect,taxon_lumped ) %>%
#   summarize( Abundance=mean(Abundance, na.rm=T))

# spread Taxon column out into many columns filled with abundance/cover data
d.comm.all  <- d.simple %>%
  select( -non.alga.flag ) %>% 
  spread( Taxon, Abundance, fill=0 )
d.comm.algae <- d.simple %>%
  filter( non.alga.flag == "Algae" ) %>% 
  select( -non.alga.flag ) %>% 
  spread( Taxon, Abundance, fill=0 )

# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
muse  <- am
splits <- strsplit( as.character(muse$UID), " " )
muse$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))
muse <- arrange(muse,UID)

# restrict to rows selected in metadata
d.comm.all <- d.comm.all[ d.comm.all$transect %in% muse$transect, ] 
d.comm.algae <- d.comm.algae[ d.comm.algae$transect %in% muse$transect, ] 


# # there is a quadrat without any metadata, remove this from the metadata
# noalgae <- anti_join( muse, d.comm )
# noalgae$UID
# mclean <- muse[ muse$UID != noalgae$UID, ]
mclean <- muse
# define levels for zones
mclean$Zone <- factor( mclean$Zone, levels = c("LOW","MID","HIGH"), ordered = F )
# define Site order
mclean$Site <- factor( mclean$Site, levels = c("Foggy Cove", "Fifth Beach", "North Beach" ))
# define Year factor
# mclean$Year <- factor( mclean$Year, ordered= TRUE )
mclean$transect <- with( mclean, paste(Site,Zone,Year,sep = " ") )
mtrans <- mclean %>%
  group_by( transect, Site, Zone, Year ) %>%
  summarize( Shore_height_cm=mean(Shore_height_cm,na.rm=T) )

##Sort metadata and community matrix to be the same order
# d.comm.order <- d.comm[ order(match(d.comm$transect, mtrans$transect)),]
# cbind( d.comm.order$transect, mtrans$transect )

# remove UID column from community data
comm.all <- as.matrix(d.comm.all[,-c(1:2)])
comm.algae <- as.matrix(d.comm.algae[,-c(1:2)])

anti_join(d.comm.all[,1:2], mclean[,c("UID", "transect")] )
anti_join(mclean[,c("UID", "transect")] , d.comm.all[,1:2]  )
# add a line to d.comm to get missing samples

# combine in a list
comm <- list(comm.all,comm.algae)

## Steps
## ENSPIE
# the function
ENSPIE <- function(prop){
  ifelse( sum(prop,na.rm=T)>0, 1 / sum(prop^2, na.rm=T), NA ) 
} 
## Evenness as defined as Evar in Smith & Wilson 1996 Oikos
Evar <- function( x ){
  S = length( x[x>0] )
  1 - 2/pi*atan( sum((log(x[ x>0 ]) - sum(log(x[ x>0 ]))/S)^2)/S ) 
}
Evar( comm.all[1,] )
# x=comm.all[1,]
# S = length( x[x>0] )
# x[ x>0 ]
# for each quadrat, calculate richness, Shannon diversity, Simpson Diversity, and ENSPIE
divcalcs <- function( z ){
  total.cover = rowSums( z )
  pa = ifelse( z>0, 1, 0)
  richness = rowSums( pa )
  shannon = diversity( z, "shannon" )
  simpson = diversity( z, "simpson" )
  prop = z / total.cover
  enspie = apply( prop, 1, ENSPIE )
  evar = apply( z, 1, Evar )
  return( data.frame(total.cover,richness,shannon,simpson,enspie,evar) )
}
divs <- lapply( comm, divcalcs )
divsdf <- bind_rows(divs, .id = "source")
divsdf$source <- factor( divsdf$source, levels=c("1","2"), labels=c("all","algae") )
ddf <- bind_rows( list(d.comm.all[,1:2],d.comm.algae[,1:2]) )
# ddf$source <- factor( ddf$source, levels=c("1","2"), labels=c("all","algae") )
dd <- bind_cols( ddf, divsdf )

mclean <- left_join( mclean, dd )
# replace NA
mclean <- mclean %>% replace_na(list(total.cover = 0, shannon = 0, simpson = 0, enspie = 0, richness = 0))

# # splom for all quadrat summaries
# psych::pairs.panels( mclean %>% ungroup() %>% select(total.cover,richness,shannon,simpson,enspie),
#               scale=F, ellipses = FALSE )


# # look at patterns over time
# ggplot( mclean, aes(y=richness,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
# ggplot( mclean, aes(y=shannon,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
# ggplot( mclean, aes(y=simpson,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
# ggplot( mclean, aes(y=enspie,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()




# make mclean longer and include both richness and ENSPIE in the same figure
mlong <- mclean %>%
  select( transect, Site, Zone, Year, source, Shore_height_cm, enspie, richness, evar ) %>%
  group_by( transect, Site, Zone, Year, source, Shore_height_cm ) %>%
  gather( Measure, species, -transect, -Site, -Zone, -Year, -Shore_height_cm, -source )

mlong$Measure[mlong$Measure=="enspie"] <- "Effective"
mlong$Measure[mlong$Measure=="richness"] <- "Total"

# read bar rock data

# consider the range of variation in estimates over time
divvar <- mclean %>%
  group_by( source, Site, Zone ) %>%
  summarise( meana=mean(total.cover), ea=sd(total.cover),
             meanr=mean(richness), er=sd(richness),
             meand=mean(shannon), ed=sd(shannon),
             means=mean(simpson), es=sd(simpson),
             meane=mean(enspie), ee=sd(enspie),
             meanv=mean(evar), ev=sd(evar),
             elev = mean(Shore_height_cm,na.rm=T) ) %>%
  mutate( cva = ea/meana, cvr = er/meanr, cvd = ed/meand, cvs = es/means, cve = ee/meane )

ggplot( divvar, aes(x=meanr,y=cva,col=source) ) + geom_point() +
  ylab( "CV( % cover )") + xlab("Mean transect species richness")
ggplot( divvar, aes(x=elev,y=cva,size=meanr,col=Site) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean transect elevation (cm)")

# collapse and plot all together
divplot <- divvar %>% 
  group_by(Site, Zone, source, cva) %>% 
  gather(key="metric",value="mean", meanr, meand, means, meane )

ggplot( data=divplot, aes(x=mean, y=1/cva)) + facet_grid(source~metric, scales="free") + 
  geom_smooth(method='lm',se=T) +
  geom_smooth(aes(group=Zone,lty=Zone),method='lm',se=F, col='black') +
  geom_point(aes(col=Site,size=Zone))

# compare variation by transect
divplot2 <- divvar %>% 
  group_by(Site, Zone, source, cva) %>% 
  gather(key="metric",value="sd", er, ed, es, ee )
ggplot( data=divplot2,aes(y=sd,x=Zone,col=Site)) +  facet_grid(source~metric, scales="free") + 
  geom_point()

## above calculations include variation across and within times, 
# exclude this variation within time becuase different quads sampled over time?
yearly <- mclean %>% 
  group_by( Site, Zone, Year, source ) %>% 
  summarize( meana=mean(total.cover), 
             meanr=mean(richness),
             meand=mean(shannon),
             means=mean(simpson),
             meane=mean(enspie),
             meanv=mean(evar),
             elev = mean(Shore_height_cm,na.rm=T) )

# summarize over time
divvar2 <- yearly %>%
  group_by( Site, Zone, source ) %>%
  summarise( gmeana=mean(meana), ea=sd(meana),
             gmeanr=mean(meanr), er=sd(meanr),
             gmeand=mean(meand), ed=sd(meand),
             gmeans=mean(means), es=sd(means),
             gmeane=mean(meane), ee=sd(meane),
             gmeanv=mean(meanv), ev=sd(meanv),
             gmeanelev=mean(elev), eelev=sd(elev)) %>%
  mutate( cva = ea/gmeana, cvr = er/gmeanr, cvd = ed/gmeand, cvs = es/gmeans, cve = ee/gmeane,
          stability=gmeana/ea )

# initial diversities
inits <- yearly %>% filter( Year==2012 )
divvar2 <- left_join(divvar2, inits)
# divvar2 <- filter( divvar2, source=="all" )

# collapse and plot all together
divplot2 <- divvar2 %>% 
  group_by(Site, Zone, source, cva) %>% 
  gather(key="metric",value="mean", gmeanr, gmeand, gmeans, gmeane, gmeanv, 
         meanr, meand, means, meane, meanv )

ggplot( data=divplot2, aes(x=mean, y=1/cva)) + facet_grid(source~metric, scales="free") + 
  geom_smooth(method='lm',se=T) +
  geom_smooth(aes(group=Zone,lty=Zone),method='lm',se=F, col='black') +
  geom_point(aes(col=Site,size=Zone))

# compare variation by transect
divplot2 <- divvar2 %>% 
  group_by(Site, Zone, source, cva) %>% 
  gather(key="metric",value="sd", er, ed, es, ee, ev )
ggplot( data=divplot2,aes(y=sd,x=Zone,col=Site)) +  
  facet_grid(metric~source, scales="free") + 
  geom_point()
divplot2$Zone <-  factor( divplot2$Zone, ordered=F, levels=c("MID","LOW","HIGH"))
library(tidyr)
library(purrr)
library(broom)
divplot2 %>%
  ungroup() %>% 
  filter(metric %in% c("er","ee")) %>% 
  nest(-source,-metric) %>% 
  mutate(
    fit = map(data, ~ (lm(sd ~ Zone, data = .x))),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)
divplot2 %>%
  ungroup() %>%
  filter(metric %in% c("er","ee")) %>% 
  nest(-metric,-source) %>% 
  mutate(
    fit = map(data, ~ anova(lm(sd ~ Zone, data = .x))),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)

with( divvar2, cor.test(gmeanr, stability) )
with( divvar2, cor.test(meanr, stability) )
divvar2 %>% group_by( Site, Zone) %>% summarize(diff=diff(stability))
summary(lm( log(stability,base=2) ~ gmeanr, data=filter(divvar2,source=="algae")))
mstab <- lm( log(stability,base=2) ~ gmeanr, data=filter(divvar2,source=="algae"))
summary(mstab)$adj.r.squared
divvar2_algae <- filter(divvar2,source=="algae")
divvar2_algae$Zone3
mzone3 <- lm( log(stability,base=2) ~ Zone3, data = divvar2_algae )
summary(mzone3)
stab <- ggplot( data=filter(divvar2,source=="algae"), aes(x=gmeanr, y=stability)) + #x=gmeanr
  # facet_wrap(~source) +
  # geom_smooth(method='glm',se=T,method.args=list(family="poisson"), col = "black", lwd = 0.5) +
  geom_smooth(method='lm',se=T, col='black', lwd=0.5) +
  geom_point(aes(fill=Zone,shape=Site),size=3) +
  ylab( expression(paste("Seaweed cover stability (",mu,"/",sigma,")")) ) + 
  xlab("Mean species richness") +
  annotate(geom = 'text', label = bquote( R^2 == .(round(summary(mstab)$adj.r.squared,2))), 
           x = max(divvar2_algae$gmeanr), y = 1.1, hjust = 1, vjust = 0, size = 5) +
  # annotate("text", label = "top", 
  #          x = 0.5*(min(mpg$hwy) + max(mpg$hwy)), y = max(mpg$cty), vjust = 1) +
  # annotate("text", label = "bottom", 
  #          x = 0.5*(min(mpg$hwy) + max(mpg$hwy)), y = min(mpg$cty), vjust = 0) +
  # annotate("text", label = "right", 
  #          x =  max(mpg$hwy), y = 0.5*(min(mpg$cty) + max(mpg$cty)), hjust = 1) +
  # annotate("text", label = "left", 
  #          x =  min(mpg$hwy), y = 0.5*(min(mpg$cty) + max(mpg$cty)), hjust = 0) +
  scale_shape_manual( values=c(21,22,24) ) +
  scale_fill_manual( values=c("black","gray50","whitesmoke"), guide=FALSE ) +
  scale_linetype_discrete(  guide=FALSE ) +
  scale_y_continuous(trans="log2") +
  coord_cartesian(ylim = c(1,16)) +
  # guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_classic() + theme( legend.position = c(0.01,.99), 
                           legend.justification = c(0,1),
                           legend.background = element_blank())

stab
ggsave( "R Code and Analysis/Figs/stability_richness_algae.svg",width = 3, height=4 )
stab2 <- ggplot( data=filter(divvar2,source=="algae"), aes(x=gmeanr, y=stability)) + #x=gmeanr
  # facet_wrap(~source) +
  geom_smooth(method='glm',method.args=list(family="poisson"),
              se=T, col = "black", lwd = 0.5) +
  # geom_smooth(aes(group=Zone,lty=Zone),method='lm',se=F, col='black') +
  geom_point(aes(fill=Site,shape=Zone),size=3) +
  ylab( expression(paste("Seaweed cover stability (",mu,"/",sigma,")")) ) + 
  xlab("Seaweed species richness") +
  scale_shape_manual( values=c(21,22,24), guide = F ) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  scale_linetype_discrete(  guide=FALSE ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_classic() + theme( legend.position = c(0.01,.99), legend.justification = c(0,1))
stab2

# model stability by zone
library(lme4)
library(lmerTest)
divvar2$Zone2 <- relevel( factor(divvar2$Zone, ordered=F), ref="LOW" )
divvar2$Zone4 <- relevel( factor(divvar2$Zone, ordered=F), ref="HIGH" )
divvar2$Zone3 <- as.character(divvar2$Zone2)
divvar2$Zone3[ divvar2$Zone3=="MID"] <- "HIGH"
divvar2$Zone3 <- factor(divvar2$Zone3)
# contrasts(divvar2$Zone2) <- matrix( c(1,-0.5,-0.5,1,0,-1),ncol=2)
anova(lm( stability ~ Zone, data=filter(divvar2,source=="all") ))
anova(lm( stability ~ Zone, data=filter(divvar2,source=="algae") ))
anova(lm( stability ~ Zone3, data=filter(divvar2,source=="algae") ))
summary(lm( stability ~ Zone4, data=filter(divvar2,source=="algae") ))
summary(lm( stability ~ Zone2, data=filter(divvar2,source=="all") ))
summary(lm( stability ~ Zone2, data=filter(divvar2,source=="algae")))
summary(lm( stability ~ gmeanr, data=filter(divvar2,source=="all") ))
summary(lm( stability ~ gmeanr, data=filter(divvar2,source=="algae")))
summary(lm( stability ~ meanr, data=filter(divvar2,source=="all") ))
summary(lm( stability ~ meanr, data=filter(divvar2,source=="algae")))


a <- ggplot( divvar2, aes(x=gmeanr,y=1/cva) ) + geom_point() + facet_wrap(~source) +
  ylab( "stability") + xlab("Mean species richness") + geom_smooth(method='lm') 
b <- ggplot( divvar2, aes(x=gmeand,y=1/cva) ) + geom_point() + facet_wrap(~source) +
  ylab( "stability") + xlab("Mean Shannon diversity") + geom_smooth(method='lm')
c <- ggplot( divvar2, aes(x=gmeans,y=1/cva) ) + geom_point() + facet_wrap(~source) +
  ylab( "stability") + xlab("Mean Simpson diversity") + geom_smooth(method='lm')
d <- ggplot( divvar2, aes(x=gmeane,y=1/cva) ) + geom_point() + facet_wrap(~source) +
  ylab( "stability") + xlab("Mean effective # species") + geom_smooth(method='lm') 
cowplot::plot_grid( a,b,c,d, ncol=1 )


# model stability by group (with or without inverts)
boxplot( stability~source, divvar2 )
t.test( stability~source, divvar2 )



## Species synchrony
# repeat calculations for each taxon
# we need to know if species abundance was zero, not just non-zero abundance
d.long.all <- d.comm.all %>% 
  gather( "taxon","Abundance",-UID,-transect)
yearly.taxon.all <- d.long.all %>% 
  group_by(transect, taxon) %>% 
  summarize(meana = mean(Abundance) )
d.long.algae <- d.comm.algae %>% 
  gather( "taxon","Abundance",-UID,-transect)
yearly.taxon.algae <- d.long.algae %>% 
  group_by(transect, taxon) %>% 
  summarize(meana = mean(Abundance) )
yearly.taxon <- bind_rows( yearly.taxon.all, yearly.taxon.algae, .id="source")
yearly.taxon$source <- factor( yearly.taxon$source, levels=c("1","2"), labels=c("all","algae") )

# 16 Sep 2020 - change synchrony to "log var ratio" as in Leps et al 2018 and Valencia et al 2020
# numerator is the variance in total cover over time
numer <- divvar2 %>% 
  mutate( Zone=factor(Zone,levels=c("LOW","MID","HIGH"),ordered=F)) %>% 
  select(Site,Zone,source,ea) %>% 
  group_by(Site, Zone, source) %>% 
  mutate( VT = ea^2 ) 

# denominator of the calculation for synchrony is the sum of the individual taxon variances
denom <- yearly.taxon  %>% 
  separate(transect, into=c("Site","blah","Zone","Year")) %>% 
  unite( col="Site" , Site, blah, sep = " " ) %>% 
  mutate( Zone=factor(Zone,levels=c("LOW","MID","HIGH"))) %>% 
  group_by( Site, Zone, source, taxon ) %>% 
  summarize( eai = sd(meana, na.rm=T)^2 ) %>%  # variance instead of SD
  mutate( eai = ifelse( eai>0,eai,NA) ) %>%
  group_by( Site, Zone ) %>% 
  mutate( Eeai=sum(eai, na.rm=T) ) %>% 
  mutate( Evi = Eeai ) # square operator now included with variance calculation


synch <- left_join( numer, denom )
synch <- synch %>% 
  mutate( logV = log(VT/Evi) )
summary( synch$logV )
summary( synch$eai )
synch$invert <- factor(synch$source, levels="algae","all")
tax.invert <- d.simple %>% ungroup() %>% select(taxon=Taxon,non.alga.flag) %>% distinct()
synch <- left_join( synch, tax.invert )
synch$source <- factor( synch$source, levels=c("algae","all") )
synch$non.alga.flag <- factor( synch$non.alga.flag, levels=c("Algae","Animal"), labels=c("algae","all") )
synch$non.alga.flag <- factor( synch$non.alga.flag, levels=c("algae","all") )
synch$eai2 <-  sqrt(synch$eai)
synch$eai2[ synch$source=="algae" & synch$non.alga.flag=="algae" ] <- NA

a <- ggplot( synch, aes(x=Zone,y=eai2, col=source, fill=source)) + facet_wrap(~Site) +
  geom_point( aes(y=ea), size=3, shape=21, col='black' ) +
  geom_point( aes(shape=non.alga.flag,size=non.alga.flag,alpha=non.alga.flag), col='black', fill="grey")+#, position = position_dodge(width=0.2) ) +
  ylab("Cover SD") +
  scale_color_manual( values=c("black","grey"), guide=F ) +
  scale_fill_manual( values=c("black","grey"), guide=F ) +
  scale_size_manual( values=(c(2,1)), guide=F ) +
  scale_shape_manual( values=(c(1,21)), guide=F ) +
  scale_alpha_manual( values=(c(0.1,1)), guide=F ) +
  theme( strip.background = element_blank(),
         strip.text.x = element_blank() )
b <- ggplot( synch, aes(x=Zone,y=logV,fill=source)) + facet_wrap(~Site) +
  geom_point( alpha=1, shape=21, size=3 ) +
  ylab( expression(paste("Synchrony (",logV,")")) ) +
  scale_fill_manual( values=c("black","grey"), guide=F ) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
cowplot::plot_grid( b,a, ncol=1, rel_heights = c(1,1.5), align = "v" )
ggsave( "R Code and Analysis/Figs/synchrony_transect.svg", width=6, height=4 )

# plot synchrony versus stability and richness
synchrony <- synch %>% 
  filter( !is.na(eai) ) %>% 
  group_by( Site, Zone, source, logV ) %>%
  summarize( richness = length(eai) )
dsynch <- left_join( divvar2, synchrony )
dsynch$source <- factor( dsynch$source, levels=c("all","algae"))
dsynch$Site <- factor( dsynch$Site, levels = c("Foggy Cove", "Fifth Beach", "North Beach" ))

mrich <- lm(logV ~ gmeanr, data = filter(dsynch,source == "algae"))
synchrich <- ggplot( filter(dsynch,source == "algae"), aes(x=gmeanr,y=logV,shape = Site, fill = Zone)) + 
  # facet_wrap(~source) +
  geom_smooth( aes(group=1), method="lm", se = T, col='black', lwd=0.5, show.legend = FALSE )+ 
  geom_point(size=3) +
  # geom_smooth( aes(group=1), method="glm", method.args=list(family=quasibinomial))+ geom_point(size=3) +
  scale_shape_manual( values=c(21,22,24), guide=F ) +
  scale_fill_manual( values=c("black","gray50","whitesmoke"), guide=FALSE ) +
  ylab( expression(paste("Seaweed pop. synchrony (",logV,")")) ) + 
  xlab( "Seaweed species richness" ) +
  annotate(geom = 'text', label = bquote( R^2 == .(round(summary(mrich)$adj.r.squared,2))), 
           x = min(divvar2_algae$gmeanr), y = -2.3, hjust = 0, vjust = 0, size = 5) +
  guides(fill=guide_legend("Site",override.aes = list(shape = 21)) ) +
  theme_classic() + theme( legend.position = "none" )
synchrich

msynch <- lm(log(stability,base=2) ~ logV, data = filter(dsynch,source == "algae"))
stabsynch <- ggplot( filter(dsynch,source == "algae"), aes(x=logV,y=(stability),shape = Site, fill = Zone)) + 
  geom_smooth( aes(group=1), method="lm", se=T, col='black', lwd=0.5, show.legend = FALSE) + 
  geom_point(size=3) +
  scale_shape_manual( values=c(21,22,24), guide=F ) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  # guides( fill = guide_legend("Site", override.aes = list(shape = 21)) ) + #,override.aes = list(shape = 21)
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab( expression(paste("Seaweed pop. synchrony (",logV,")")) ) + 
  ylab( expression(paste("Seaweed cover stability (",mu,"/",sigma,")"))  ) +
  annotate(geom = 'text', label = bquote( R^2 == .(round(summary(msynch)$adj.r.squared,2))), 
           x = min(dsynch$logV), y = 1.1, hjust = 0, vjust = 0, size = 5) +
  scale_y_continuous(trans="log2") +
  coord_cartesian( ylim = c(1,16)) +
  theme_classic() +
  theme( legend.position = c(0.99,0.99),legend.justification = c(1,1),
         legend.background = element_blank() )
stabsynch


cowplot::plot_grid( stabsynch, synchrich, stab, 
                    nrow = 1, align = "hv", 
                    labels = "AUTO", vjust = 1.01 )
ggsave( "R/Figs/stability_synchrony_richness.svg", width = 9, height = 3 )

ggplot( dsynch, aes(x=logV,y=(stability),shape=Zone, fill=Site)) + facet_wrap(~source) +
  # geom_smooth( aes(group=1), method="glm",method.args=list(family=quasipoisson),se=F) + 
  geom_smooth( aes(group=1), method="lm", se=F) + 
  geom_point() +
  scale_shape_manual( values=values=c(21,22,24) ) +
  scale_fill_manual( values=c("black","gray50","whitesmoke"), guide=F ) +
  guides( fill=guide_legend("Site") ) +  #,override.aes = list(shape = 21)
  xlab( expression(paste("Species synchrony (",logV,")")) ) + ylab( expression(paste("Algal cover stability (",mu,"/",sigma,")"))  ) +
  scale_y_continuous(trans="log2") +
  theme_classic() #+ theme( legend.position = c(0.99,0.99), legend.justification = c(1,1) ) 
ggsave( "R Code and Analysis/Figs/stability~synchrony.svg", width=6, height=3 )

#



###
# Richness, Evenness, Stability, Synchrony
responses <- dsynch %>% 
  ungroup() %>% 
  select( source, gmeanelev,gmeanr, gmeanv, logV, stability )

responses %>% 
  filter( source=="all" ) %>% 
  mutate( stability=log(stability), gmeanr=log(gmeanr), gmeanv=log(gmeanv) ) %>% 
  select( -source ) %>% 
  pairs.panels( )


responses %>%
  nest(-source) %>% 
  mutate(
    fitstab = map(data, ~ lm( log(stability) ~ log(gmeanr)+gmeanelev+gmeanv+logV, data = .x)),
    tidied = map(fitstab, tidy)
  ) %>% 
  unnest(tidied)

responses %>%
  nest(-source) %>% 
  mutate(
    fitstab = map(data, ~ lm( logV ~ log(gmeanr), data = .x)),
    tidied = map(fitstab, tidy)
  ) %>% 
  unnest(tidied)
###







## resistance versus resilience
## compare good years and bad years
initial <- 2012:2013
maybe_normal <- 2018
heatwave <- 2014:2017
maybe_hot <- 2019

yearly$event <- "heatwave"
yearly$event[yearly$Year %in% c(initial)] <- "normal"
yearly$event[yearly$Year %in% maybe_hot] <- "heatwave2"

  
# function to calculate resistance and resilience
omega <- function( z, window=1 ){
  Yn  = mean(z$meana[z$event=="normal"])
  Ye  = mean(z$meana[z$event=="heatwave"][1:window])
  omega = Yn / abs(Ye-Yn) 
  return(omega)
}

delta <- function( z,lag=1 ){
  Yn  = mean(z$meana[z$Year %in% initial ])
  Ye  = mean(z$meana[z$event=="heatwave"])
  Ye1 = mean(z$meana[z$Year==c(2018,2019)[1:lag]])
  delta = abs( (Ye-Yn) / (Ye1-Yn) )
  return( delta )
}


resistance1 <- by( yearly, list( factor(yearly$Zone),factor(yearly$Site),factor(yearly$source)), omega, window=1 )
resistance2 <- by( yearly, list( factor(yearly$Zone),factor(yearly$Site),factor(yearly$source)), omega, window=2  )
resistance3 <- by( yearly, list( factor(yearly$Zone),factor(yearly$Site),factor(yearly$source)), omega, window=3  )
resistance4 <- by( yearly, list( factor(yearly$Zone),factor(yearly$Site),factor(yearly$source)), omega, window=4  )

resilience1 <- by( yearly, list( factor(yearly$Zone),factor(yearly$Site),factor(yearly$source)), delta, lag=1 )
resilience2 <- by( yearly, list( factor(yearly$Zone),factor(yearly$Site),factor(yearly$source)), delta, lag=2 )

ress <- bind_cols( yearly %>% ungroup() %>% select(Site,Zone,source) %>% distinct() %>% arrange(source,Site,Zone),
           data.frame( O1=c(resistance1), O2=c(resistance2), 
            O3=c(resistance3),O4=c(resistance4), D1=c(resilience1), D2=c(resilience2) ) )


psych::pairs.panels( ress[,-c(1:2)])
psych::pairs.panels( log(ress[,-c(1:2)]) )


ress_long <- ress %>%
  gather( "measure","value", -Site, -Zone, -source )

ress_long <- left_join(ress_long,dsynch)
# ress_long <- filter(ress_long, source=="all")

# ggplot( ress_long, aes(x=logV,y=(value),fill=Site)) +
ggplot( filter(ress_long,measure %in% c("D1","O4")), aes(x=logV,y=(value),fill=Site)) +
  facet_grid(measure~source,scales="free_y") + geom_point(size=3, alpha=1, pch=21) +
  scale_fill_manual(values=c("black","gray50","whitesmoke") ) +
  geom_smooth(aes(group=1), method='lm') +
  geom_smooth( method='lm', se=F, col='black') +
  scale_y_continuous(trans="log2") +
  theme_bw()
ggsave( "R Code and Analysis/Figs/resist_resil_logV.svg", width=6, height=5 )

# ggplot( ress_long, aes(x=log(gmeanr),y=(value),fill=Site)) +
ggplot( filter(ress_long,measure %in% c("D1","O4")), aes(x=log(gmeanr),y=(value),fill=Site)) +
  facet_grid(measure~source,scales="free_y") + geom_point(size=3, alpha=1, pch=21) +
  scale_fill_manual(values=c("black","gray50","whitesmoke") ) +
  geom_smooth(aes(group=1), method='lm') +
  geom_smooth( method='lm', se=F, col='black') +
  scale_y_continuous(trans="log2") +
  theme_bw()
ggsave( "R Code and Analysis/Figs/resist_resil_rich.svg", width=6, height=5 )

resist <- ggplot( filter(ress_long,measure %in% c("O1","O2","O3","O4")), 
                  aes(x=as.numeric(as.factor(measure)),y=(value),fill=Zone,shape=Site)) +
  geom_smooth( aes(group=Zone, lty=Zone), method='glm', se=F, col="black", method.args = list(family="quasipoisson")  ) +
  geom_point(size=3) + 
  facet_wrap(~source)+
  # geom_smooth( method='glm', method.args=list(family="quasipoisson")) +
  scale_y_continuous(trans="log2") +
  scale_shape_manual( values=21:23, guide=F ) +
  scale_fill_manual(values=c("black","gray50","whitesmoke") ) +
  # scale_linetype_manual(guide=F ) +
  # guides(fill=guide_legend("Zone",override.aes = list(shape = 21))) +
  xlab("Year of heatwave") + ylab(expression(paste("Resistance (",Omega,")"))) +
  theme_classic() + theme( legend.position = "none" )
resist.df <- filter(ress_long,measure %in% c("O1","O2","O3","O4"))
resist.df$zone <- factor(resist.df$Zone, ordered=F, levels=c("MID","LOW","HIGH") )
resist.df$year <- scale( as.numeric( factor(resist.df$measure, ordered=F) ) )
resist.lm <- filter(resist.df, source=="algae")
anova( lm( log(value) ~ year*zone,  data = resist.lm ) )
summary( lm( log(value) ~ year*zone,  data = resist.lm ) )
plot( lm( log(value) ~ year*zone,  data = resist.lm ) )
anova( lm( log(value) ~ gmeanr,  data = resist.lm ) )
summary( lm( log(value) ~ gmeanr,  data = resist.lm ) )
# summary(glm(value~as.numeric(as.factor(measure))+factor(Zone,ordered=F, levels=c("MID","LOW","HIGH")), 
#             data=filter(ress_long,measure %in% c("O1","O2","O3","O4"), source=="algae"),
#             family="quasipoisson"))
# summary(glm(value~factor(Zone,ordered=F, levels=c("MID","LOW","HIGH"))*source, 
#             data=filter(ress_long,measure %in% c("O1","O2","O3","O4")),
#             family="quasipoisson"))
# summary(glm(value~as.numeric(as.factor(measure))*source, 
#             data=filter(ress_long,measure %in% c("O1","O2","O3","O4")),
#             family="quasipoisson"))
# summary(glm(value~as.numeric(as.factor(measure)), 
#             data=filter(ress_long,measure %in% c("O1","O2","O3","O4"), source=="algae"),
#             family="poisson"))
resil <-  ggplot( filter(ress_long,measure %in% c("D1")), aes(x=gmeanr,y=(value))) +
  # geom_smooth(aes(group=1), method='lm') +
  geom_smooth(aes(group=1), method='glm',method.args=list(family='quasipoisson')) +
  scale_y_continuous(trans="log2") +
  geom_point( aes(shape=Site,fill=Zone), size=3 ) + 
  facet_wrap(~source)+
  xlab("Mean species richness") + ylab(expression(paste("Resilience (",Delta,")"))) +
  scale_shape_manual( values=21:23 ) +
  scale_fill_manual(values=c("black","gray50","whitesmoke") ) +
  theme_classic() + theme( legend.position = "none" ) #theme( legend.position = c(0.01,0.99),legend.justification = c(0,1)) 
resil.df <- filter(ress_long,measure %in% c("D1"))
resil.df$zone <- factor(resil.df$Zone, ordered=F, levels=c("MID","LOW","HIGH") )
by(resil.df,resil.df$source, function(z) with( z, cor.test( value,gmeanr,method = 'spearman' )))
anova( glm( value ~ gmeanr,  data = resil.df, family="quasipoisson" ) )
summary( glm( value ~ gmeanr,  data = resil.df, family="quasipoisson" ) )
plot( glm( value ~ gmeanr*zone,  data = resil.df, family="quasipoisson" ) )

cowplot::plot_grid(resist,resil,ncol=1,rel_widths = c(1,1))
ggsave( "R Code and Analysis/Figs/stability_resist_resil.svg", width=6, height=5 )
summary(lm(value~gmeanr, data=filter(ress_long,measure %in% c("D2")) ))
summary(glm(value~gmeanr, 
            data=filter(ress_long,measure %in% c("D2")),
            family="quasipoisson"))

stabresil <- cowplot::plot_grid( stab, resil, ncol=2, labels = "AUTO" )
# cowplot::plot_grid( stabresil, resist, ncol=1 )
resistsynch <- cowplot::plot_grid( resist, synchplot, ncol=2, labels = c("C","D") )
cowplot::plot_grid( stab, resil, resist, synchplot, ncol=1 )
cowplot::plot_grid( stabresil, resistsynch, ncol=1 )
ggsave( "R Code and Analysis/Figs/stability+resist_resil_synch.svg", width=6, height=6 )





library(purrr)
library(broom)
ress_long %>%
  nest(-measure) %>% 
  mutate(
    fit = map(data, ~ lm( log(value) ~ gmeanr, data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)

ress_long2 <- ress %>%
  gather( "resistance","value", -Site, -Zone, -source, -D1, -D2 )

# show correlations between resistance and resilience

ress_long2 %>%
  nest(-resistance,-source) %>% 
  mutate(
    fit = map(data, ~ cor.test( (.x$value), (.x$D1) )),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied, .drop=TRUE )

ggplot( ress_long2, aes(x=(value),y=(D1),fill=Site)) +
  facet_grid(source~resistance, scales="free") + geom_point(size=3, alpha=0.5, pch=21) +
  scale_fill_manual(values=c("black","gray50","whitesmoke") ) +
  geom_smooth(aes(group=1),method='lm') +
  theme_bw()
