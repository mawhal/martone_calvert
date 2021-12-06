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

# read in data from Community_rda.R
d.simple <- read_csv("R/output/data_select_rda_HMSC.csv")

# define quadrats
d.simple <- d.simple %>% 
  mutate( quadrat = paste( Site, Zone, Meter.point) )

# read in quadrat overlap data pre and post heatwave
quadsprepost <- read_csv("R/output/quadrats_pre_post_heatwave.csv")
d.select <- d.simple %>% 
  filter( quadrat %in% quadsprepost$quadrat )
unique(quadsprepost$quadrat)
unique(d.simple$quadrat)

# only take data from before and after the heatwaves started
d.select <- d.select %>% filter( Year %in% c(2012,2013,2018,2019))

# # calculate average abundance by transect
# d.trans <- d.simple %>%
#   group_by( Year, Site, Zone,taxon  ) %>%
#   summarize( Abundance=mean(Abundance, na.rm=T))


# spread Taxon column out into many columns filled with abundance/cover data
# d.comm.all  <- d.simple %>%
  # select(-funct_2021) %>% 
  # spread( taxon, Abundance, fill=0 )
d.comm.algae <- d.select %>%
  filter( funct_2021 != "animal" ) %>% 
  select(-funct_2021) %>% 
  spread( taxon, Abundance, fill=0 )
# d.comm.fun.algae <- d.simple %>%
#   filter( funct_2021 != "animal" ) %>% 
#   group_by( UID, Year, Site, Zone, Quadrat, funct_2021 ) %>% 
#   summarize( Abundance = sum(Abundance)) %>% 
#   ungroup() %>% 
#   spread( funct_2021, Abundance, fill=0 )






# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
am.select  <- am %>% 
  mutate( quadrat = paste( Site, Zone, Meter.point) ) %>% 
  select( quadrat, Site, Zone, Meter.point, Shore_height_cm )
muse <- left_join( quadsprepost, distinct(am.select) )
splits <- strsplit( as.character(muse$UID), " " )
muse$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))
muse <- arrange(muse,UID)



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

##Sort metadata and community matrix to be the same order
# d.comm.order <- d.comm[ order(match(d.comm$transect, mtrans$transect)),]
# cbind( d.comm.order$transect, mtrans$transect )

# remove UID column from community data
# comm.all <- as.matrix(d.comm.all[,-c(1:6)])
comm.algae <- as.matrix(d.comm.algae[,-c(1:7)])
# fun.algae <- as.matrix(d.comm.fun.algae[,-c(1:5)])

# add a line to d.comm to get missing samples

# combine in a list
# comm <- list(comm.all,comm.algae,fun.algae)




## Steps
## ENSPIE
# the function
ENSPIE <- function(prop){
  ifelse( sum(prop,na.rm=T)>0, 1 / sum(prop^2, na.rm=T), NA ) 
} 
# Hill-Shannon
Hill_Shannon <- function(prop){
  exp( -sum(prop*log(prop),na.rm=T) )
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
  hillshan = apply( prop, 1, Hill_Shannon )
  evar = apply( z, 1, Evar )
  return( data.frame(total.cover,richness,shannon,simpson,enspie,hillshan,evar) )
}
# divs <- lapply( comm, divcalcs )
divs <- divcalcs( comm.algae )
# divsdf <- bind_rows(divs, .id = "source")
# divsdf$source <- factor( divsdf$source, levels=c("1","2","3"), labels=c("all","algae","functional") )
ddf <- d.comm.algae[,1:7]
# ddf$source <- factor( ddf$source, levels=c("1","2"), labels=c("all","algae") )
dd <- bind_cols( ddf, divs )

# add other metadata
dd <- left_join( dd, muse)

# define
dd$prepost <- ifelse( dd$Year %in% c(2012,2013), "pre", "post")

# take average of pre and post conditions
ddmean <- dd %>% 
  select(-UID, -Year, -Site, -Zone, -Quadrat, -Meter.point, -post, -pre, -both, -Shore_height_cm) %>% 
  group_by(quadrat, prepost) %>% 
  summarize_all( .funs = mean)

ddratio <- ddmean %>%
  group_by( quadrat ) %>% 
  summarize( cover.ratio = total.cover[prepost=="post"]/total.cover[prepost=="pre"],
             richness.ratio = richness[prepost=="post"]/richness[prepost=="pre"],
             enspie.ratio = enspie[prepost=="post"]/enspie[prepost=="pre"],
             hillshan.ratio = hillshan[prepost=="post"]/hillshan[prepost=="pre"],
             cover.var = var( total.cover),
             cover.stab = mean(total.cover)/(var(total.cover)+1), # adding one to keep a single quadrat in the analysis that had ZERO change in total seaweed cover, despite variation across taxa
             cover.prop = (total.cover[prepost=="post"]-total.cover[prepost=="pre"])/ ((total.cover[prepost=="post"]+total.cover[prepost=="pre"])/2)  )
pairs( select(ddratio, cover.ratio, cover.var, cover.stab) )


ddratiolog <- ddratio %>% 
  mutate( cover.ratio = (log(cover.ratio)),
          richness.ratio = (log(richness.ratio)),
          enspie.ratio = (log(enspie.ratio)),
          hillshan.ratio = (log(hillshan.ratio)),
          # cover.var = log(cover.var),
          cover.stab = log(cover.stab) )

pairs( select(ddratiolog, cover.ratio, cover.var, cover.stab) )

ddratiolog$posneg <- sign(ddratiolog$cover.ratio)
ddratiologabs <- ddratiolog
ddratiologabs[2:5] <- abs(ddratiolog[2:5])


ddmeanpre <- ddmean %>% 
  filter( prepost == "pre" ) %>% 
  select( quadrat, richnesspre = richness, enspiepre = enspie, 
          hillshanpre = hillshan)

quadstability <- left_join( left_join( ddratio, ddmeanpre), select(muse, quadrat, Site, Zone, Shore_height_cm))
quadstability <- quadstability %>%
  unite( transect, Site, Zone, remove = F ) %>% 
  mutate( log.cover.ratio = log(cover.ratio),
          log.cover.stab = log(cover.stab))




### analogous measure of synchrony from Leps et al. 2018 Ecology
# numerator is variance in total cover
# denominator is the sum variance of each species
# ratio is log transformed, hence 'logV'

# because we are using two time points in this case, we use the log ratio of cover as a metric for change in total cover or the mean/variance method
# summarize over time
d.algae <- d.comm.algae %>% 
  mutate(prepost = ifelse( d.comm.algae$Year %in% c(2012,2013), "pre", "post") ) %>% 
  pivot_longer( Acrosiphonia:Unknown.crust, names_to = "taxon", values_to = "cover" ) #, -Site, -Zone, -Quadrat, -Meter.point, -quadrat )

d.algae.mean <- d.algae %>% 
  select(-UID, -Year, -Site, -Zone, -Quadrat, -Meter.point ) %>% 
  group_by(quadrat, prepost, taxon) %>% 
  summarize( meancover = mean(cover) ) 

# need to remove taxa that did not appear in either year
taxon.years <- d.algae.mean %>% 
  mutate( presence = ifelse( meancover>0, T, F) ) %>% 
  group_by(quadrat, taxon) %>% 
  summarize( years = sum(presence) ) 


d.algae.both <- left_join( d.algae.mean, taxon.years ) %>% 
  filter( years > 0  )

d.algae.both %>% filter( quadrat == "Foggy Cove HIGH 26" ) %>% 
  group_by(quadrat, prepost) %>% 
  summarize( meancover = sum(meancover) )

denom <- d.algae.both  %>% 
  group_by( quadrat, taxon ) %>%
  summarize( eai = var(meancover, na.rm=T), # variance instead of SD
             cover.ratio = meancover[prepost=="post"]/meancover[prepost=="pre"]) %>%  
  mutate( eai = ifelse( eai>0,eai,NA) ) %>%
  mutate( cover.ratio = ifelse( is.infinite(cover.ratio) | cover.ratio == 0, NA, cover.ratio) ) %>%
  group_by( quadrat ) %>% 
  summarize( Evi=sum(eai, na.rm=T),
             cover.ratio.sum = sum( cover.ratio, na.rm=T ) ) 





synch <- left_join( ddratio, denom )
synch <- synch %>% 
  mutate( logV = log(cover.var+1/Evi+1),# adding 1 to keep in one quadrat (see cover.stab calculation above) 
          synch.cover.ratio = log(cover.ratio / cover.ratio.sum) ) %>%  
  select( quadrat, logV, synch.cover.ratio)

quadstability <- left_join( quadstability, synch )



quadstability$abs.log.cover.ratio <- abs(quadstability$log.cover.ratio )
psych::pairs.panels( select(quadstability, richnesspre, Shore_height_cm, log.cover.ratio, log.cover.stab))

windows(5,5)
par( mar = c(5,4,2,2)+0.1 )
psych::pairs.panels( select(quadstability, richnesspre, Shore_height_cm, logV, log.cover.stab ),
                     hist.col = "whitesmoke" )









a <- ggplot( data = quadstability, 
        aes(x = richnesspre, y = cover.ratio, col = Zone, shape = Site)) + 
  geom_hline( yintercept = 1) +
  # geom_smooth(aes(group=Site),method = 'lm', se = F) +
  geom_smooth(aes(group=1),method = 'lm', se = T) +
  geom_point() +
  scale_y_log10() +
  theme_classic()

ggplot( data = quadstability, 
        aes(x = richnesspre, y = cover.ratio, col = Zone, shape = Site)) + 
  geom_hline( yintercept = 1) +
  geom_smooth(aes(group=transect),method = 'lm', se = F) +
  # geom_smooth(aes(group=1),method = 'lm', se = T) +
  geom_point() +
  scale_y_log10() +
  theme_classic()
# ggsave("R/Figs/quadrat_stability_richness_transect.svg", height = 3, width = 4)

b <- ggplot( data = quadstability, 
        aes(x = Shore_height_cm, y = richnesspre, col = Zone, shape = Site)) + 
  geom_smooth(aes(group=1),method = 'lm', se = T) +
  geom_point() +
  theme_classic()

c <- ggplot( data = quadstability, 
             aes(x = Shore_height_cm, y = cover.ratio, col = Zone, shape = Site)) + 
  geom_hline( yintercept = 1) +
  geom_smooth(aes(group=1),method = 'lm', se = T) +  
  geom_point() +
  scale_y_log10() +
  theme_classic()

d <- ggplot( data = quadstability, 
             aes(x = Shore_height_cm, y = cover.stab, col = Zone, shape = Site)) + 
  geom_smooth(aes(group=1),method = 'lm', se = T) +  
  geom_point() +
  scale_y_log10( ) +
  # coord_cartesian( ylim = c(0.0001,100) ) +
  theme_classic()

cowplot::plot_grid( b, c, d, ncol=1)
ggsave("R/Figs/quadrat_stability_richness_elevation.svg", height = 9, width = 4)



### linear (mixed) models
library(partR2)
library(lme4)
library(lmerTest)
library(car)

# scale predictors
quadstability$richscale <- scale( quadstability$richnesspre )
quadstability$shanscale <- scale( quadstability$hillshanpre )
quadstability$divpick   <- quadstability$richscale
quadstability$elevscale <- scale( quadstability$Shore_height_cm )



lme0 <- lmer( log(cover.stab) ~ 1 + (1 | transect), data = quadstability)
summary(lme0)
lme1 <- lmer( log(cover.stab) ~ divpick + (1 | transect), data = quadstability)
summary(lme1)
lme2 <- lmer( log(cover.stab) ~ divpick+elevscale + (1 | transect), data = quadstability)
summary(lme2)
lme3 <- lmer( log(cover.stab) ~ elevscale + (1 | transect), data = quadstability)
summary(lme3)
lme4 <- lmer( log(cover.stab) ~ divpick+elevscale+logV + (1 | transect), data = quadstability)
summary(lme4)
lme5 <- lmer( log(cover.stab) ~ logV + (1 | transect), data = quadstability)
summary(lme5)
lme6 <- lmer( log(cover.stab) ~ logV+elevscale + (1 | transect), data = quadstability)
summary(lme6)
lme7 <- lmer( log(cover.stab) ~ logV+divpick + (1 | transect), data = quadstability)
summary(lme7)

vif(lme4)
aictable <- bbmle::AICctab( lme0,lme1, lme2, lme3, lme4, lme5, lme6, lme7, nobs = nrow(quadstability), weights = T, delta = T, base = T  )
aictable
write_csv(as.data.frame(aictable), "R/output/stability_diversity_synchrony_aic_quadrat.csv")
anova(lme4)
R2_lme2 <- partR2( lme4, partvars = c("divpick","elevscale","logV"), R2_type = "marginal", nboot = 100)
summary(R2_lme2)
## results (takes a while with so many bootstraps -- 100 used)
# R2 (marginal) and 95% CI for the full model: 
#   R2     CI_lower CI_upper ndf
# 0.9557 0.9204   0.9717   4  
# 
# ----------
#   
#   Part (semi-partial) R2:
#   Predictor(s)           R2     CI_lower CI_upper ndf
# Model                  0.9557 0.9204   0.9717   4  
# divpick                0.0000 0.0000   0.0740   3  
# elevscale              0.1100 0.0333   0.1778   3  
# logV                   0.7320 0.7026   0.7606   3  
# divpick+elevscale      0.1107 0.0341   0.1784   2  
# divpick+logV           0.7327 0.7032   0.7612   2  
# elevscale+logV         0.9455 0.9105   0.9618   2  
# divpick+elevscale+logV 0.9557 0.9204   0.9717   1  
# 
# ----------
#   
#   Inclusive R2 (SC^2 * R2):
#   Predictor IR2    CI_lower CI_upper
# divpick   0.0188 0.0114   0.0357  
# elevscale 0.2238 0.1738   0.2810  
# logV      0.9231 0.8805   0.9540  
# 
# ----------
#   
#   Structure coefficients r(Yhat,x):
#   Predictor SC      CI_lower CI_upper
# divpick    0.1403  0.1091   0.1943 
# elevscale -0.4840 -0.5439  -0.4236 
# logV      -0.9828 -0.9927  -0.9675 
# 
# ----------
#   
#   Beta weights (standardised estimates)
# Predictor BW      CI_lower CI_upper
# divpick    0.0099 -0.0191   0.0514 
# elevscale -0.1898 -0.2618  -0.1112 
# logV      -0.9238 -0.9564  -0.8380 
# 
# ----------
# R2_lme2cond <- partR2( lme2, partvars = c("richscale","elevscale"), R2_type = "conditional", nboot = 10)
# summary(R2_lme2cond)


lme_synch <- lmer( logV ~ divpick+elevscale + (1 | transect), data = quadstability)
summary( lme_synch )
vif( lme_synch )


# write the main data.frame to file for use in another script - transect-level stability
write_csv(quadstability.synch, "R/output/stability_diversity_quadrat.csv")
