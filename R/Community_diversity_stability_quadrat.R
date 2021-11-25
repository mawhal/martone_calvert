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
             hillshan.ratio = hillshan[prepost=="post"]/hillshan[prepost=="pre"] )


ddratiolog <- ddratio %>% 
  mutate( cover.ratio = (log(cover.ratio)),
          richness.ratio = (log(richness.ratio)),
          enspie.ratio = (log(enspie.ratio)),
          hillshan.ratio = (log(hillshan.ratio)))

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
  mutate( log.cover.ratio = log(cover.ratio))


psych::pairs.panels( select(quadstability, richnesspre, Shore_height_cm, cover.ratio))
psych::pairs.panels( select(quadstability, richnesspre, Shore_height_cm, log.cover.ratio))

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
ggsave("R/Figs/quadrat_stability_richness_transect.svg", height = 3, width = 4)

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

cowplot::plot_grid( b, c, a, ncol=1)
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

lm1 <- lm( log10(cover.ratio) ~ divpick, data = quadstability)
summary(lm1)
lm2 <- lm( log10(cover.ratio) ~ elevscale, data = quadstability)
summary(lm2)
lm3 <- lm( log10(cover.ratio) ~ divpick+elevscale, data = quadstability)
summary(lm3)
vif(lm3)
anova(lm3)
modEvA::varPart( A = 0.102, B = 0.2812, AB = 0.2906, A.name = "richness", B.name = "elevation")

lme0 <- lmer( log10(cover.ratio) ~ 1 + (1 | transect), data = quadstability)
summary(lme0)
lme1 <- lmer( log10(cover.ratio) ~ divpick + (1 | transect), data = quadstability)
summary(lme1)
lme2 <- lmer( log10(cover.ratio) ~ divpick+elevscale + (1 | transect), data = quadstability)
summary(lme2)
lme3 <- lmer( log10(cover.ratio) ~ elevscale + (1 | transect), data = quadstability)
summary(lme3)
vif(lme2)
bbmle::AICctab( lme0,lme1, lme2, lme3, nobs = nrow(quadstability) )
anova(lme2)
R2_lme2 <- partR2( lme2, partvars = c("divpick","elevscale"), R2_type = "marginal", nboot = 10)
summary(R2_lme2)
## results (takes a while with so many bootstraps -- 1000 used)
# R2 (marginal) and 95% CI for the full model: 
#   R2     CI_lower CI_upper ndf
# 0.3146 0.1447   0.4987   3  
# 
# ----------
#   
#   Part (semi-partial) R2:
#   Predictor(s)        R2     CI_lower CI_upper ndf
# Model               0.3146 0.1447   0.4987   3  
# richscale           0.0000 0.0000   0.2281   2  
# elevscale           0.2765 0.0977   0.4673   2  
# richscale+elevscale 0.3146 0.1447   0.4987   1  
# 
# ----------
#   
#   Inclusive R2 (SC^2 * R2):
#   Predictor IR2    CI_lower CI_upper
# richscale 0.0896 0.0159   0.2195  
# elevscale 0.3109 0.1288   0.4890  
# 
# ----------
#   
#   Structure coefficients r(Yhat,x):
#   Predictor SC      CI_lower CI_upper
# richscale  0.5336  0.2533   0.7999 
# elevscale -0.9940 -1.0000  -0.8897 
# 
# ----------
#   
#   Beta weights (standardised estimates)
# Predictor BW      CI_lower CI_upper
# richscale  0.0712 -0.1262   0.2908 
# elevscale -0.5504 -0.7107  -0.3159 
# 
# ----------
#   
#   Parametric bootstrapping resulted in warnings or messages:
#   Check r2obj$boot_warnings and r2obj$boot_messages.


# R2_lme2cond <- partR2( lme2, partvars = c("richscale","elevscale"), R2_type = "conditional", nboot = 10)
# summary(R2_lme2cond)


# write the main data.frame to file for use in another script - transect-level stability
write_csv(quadstability, "R/output/stability_diversity_quadrat.csv")
