# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen 
# created 2 Dec 2021

# Using Wang et al. 2019 framework for stability and synchrony
# https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.04290

# 

# 
library(tidyverse)
library(codyn)
library(str2str)
# var.partition function
source( "R/var.partition.R" )

# 

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



# just get algae
d.comm.algae <- d.select %>%
  filter( funct_2021 != "animal" ) %>% 
  select(-funct_2021) %>% 
  spread( taxon, Abundance, fill=0 )





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
# define transects
mclean <- mclean %>% 
  unite( "transect", Site, Zone, remove = F)

##Sort metadata and community matrix to be the same order
# d.comm.order <- d.comm[ order(match(d.comm$transect, mtrans$transect)),]
# cbind( d.comm.order$transect, mtrans$transect )

# remove UID column from community data
# comm.all <- as.matrix(d.comm.all[,-c(1:6)])
comm.algae <- as.matrix(d.comm.algae[,-c(1:7)])



# define
d.comm.algae$prepost <- ifelse( d.comm.algae$Year %in% c(2012,2013), "pre", "post")

dd <- left_join(mclean, d.comm.algae)

# take average of pre and post conditions
ddmean <- dd %>% 
  select(-UID, -Year, -Site, -Zone, -Quadrat, -Meter.point, -post, -pre, -both, -Shore_height_cm, -transect) %>% 
  group_by(quadrat, prepost) %>% 
  summarize_all( .funs = mean)





# remove taxa that do not appear
colSums(ddmean[,-c(1,2)])
ddmean.transect <- ddmean %>% 
  separate( quadrat, c('Site', 'blah', 'Zone', 'meter.point'), sep=" " ) %>% 
  unite( "transect", Site, blah, Zone) %>% 
  mutate( transect = factor(transect) )





### var.partition function requires a species x time x community array
ddcom <- ddmean[,-c(1,2)]



# how to define the metacommunities?
# quadrat level data, consider transects as metacommunities? or sites?

## transect as metacommunity
# independently calculate partitioning 
# split into list of transect data
transect <- ddmean.transect$transect
# quadrat <- with( ddmean.transect, paste( transect, meter.point ))
dmt.comm <- ddmean.transect[,-c(1:3)]
dmt.split <- split( ddmean.transect, f = transect)
# abind::abind( dmt.split, along = 3 )
# 


# 
array.list <- lapply( dmt.split, function(z){
  tmp <- z %>% 
    unite( "quadrat", transect, meter.point)
  test.comm <- tmp[, -c(1,2)]
  test.split <- split( test.comm, f= factor(tmp$quadrat) )
  test.array <- ld2a( test.split )
  return( test.array )
} )

# inspect the object
names( array.list )
str( array.list )  # a list of arrays
# each list element in a transect containing a 1 x 93 x Q arrary, 
#    where Q is the number of quadrats sampled in the end-member years of the survey as we defined (2012/2013 and 2018/2019)

# how many total quads in each transect?
quad.number <- unlist(lapply( array.list, function(z) dim(z)[3] ))
quad.number


# var.partition accepts the array so
# provide arrays to var.partition in another lapply call
var.synch.quads.in.transects <- lapply( array.list, function(z) {
  var.partition(z)
})

var.synch.quads.in.transects <- do.call( rbind, var.synch.quads.in.transects )

psych::pairs.panels( var.synch.quads.in.transects, scale = T ) 
psych::pairs.panels( var.synch.quads.in.transects[, c('CV_S_L', 'phi_S2C_L', 'CV_C_L')], scale = T ) # left-hand side of Fig. 1 in Wang et al. 2019



#### Diversity - consider adding element of this to the Community_diversity_elevtime script to consider mean local richness in each year with regional (totals on transects) richness trends
# combine this with comparisons of diversity
# local richness - mean species richness per quad
# regional richness - total species sampled across quads in each transect
# beta = alpha/gamma
ddmean.transect$localrich <- rowSums( ifelse( ddmean.transect[, -c(1:3)] > 0, 1, 0 ) )
# plot changes at quad level, categorized by site and zone
ddmean.transect <- ddmean.transect %>% 
  separate( transect, c("Site","Blah","Zone"), remove = F ) %>% 
  unite( 'Site', Site, Blah, sep = " ") %>% 
  unite( 'quadrat', transect, meter.point, sep = "_", remove = F)
ddmean.transect$Zone <- factor( ddmean.transect$Zone, levels = c("LOW","MID","HIGH") )
ddmean.transect$prepost <- factor( ddmean.transect$prepost, levels = c("pre","post") )
ggplot( ddmean.transect, aes(x = prepost, y = localrich, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = quadrat)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 )


regional <- ddmean.transect %>% 
  select( -meter.point, -quadrat ) %>% 
  group_by( transect, Site, Zone, prepost ) %>% 
  summarize_all( .funs = sum )


regional$regionalrich <- rowSums( ifelse( regional[, -c(1:4)] > 0, 1, 0 ) )

# combine to get a comparison of mean quad local richness and regional richness before and after the heatwaves

local.regional.prepost <- left_join( ddmean.transect, select( regional, transect, prepost, regionalrich) ) %>% 
  mutate( beta = regionalrich/localrich )
local.regional.mean <- local.regional.prepost %>% 
  select( transect, Site, Zone, localrich, regionalrich, beta, prepost ) %>% 
  group_by( transect, Site, Zone, prepost ) %>% 
  summarize( meanlocalrich = mean(localrich), regionalrich = mean(regionalrich), meanbeta = mean(beta) )

ggplot( local.regional.mean, aes( x = prepost, y = meanlocalrich, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = transect)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 ) 
ggplot( local.regional.mean, aes( x = prepost, y = regionalrich, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = transect)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 ) 
ggplot( local.regional.mean, aes( x = prepost, y = meanbeta, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = transect)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 ) 

# average pre and post
local.regional <- local.regional.mean %>% 
  select( -prepost ) %>% 
  group_by( transect, Site, Zone ) %>% 
  summarize_all( .funs = mean ) 

local.regional.pre <- filter( local.regional.mean, prepost == "pre" )
names(local.regional.pre)[ 5:7 ] <- c( 'meanlocalrich.pre', 'regionalrich.pre', 'meanbeta.pre' )

### combine diversity with variance and asynchrony measures
pre.mean <-  left_join(local.regional, local.regional.pre) 

var.div <- bind_cols( pre.mean, as.data.frame(var.synch.quads.in.transects) )
var.div <- var.div %>% 
  mutate( lambda = phi_C_L2R/phi_S_L2R ) %>% 
  mutate( transect = as.character(transect), Zone = as.character(Zone) ) %>% 
  unite( "transect", Site, Zone, remove = F, sep = "_")

# add elevation
transect.elev <- mclean %>% 
  select( transect, Site, Zone, Shore_height_cm ) %>% 
  group_by( transect, Site, Zone ) %>% 
  summarize( mean_elev = mean(Shore_height_cm)) %>% 
  mutate( Site = as.character(Site), Zone = as.character(Zone))
var.div <- left_join( transect.elev, var.div )
var.div$Zone <- factor( var.div$Zone, levels = c("LOW","MID","HIGH") )


# linear models and plots
summary(lm(CV_C_L ~ meanlocalrich, var.div))
a <- ggplot( var.div, aes( x = meanlocalrich, y = CV_C_L ) ) + 
  geom_smooth( method = 'lm', se = F ) +
  geom_point()
summary(lm(phi_S2C_L ~ meanlocalrich, var.div))
b <- ggplot( var.div, aes( x = meanlocalrich, y = phi_S2C_L ) ) + 
  geom_smooth( method = 'lm', se = F ) +
  geom_point()
summary(lm(CV_S_L ~ meanlocalrich, var.div))
c <- ggplot( var.div, aes( x = meanlocalrich, y = CV_S_L ) ) + 
  geom_smooth( method = 'lm', se = F ) +
  geom_point()
summary(lm(CV_C_R ~ regionalrich, var.div))
d <- ggplot( var.div, aes( x = regionalrich, y = CV_C_R ) ) + 
  geom_smooth( method = 'lm', se = F ) + 
  geom_point()
summary(lm(phi_S2C_R ~ regionalrich, var.div))
e <- ggplot( var.div, aes( x = regionalrich, y = phi_S2C_R ) ) + 
  geom_smooth( method = 'lm', se = F ) + 
  geom_point()
summary(lm(CV_S_R ~ regionalrich, var.div))
f <- ggplot( var.div, aes( x = regionalrich, y = CV_S_R ) ) + 
  geom_smooth( method = 'lm', se = F ) + 
  geom_point()


cowplot::plot_grid( a, d, b, e, c, f,  nrow = 3,
                    labels = "auto")
ggsave( "R/Figs/meta_stability_synchrony_richness.svg", width = 6, height = 9 )


summary(lm(phi_C_L2R ~ meanbeta, var.div))
a <- ggplot( var.div, aes( x = meanbeta, y = phi_C_L2R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) + 
  geom_point()
summary(lm(lambda ~ meanbeta, var.div))
b <- ggplot( var.div, aes( x = meanbeta, y = lambda ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) + 
  geom_point()
cowplot::plot_grid( a, b,  nrow = 1,
                    labels = "auto")
ggsave( "R/Figs/meta_stability_synchrony_beta.svg", width = 6, height = 3 )



# repeat for elevation
summary(lm(CV_C_L ~ mean_elev, var.div))
a <- ggplot( var.div, aes( x = mean_elev, y = CV_C_L ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) +
  geom_point()
summary(lm(phi_S2C_L ~ mean_elev, var.div))
b <- ggplot( var.div, aes( x = mean_elev, y = phi_S2C_L ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) +
  geom_point()
summary(lm(CV_S_L ~ mean_elev, var.div))
c <- ggplot( var.div, aes( x = mean_elev, y = CV_S_L ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 1 ) +
  geom_point()
summary(lm(CV_C_R ~ mean_elev, var.div))
d <- ggplot( var.div, aes( x = mean_elev, y = CV_C_R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) + 
  geom_point()
summary(lm(phi_S2C_R ~ mean_elev, var.div))
e <- ggplot( var.div, aes( x = mean_elev, y = phi_S2C_R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) + 
  geom_point()
summary(lm(CV_S_R ~ mean_elev, var.div))
f <- ggplot( var.div, aes( x = mean_elev, y = CV_S_R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) + 
  geom_point()


cowplot::plot_grid( a, d, b, e, c, f,  nrow = 3,
                    labels = "auto")
ggsave( "R/Figs/meta_stability_synchrony_elevation.svg", width = 6, height = 9 )


# compare synchrony at the different levels


my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs( select(ungroup(var.div), CV_S_L, CV_C_L, CV_S_R, CV_C_R), lower.panel = my_line, upper.panel = panel.cor, xlim = c(2,6), ylim = c(2,6) )
# metapopulation variability explains metacommunity variability (or, population variability explains transect variability), they are the same
# community variability of seaweed cover in transects is less that population variability 
# metacommunity variability is much lower than population variability


pairs( select(ungroup(var.div), phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), lower.panel = my_line, upper.panel = panel.cor, xlim = c(0.5,1), ylim = c(0.5,1) )
# community level spatial synchrony was typically a bit higher than species-level spatial synchrony, 8 of 9 transects
# regional-scale species synchrony was typically higher than local-scale species synchrony; both we pretty high

# species synchrony is generally higher than spatial synchrony, even though these are pretty small areas (transects are small metacommunities)
# suggests that spatial asynchrony reduces community and metacommunity variability despite high species synchrony (synchrony of species in communities at both local to regional levels)
# indeed species synchrony is very high
# suggests a strong role for spatial heterogeneity?

# these are sort of the opposite of what Wnag et al. 2019 found in grassland systems

# the relative strength of local richness and regional richness to species variability, suggest that species-level variability is what drives variabilty in community  cover
# but the high synchrony suggest that heatwave effects were strong?
