# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen 
# created 2 Dec 2021

# Using Wang et al. 2019 framework for stability and synchrony
# https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.04290

# 

# 
library(tidyverse)
library(str2str)
## var.partition function from Wang et al. 2019
source( "R/var.partition.R" )

## Evenness as defined as Evar in Smith & Wilson 1996 Oikos
Evar <- function( x ){
  S = length( x[x>0] )
  1 - 2/pi*atan( sum((log(x[ x>0 ]) - sum(log(x[ x>0 ]))/S)^2)/S ) 
}


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
# unique(quadsprepost$quadrat)
# unique(d.simple$quadrat)

# only take data from before and after the heatwaves started
d.select <- d.select %>% filter( Year %in% c(2012,2013,2018,2019))






# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
am.select  <- am %>% 
  mutate( quadrat = paste( Site, Zone, Meter.point, sep = " ") ) %>% 
  select( quadrat, Site, Zone, Meter.point, Shore_height_cm )
muse <- left_join( quadsprepost, distinct(am.select) )



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
  unite( "transect", Site, Zone, remove = F, sep = " ")

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
  unite( "transect", Site, blah, Zone, sep = " ") %>% 
  mutate( transect = factor(transect) )






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


ddmean.transect$total.cover <- rowSums(dmt.comm)

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
ddmean.transect$evenness <- apply( ddmean.transect[, -c(1:3)], 1, Evar )
boxplot( ddmean.transect$evenness) # pretty low evenness on average, but every possibility
# plot changes at quad level, categorized by site and zone
ddmean.transect <- ddmean.transect %>% 
  separate( transect, c("Site","Blah","Zone"), remove = F, sep = " " ) %>% 
  unite( 'Site', Site, Blah, sep = " ") %>% 
  unite( 'quadrat', transect, meter.point, sep = " ", remove = F) %>% 
  mutate( transect = as.character(transect) )
ddmean.transect$Zone <- factor( ddmean.transect$Zone, levels = c("LOW","MID","HIGH") )
ddmean.transect$prepost <- factor( ddmean.transect$prepost, levels = c("pre","post") )
ggplot( ddmean.transect, aes(x = prepost, y = localrich, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = quadrat)) +
  geom_smooth( aes( group = Zone ), method = 'lm', lwd = 3 )

# add elevation
quad.elev <- mclean %>% 
  select( quadrat, transect, Site, Zone, Shore_height_cm ) %>% 
  # group_by( transect, Site, Zone ) %>% 
  # summarize( meanelev = mean(Shore_height_cm)) %>% 
  mutate( Site = as.character(Site), Zone = as.character(Zone))
quad.elev$Zone <- factor( quad.elev$Zone, levels = c("LOW","MID","HIGH") )

ddmean.transect <- left_join( ddmean.transect, quad.elev )
ggplot( ddmean.transect, aes(x = localrich, y = evenness, col = Zone, shape = Site)) + 
  geom_point() + geom_path( aes(group = quadrat)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 )
ggplot( ddmean.transect, aes(x = Shore_height_cm, y = evenness, col = Zone, shape = Site)) + 
  facet_wrap(~prepost) +
  geom_point() + geom_path( aes(group = quadrat)) +
  geom_smooth( aes( group = prepost ), method = 'lm', lwd = 3 )
ggplot( ddmean.transect, aes(x = Shore_height_cm, y = localrich, col = Zone, shape = Site)) + 
  facet_wrap(~prepost) +
  geom_point() + geom_path( aes(group = quadrat)) +
  geom_smooth( aes( group = prepost ), method = 'lm', lwd = 3 )
lme1 <- lme4::lmer( localrich ~ Shore_height_cm*prepost + (1|transect), data = ddmean.transect )
summary(lme1)

# temporal beta diversity using the same array format
# bray-curtis
braydiss <- lapply( array.list, function(z) {
  apply( z, 3, vegan::vegdist )
})
names(braydiss) <- NULL
braydiss <- unlist( braydiss )
braydiss <- data.frame( quadrat = names(braydiss), braydiss )
braydiss$quadrat <- gsub( "_", " ", braydiss$quadrat)
ddmean.transect <- right_join( ddmean.transect, braydiss )

summary( lm( braydiss ~ Shore_height_cm, data = ddmean.transect) )
a <- ggplot( ddmean.transect, aes( x = Shore_height_cm, y = braydiss, col = Zone )) + 
  geom_point() + 
  geom_smooth( aes(group = Zone), method = 'lm', se = F ) +
  geom_smooth( aes(group = 1), method = 'lm') 
summary( lm( braydiss ~ localrich, data = ddmean.transect) )
b <- ggplot( ddmean.transect, aes( x = localrich, y = braydiss, col = Zone )) + geom_point()
ggplot( ddmean.transect, aes( x = Shore_height_cm, y = localrich, col = Zone )) + geom_point()

cowplot::plot_grid( a, b )

#



regional <- ddmean.transect %>% 
  select( -meter.point, -quadrat ) %>% 
  group_by( transect, Site, Zone, prepost ) %>% 
  summarize_all( .funs = sum ) %>% 
  ungroup()


regional$regionalrich <- rowSums( ifelse( select(regional, Acrosiphonia:Unknown.crust) > 0, 1, 0 ) )

# combine to get a comparison of mean quad local richness and regional richness before and after the heatwaves

local.regional.prepost <- left_join( ddmean.transect, select( regional, transect, prepost, regionalrich) ) %>% 
  mutate( beta = regionalrich/localrich )
local.regional.mean <- local.regional.prepost %>% 
  select( transect, Site, Zone, localrich, regionalrich, beta, evenness, prepost, Shore_height_cm ) %>% 
  group_by( transect, Site, Zone, prepost ) %>% 
  summarize( meanlocalrich = mean(localrich), regionalrich = mean(regionalrich), meanbeta = mean(beta), meanelev = mean(Shore_height_cm), meaneven = mean(evenness) )

ggplot( local.regional.mean, aes( x = prepost, y = meanlocalrich, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = transect)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 ) 
ggplot( local.regional.mean, aes( x = prepost, y = regionalrich, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = transect)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 ) 
ggplot( local.regional.mean, aes( x = prepost, y = meanbeta, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = transect)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 ) 
ggplot( local.regional.mean, aes( x = prepost, y = meaneven, col = Zone, shape = Site)) +  geom_point() + geom_path( aes(group = transect)) +
  geom_smooth( aes( group = 1 ), method = 'lm', lwd = 3 ) 


## look at change in variability in elevation optima before and after
cover_quads_all <-  read_csv( "R/output/cover_quads_all_elevvar.csv" ) %>% 
  select( transect, quadrat, Year, elev.var) %>% 
  filter( Year %in% c(2012,2013,2018,2019))
cover_quads_all$prepost <- ifelse( cover_quads_all$Year %in% c(2012,2013), "pre", "post" )
cover_quads_all$prepost <- factor(cover_quads_all$prepost, levels = c('pre','post'))
cover_quads_prepost <- cover_quads_all %>% 
  group_by( transect, quadrat, prepost) %>% 
  summarize( elev.var = mean(elev.var)) %>% 
  separate( quadrat, c("Site","Blah","Zone","meter.point"), remove = F) %>% 
  unite( "Site", Site:Blah)
ggplot( cover_quads_prepost, aes( x = prepost, y = elev.var, group = quadrat, col = Zone, shape = Site)) +
  geom_point() + 
  # geom_line() + 
  geom_smooth( aes( group = Zone), method = 'lm' )
summary(lme4::lmer( elev.var ~ prepost*Zone + (1|transect), data = cover_quads_prepost ))
# 

# average pre and post
local.regional <- local.regional.mean %>% 
  select( -prepost ) %>% 
  group_by( transect, Site, Zone ) %>% 
  summarize_all( .funs = mean ) 

local.regional.pre <- filter( local.regional.mean, prepost == "pre" )
names(local.regional.pre)[ 5:9 ] <- c( 'meanlocalrich.pre', 'regionalrich.pre', 'meanbeta.pre', 'meanelev.pre', 'meaneven.pre' )

### combine diversity with variance and asynchrony measures
pre.mean <-  left_join(local.regional, local.regional.pre) 

var.div <- bind_cols( pre.mean, as.data.frame(var.synch.quads.in.transects) )
var.div <- var.div %>% 
  mutate( lambda = phi_C_L2R/phi_S_L2R ) %>% 
  mutate( transect = as.character(transect), Zone = as.character(Zone) ) %>% 
  unite( "transect", Site, Zone, remove = F, sep = " ")




# look at variability as a function of mean total cover along transects before the heatwave
transect.cover <- ddmean.transect %>% 
  filter(prepost == "pre") %>% 
  group_by(transect) %>% 
  summarize( mean.total.cover = mean(total.cover))
var.div <- left_join( var.div, transect.cover )

ggplot( var.div, aes( x = mean.total.cover, y = CV_S_L) ) + geom_point()
ggplot( var.div, aes( x = mean.total.cover, y = CV_C_L) ) + geom_point()

# linear models and plots
summary(lm(CV_C_L ~ meanlocalrich, var.div))
a <- ggplot( var.div, aes( x = meanlocalrich, y = CV_C_L ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black' ) +
  geom_point() +
  xlab("Mean local species richness") +
  ylab( expression(atop( "Local-scale community",paste("variability (",italic(CV["C,L"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()
summary(lm(phi_S2C_L ~ meanlocalrich, var.div))
b <- ggplot( var.div, aes( x = meanlocalrich, y = phi_S2C_L ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black' ) +
  geom_point() +
  xlab("Mean local species richness") +
  ylab( expression(atop( "Local-scale species",paste("synchrony (",italic(phi[S %->% "C,L"]),")" ))) ) +
  ylim( c(0.79,1) ) +
  theme_classic()
summary(lm(CV_S_L ~ meanlocalrich, var.div))
c <- ggplot( var.div, aes( x = meanlocalrich, y = CV_S_L ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black' ) +
  geom_point()  +
  xlab("Mean local species richness") +
  ylab( expression(atop( "Local-scale species",paste("variability (",italic(CV["S,L"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()
summary(lm(CV_C_R ~ regionalrich, var.div))
d <- ggplot( var.div, aes( x = regionalrich, y = CV_C_R ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black' ) + 
  geom_point()  +
  xlab("Regional species richness") +
  ylab( expression(atop( "Metacommunity variability",paste("(",italic(CV["C,R"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()
summary(lm(phi_S2C_R ~ regionalrich, var.div))
e <- ggplot( var.div, aes( x = regionalrich, y = phi_S2C_R ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black' ) + 
  geom_point()  +
  xlab("Regional species richness") +
  ylab( expression(atop( "Regional-scale species",paste("synchrony (",italic(phi[S %->% "C,R"]),")" ))) ) +
  ylim( c(0.79,1) ) +
  theme_classic()
summary(lm(CV_S_R ~ regionalrich, var.div))
f <- ggplot( var.div, aes( x = regionalrich, y = CV_S_R ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black' ) + 
  geom_point()  +
  xlab("Regional species richness") +
  ylab( expression(atop( "Metapopulation variability",paste("(",italic(CV["S,R"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()


cowplot::plot_grid( a, d, b, e, c, f,  nrow = 3,
                    labels = "auto", align = 'hv' )
ggsave( "R/Figs/meta_stability_synchrony_richness.svg", width = 6.5, height = 9 )


summary(lm(phi_C_L2R ~ meanbeta, var.div))
a <- ggplot( var.div, aes( x = meanbeta, y = phi_C_L2R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2, col = 'black' ) + 
  geom_point() +
  xlab( expression( atop("Spatial community dissimilarity", "(beta diversity)") )) +
  ylab( expression(atop( "Community-level spatial",paste("synchrony (",italic(phi["C,L" %->% R]),")" ))) ) +
  theme_classic()
summary(lm(lambda ~ meanbeta, var.div))
b <- ggplot( var.div, aes( x = meanbeta, y = lambda ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2, col = 'black' ) + 
  geom_point() +
  xlab( expression( atop("Spatial community dissimilarity", "(beta diversity)") )) +
  ylab( expression( atop("Synchrony ratio", paste( "across levels (",italic(lambda),")" ))) ) +
  theme_classic()
cowplot::plot_grid( a, b,  nrow = 1,
                    labels = "auto")
ggsave( "R/Figs/meta_stability_synchrony_beta.svg", width = 6, height = 3 )



# repeat for elevation
summary(lm(CV_C_L ~ meanelev, var.div))
a <- ggplot( var.div, aes( x = meanelev, y = CV_C_L ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black', lty = 2 ) +
  geom_point() +
  xlab("Mean transect elevation (cm)") +
  ylab( expression(atop( "Local-scale community",paste("variability (",italic(CV["C,L"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()
summary(lm(phi_S2C_L ~ meanelev, var.div))
b <- ggplot( var.div, aes( x = meanelev, y = phi_S2C_L ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black', lty = 2 ) +
  geom_point() +
  xlab("Mean transect elevation (cm)") +
  ylab( expression(atop( "Local-scale species",paste("synchrony (",italic(phi[S %->% "C,L"]),")" ))) ) +
  ylim( c(0.79,1) ) +
  theme_classic()
summary(lm(CV_S_L ~ meanelev, var.div))
c <- ggplot( var.div, aes( x = meanelev, y = CV_S_L ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black', lty = 2 ) +
  geom_point()  +
  xlab("Mean transect elevation (cm)") +
  ylab( expression(atop( "Local-scale species",paste("variability (",italic(CV["S,L"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()
summary(lm(CV_C_R ~ meanelev, var.div))
d <- ggplot( var.div, aes( x = meanelev, y = CV_C_R ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black', lty = 2 ) + 
  geom_point()  +
  xlab("Mean transect elevation (cm)") +
  ylab( expression(atop( "Metacommunity variability",paste("(",italic(CV["C,R"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()
summary(lm(phi_S2C_R ~ meanelev, var.div))
e <- ggplot( var.div, aes( x = meanelev, y = phi_S2C_R ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black', lty = 2 ) + 
  geom_point()  +
  xlab("Mean transect elevation (cm)") +
  ylab( expression(atop( "Regional-scale species",paste("synchrony (",italic(phi[S %->% "C,R"]),")" ))) ) +
  ylim( c(0.79,1) ) +
  theme_classic()
summary(lm(CV_S_R ~ meanelev, var.div))
f <- ggplot( var.div, aes( x = meanelev, y = CV_S_R ) ) + 
  geom_smooth( method = 'lm', se = F, col = 'black', lty = 2 ) + 
  geom_point()  +
  xlab("Mean transect elevation (cm)") +
  ylab( expression(atop( "Metapopulation variability",paste("(",italic(CV["S,R"]),")" ))) ) +
  ylim( c(2.3,5.81) ) +
  theme_classic()
cowplot::plot_grid( a, d, b, e, c, f,  nrow = 3,
                    labels = "auto", align = 'hv' )
ggsave( "R/Figs/meta_stability_synchrony_elevation.svg", width = 6, height = 9 )

# richness and elevation
summary(lm(meanlocalrich ~ meanelev, var.div))
ggplot( var.div, aes( x = meanelev, y = meanlocalrich ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 1 ) + 
  geom_point()
summary(lm(regionalrich ~ meanelev, var.div))
ggplot( var.div, aes( x = meanelev, y = regionalrich ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) + 
  geom_point()


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

vars.pairs <- select(ungroup(var.div), CV_S_L, CV_C_L, CV_S_R, CV_C_R)
names( vars.pairs ) <- c( "species\nvariability", "community\nvariability", "metapopulation\nvariability", "metacommunity\nvariability" )
svg( "R/Figs/variability_pairs.svg", width = 4.2, height = 4.2 )
pairs( vars.pairs, lower.panel = my_line, upper.panel = panel.cor, xlim = c(2,6), ylim = c(2,6) )
dev.off()
# metapopulation variability explains metacommunity variability (or, population variability explains transect variability), they are the same
# community variability of seaweed cover in transects is less that population variability 
# metacommunity variability is much lower than population variability

phi.pairs <- select(ungroup(var.div), phi_S2C_L, phi_S2C_R, phi_S_L2R, phi_C_L2R)
names( phi.pairs ) <- c( "local\nspecies\nsynchrony", "regional\nspecies\nsynchrony", "species\nspatial\nsynchrony", "community\nspatial\nsynchrony" )
svg( "R/Figs/synchony_pairs.svg", width = 4.2, height = 4.2 )
pairs( phi.pairs, lower.panel = my_line, upper.panel = panel.cor, xlim = c(0.5,1), ylim = c(0.5,1) )
dev.off()
# community level spatial synchrony was typically a bit higher than species-level spatial synchrony, 8 of 9 transects
# regional-scale species synchrony was typically higher than local-scale species synchrony; both we pretty high

# species synchrony is generally higher than spatial synchrony, even though these are pretty small areas (transects are small metacommunities)
# suggests that spatial asynchrony reduces community and metacommunity variability despite high species synchrony (synchrony of species in communities at both local to regional levels)
# indeed species synchrony is very high
# suggests a strong role for spatial heterogeneity?

# these are sort of the opposite of what Wnag et al. 2019 found in grassland systems

# the relative strength of local richness and regional richness to species variability, suggest that species-level variability is what drives variabilty in community  cover
# but the high synchrony suggest that heatwave effects were strong?




# Now come back to HMSC analysis. Let's consider the peak elevation of a species at the start of the survey (first two years?)
elev_abun_shifts <- read_csv( "R/output/shifts_predicted.csv" )

# compare variability of elevation optima within transects to stability and synchrony
# replace occurrences with initial elevation peaks
# variability as a mean-variance relationship?
initial_peaks_taxon <- elev_abun_shifts %>% 
  mutate( elev.final.med = elev.init.med + elev.shifts.med ) %>% 
  select( taxon, elev.init.med, elev.final.med ) 



# compare variability using composition data from first two years
cover_quads_pre <- ddmean.transect %>% filter( prepost == 'pre')

cover_quads_pre_comm <- select( cover_quads_pre, Alaria.marginata:Ulva )
cover_quads_pre_pa <- ifelse( cover_quads_pre_comm > 0, 1, 0)

# only select columns for which we have results from HMSC
pa_hmsc <- cover_quads_pre_pa[, colnames(cover_quads_pre_pa) %in% initial_peaks_taxon$taxon ]
# filter out peak elevations missing from surveyed quadrats
# interesting that we are missing 3 taxa from 85 quadrats sampled over four years...but, goes to show ya, eh?
initial_peaks_taxon <- filter( initial_peaks_taxon, taxon %in% colnames(pa_hmsc) )

# get things in same order
initial_peaks_taxon <- initial_peaks_taxon %>% 
  arrange(taxon)


# for loop
# add elevation optima for before and after - noting that post heatwave rows are odd and pre heatwave is even
all( initial_peaks_taxon$taxon == colnames(pa_hmsc) )
for(i in 1:ncol(pa_hmsc) ){
    pa_hmsc[,i] <-  ifelse( pa_hmsc[,i] == 1, initial_peaks_taxon$elev.init.med[i], NA )
  }


# calculate mean and sd of peak elevation for each quadrat
peaks.mean <- apply( pa_hmsc, 1, mean, na.rm = T)
peaks.sd <- apply( pa_hmsc, 1, sd, na.rm = T)
boxplot(peaks.sd/peaks.mean)


var.div


cover_quads_pre$elev.var <- peaks.sd/peaks.mean
psych::pairs.panels( select(cover_quads_pre, Shore_height_cm, localrich, evenness, elev.var) )
elev.var.mean <- cover_quads_pre %>% 
  group_by( transect ) %>% 
  summarize( mean.elev.var = mean(elev.var, na.rm = T) )

var.div$transect <- gsub( "_", " ", x = var.div$transect  )
elev.var.mean$transect

# could apply this to every quadrat, but could not relate to stability/synchrony

var.div.elev <- left_join( elev.var.mean, var.div )


var.div.elev$Zone <- factor( var.div.elev$Zone, levels = c("LOW","MID","HIGH"))
var.div.elev$Site <- factor( var.div.elev$Site, levels = c("Foggy Cove", "Fifth Beach", "North Beach" ))

summary( lm( CV_S_L ~ mean.elev.var, data = var.div.elev))
ggplot( var.div.elev, aes( x = mean.elev.var, CV_S_L, col = Zone, shape = Site ) ) +
  geom_point( size = 3 ) +
  geom_smooth( aes(group = 1), method = "lm", se = F )

# higher variability in peak elevation across species (mean aross quads in a transect) is associated with less variability in local conver 

summary(lm(CV_C_L ~ mean.elev.var, var.div.elev))
a <- ggplot( var.div.elev, aes( x = mean.elev.var, y = CV_C_L ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 1 ) +
  geom_point()
summary(lm(phi_S2C_L ~ mean.elev.var, var.div.elev))
b <- ggplot( var.div.elev, aes( x = mean.elev.var, y = phi_S2C_L ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) +
  geom_point()
summary(lm(CV_S_L ~ mean.elev.var, var.div.elev))
c <- ggplot( var.div.elev, aes( x = mean.elev.var, y = CV_S_L ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 1 ) +
  geom_point()
summary(lm(CV_C_R ~ mean.elev.var, var.div.elev))
d <- ggplot( var.div.elev, aes( x = mean.elev.var, y = CV_C_R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 1 ) + 
  geom_point()
summary(lm(phi_S2C_R ~ mean.elev.var, var.div.elev))
e <- ggplot( var.div.elev, aes( x = mean.elev.var, y = phi_S2C_R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 2 ) + 
  geom_point()
summary(lm(CV_S_R ~ mean.elev.var, var.div.elev))
f <- ggplot( var.div.elev, aes( x = mean.elev.var, y = CV_S_R ) ) + 
  geom_smooth( method = 'lm', se = F, lty = 1 ) + 
  geom_point()


cowplot::plot_grid( a, d, b, e, c, f,  nrow = 3,
                    labels = "auto")
ggsave( "R/Figs/meta_stability_synchrony_peak_elev_var.svg", width = 6, height = 9 )



summary( lm( meanlocalrich ~ mean.elev.var, data = var.div.elev))
ggplot( var.div.elev, aes( x = mean.elev.var, meanlocalrich, col = Zone, shape = Site ) ) +
  geom_point( size = 3 ) +
  geom_smooth( aes(group = 1), method = "lm", se = F )



psych::pairs.panels( select(var.div.elev, meanelev, meanlocalrich, mean.elev.var ) )


# models with both variabilty in optima and richness
summary(lm(CV_S_L ~ meanlocalrich+mean.elev.var, var.div.elev))
summary(lm(CV_C_L ~ meanlocalrich+mean.elev.var, var.div.elev))
summary(lm(CV_S_R ~ regionalrich+mean.elev.var, var.div.elev))
summary(lm(CV_C_R ~ regionalrich+mean.elev.var, var.div.elev))

# models of optima variability per species
var.div.elev <- var.div.elev %>% 
  mutate( mean.elev.var.per.sp = mean.elev.var/meanlocalrich )
psych::pairs.panels( select(var.div.elev, meanelev, meanlocalrich, mean.elev.var.per.sp ) )

summary(lm(CV_S_L ~ mean.elev.var.per.sp, var.div.elev))
summary(lm(CV_C_L ~ mean.elev.var.per.sp, var.div.elev))




### plots for manuscript
var.div.elev$Zone <- factor( var.div.elev$Zone, levels = c("LOW","MID","HIGH") )
var.div.elev$Site <- factor( var.div.elev$Site, levels = c("Foggy Cove", "Fifth Beach", "North Beach" ))


localrich_speciesstab <- ggplot( var.div.elev, aes(x = meanlocalrich,y = CV_S_L, shape = Site, fill = Zone)) + 
  # geom_smooth( data = filter( qt, level == "quadrat"), aes( group = level ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  # geom_line( data = filter( qt, level == "transect"), aes(group = Site), size = 0.5, lty = 2) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  ylab( "Local-scale species\nvariability (CV_S,L)" ) +
  xlab( "Local species richness" ) + 
  # coord_cartesian( ylim = c(0.005,50), xlim = c(-2,10)) +
  # scale_y_log10() +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
localrich_speciesstab

localrich_commstab <- ggplot( var.div.elev, aes(x = meanlocalrich,y = CV_C_L, shape = Site, fill = Zone)) + 
  # geom_smooth( data = filter( qt, level == "quadrat"), aes( group = level ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  # geom_line( data = filter( qt, level == "transect"), aes(group = Site), size = 0.5, lty = 2) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  ylab( "Local-scale community\nvariability (CV_C,L)" ) +
  xlab( "Local species richness" ) + 
  # coord_cartesian( ylim = c(0.005,50), xlim = c(-2,10)) +
  # scale_y_log10() +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
localrich_commstab

elevvar_speciesstab <- ggplot( var.div.elev, aes(x = mean.elev.var,y = CV_S_L, shape = Site, fill = Zone)) + 
  # geom_smooth( data = filter( qt, level == "quadrat"), aes( group = level ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  # geom_line( data = filter( qt, level == "transect"), aes(group = Site), size = 0.5, lty = 2) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  ylab( "Local-scale species\nvariability (CV_C,L)" ) +
  xlab( "Variability in species\nelevation optima (CV_elev)" ) + 
  # coord_cartesian( ylim = c(0.005,50), xlim = c(-2,10)) +
  # scale_y_log10() +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elevvar_speciesstab

elevvar_localrich <- ggplot( var.div.elev, aes(x = meanlocalrich,y = mean.elev.var, shape = Site, fill = Zone)) + 
  # geom_smooth( data = filter( qt, level == "quadrat"), aes( group = level ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  # geom_line( data = filter( qt, level == "transect"), aes(group = Site), size = 0.5, lty = 2) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab( "Local species richness" ) + 
  ylab( "Variability in elevation optima" ) + 
  # coord_cartesian( ylim = c(0.005,50), xlim = c(-2,10)) +
  # scale_y_log10() +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elevvar_localrich

cowplot::plot_grid( localrich_speciesstab, elevvar_localrich, elevvar_speciesstab, 
                    nrow = 1, labels = "auto", align = 'hv' )
ggsave( "R/Figs/meta_stability_richness_elevvar.svg", width = 8.5, height = 3)



localrich_localstab <- ggplot( var.div.elev, aes(x = meanlocalrich,y = CV_S_L, shape = Site, fill = Zone)) + 
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( y = CV_C_L, group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE, lty = 2) +
  geom_point( show.legend = TRUE, size = 4 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke"), guide = "none" ) +
  # guides(fill = guide_legend(override.aes = list(shape = 21))) +
  guides(shape = guide_legend(override.aes = list(size=3))) +
  ylab( "Quadrat cover variability (CV)" ) +
  xlab( "Mean quadrat species richness" ) + 
  coord_cartesian( ylim = c(2,6)) +
  annotate( "text", x = 11.5, y = 4.9, label = "total cover\nvariability", 
            hjust = 1, vjust = 1, size = 3 ) +
  annotate( "text", x = 13, y = 5.1, label = "species\nvariability", 
            hjust = 0, vjust = -0.2, size = 3 ) +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank(),
                           legend.text = element_text(size=8),
                           legend.title = element_text(size=9)) 
localrich_localstab

localrich_synchrony <- ggplot( var.div.elev, aes(x = meanlocalrich, y = phi_S2C_L, shape = Site, fill = Zone)) + 
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( y = phi_S2C_R, group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE, lty = 2) +
  geom_point( show.legend = TRUE, size = 4 ) +
  scale_shape_manual( values=c(21,22,24), guide = "none") +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21,size=3)) ) +
  ylab( "Species synchrony" ) +
  xlab( "Mean quadrat species richness" ) + 
  coord_cartesian( ylim = c(0.75,1)) +
  annotate( "text", x = 13, y = 0.9, label = "quadrat\nscale", 
            hjust = 1, vjust = 1, size = 3 ) +
  annotate( "text", x = 15, y = 0.95, label = "transect\nscale", 
            hjust = 0, vjust = -0.2, size = 3 ) +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank(),
                           legend.text = element_text(size=8),
                           legend.title = element_text(size=9)) 
localrich_synchrony

regionalrich_regionalstab <- ggplot( var.div.elev, aes(x = regionalrich,y = CV_S_R, shape = Site, fill = Zone)) + 
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( y = CV_C_R, group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE, lty = 2) +
  geom_point( show.legend = FALSE, size = 4 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  ylab( "Regional-scale variability (CV)" ) +
  xlab( "Regional species richness" ) + 
  coord_cartesian( ylim = c(2,6)) +
  annotate( "text", x = 30, y = 3.9, label = "metacommunity\nvariability", 
            hjust = 1, vjust = 1, size = 3 ) +
  annotate( "text", x = 32, y = 4.1, label = "metapopulation\nvariability", 
            hjust = 0, vjust = -0.2, size = 3 ) +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
regionalrich_regionalstab
cowplot::plot_grid( localrich_localstab, regionalrich_regionalstab)


elev_rich <- ggplot( var.div.elev, aes(x = meanelev,y = meanlocalrich, shape = Site, fill = Zone)) + 
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab( "Elevation (cm)" ) + 
  ylab( "Species richness" ) + 
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elev_rich

elev_stab <- ggplot( var.div.elev, aes(x = meanelev,y = CV_S_L, shape = Site, fill = Zone)) + 
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab( "Elevation (cm)" ) + 
  ylab( "Local-scale variability (CV)" ) + 
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elev_stab

elev_var <- ggplot( var.div.elev, aes(x = meanelev,y = mean.elev.var, shape = Site, fill = Zone)) + 
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab( "Elevation (cm)" ) + 
  ylab( "Variability in elevation optima" ) + 
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elev_var


cowplot::plot_grid(elev_rich, elev_var, elev_stab, ncol = 1)
ggsave( "R/Figs/elevation_richness_stability_var.svg", width = 3, height = 9)

cowplot::plot_grid(elev_rich, elev_var, elev_stab, 
                   NULL,elevvar_localrich,localrich_localstab,
                   ncol = 2, byrow = F )
ggsave( "R/Figs/elevation_richness_stability_var.svg", width = 6, height = 9)


#


elevvar_localrich <- ggplot( var.div.elev, aes(y = meanlocalrich,x = mean.elev.var, shape = Site, fill = Zone)) + 
  # geom_smooth( data = filter( qt, level == "quadrat"), aes( group = level ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  # geom_line( data = filter( qt, level == "transect"), aes(group = Site), size = 0.5, lty = 2) +
  geom_point( show.legend = FALSE, size = 3 ) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  ylab( "Local species richness" ) + 
  xlab( "Variability in elevation optima" ) + 
  # coord_cartesian( ylim = c(0.005,50), xlim = c(-2,10)) +
  # scale_y_log10() +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elevvar_localrich

summary( lm( meanlocalrich ~ mean.elev.var + meanelev, data = var.div.elev ) )







summary( lm( CV_S_L ~ meanlocalrich + mean.elev.var, data = var.div.elev ) )













summary( lm( phi_S2C_L ~ meanlocalrich, data = var.div.elev ) )




summary( lm( CV_C_L ~ CV_S_L + meanlocalrich + phi_S2C_L, data = var.div.elev ) )
car::vif( lm( CV_C_L ~ meanlocalrich + phi_S2C_L, data = var.div.elev ) )
lm1 <-  lm( CV_C_L ~ CV_S_L + meanlocalrich + phi_S2C_L, data = var.div.elev )
lm1 <-  lm( CV_C_L ~ meanlocalrich, data = var.div.elev )
summary( lm( CV_S_R ~ CV_S_L + phi_S_L2R, data = var.div.elev ) )
summary( lm( CV_S_R ~ phi_S_L2R, data = var.div.elev ) )
         
ggplot( data = var.div.elev, aes( x = meanelev, y = CV_S_L, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanelev, y = CV_C_L, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanelev, y = CV_S_R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanelev, y = CV_C_R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanelev, y = phi_S_L2R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanelev, y = phi_C_L2R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanelev, y = phi_S2C_L, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanelev, y = phi_S2C_R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanlocalrich, y = phi_S2C_L, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanlocalrich, y = phi_S2C_R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanlocalrich, y = phi_S_L2R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = meanlocalrich, y = phi_C_L2R, col=Zone)) + geom_point()
ggplot( data = var.div.elev, aes( x = phi_S_L2R, y = CV_S_R)) + geom_point()
ggplot( data = var.div.elev, aes( x = CV_S_L, y = CV_S_R)) + geom_point()


# compare with quadrat-level estimates of elevation, niche (elevation) variability, and richness
# read output data from script "community_diversity_peak_elevation_allquads.R"
cover_quads_all <- read_csv( "R/output/cover_quads_all_elevvar.csv" )
cover_quads_all <- cover_quads_all %>% 
  separate( UID, c("Site", "Blah", "Zone", "Year", "meter.point") ) %>% 
  unite( "Site", Site, Blah, sep = " ") %>% 
  mutate( Site = factor(Site, levels =  c("Foggy Cove", "Fifth Beach", "North Beach" ) )) %>% 
  mutate( Zone = factor(Zone, levels = c("LOW","MID","HIGH")))

cover_quads_all_pre <- cover_quads_all %>% filter( Year == 2012 )
cover_quads_all_post <- cover_quads_all %>% filter( Year == 2019 )

elevvar_localrich <- ggplot( var.div.elev, aes(y = meanlocalrich,x = mean.elev.var, shape = Site, fill = Zone)) + 
  geom_smooth( data = cover_quads_all_pre, aes( y = localrich, x = elev.var, group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_point( data = cover_quads_all_pre, aes( y = localrich, x = elev.var), show.legend = FALSE, alpha = 1, size = 0.75 ) +
  # geom_smooth( data = cover_quads_all_post, aes( y = localrich, x = elev.var, group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  # geom_point( data = cover_quads_all_post, aes( y = localrich, x = elev.var), show.legend = FALSE, alpha = 0.5, size = 0.75 ) +
  # geom_point( show.legend = FALSE, size = 4 ) +
  scale_shape_manual( values=c(21,22,24)) +
  # geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE, lty =2) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  ylab( "Local species richness" ) + 
  xlab( "Variability in elevation optima" ) + 
  # coord_cartesian( ylim = c(0.005,50), xlim = c(-2,10)) +
  # scale_y_log10() +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elevvar_localrich

cor.test( ~ elev.var + localrich, data = cover_quads_all_pre )
localrich_elevvar <- ggplot( cover_quads_all_pre, aes(x = localrich,y = elev.var, shape = Site, fill = Zone)) + 
  geom_point( show.legend = FALSE, alpha = 1, size = 1.5 ) +
  stat_ellipse( aes( group = 1), level = 0.95) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab( "Quadrat species richness" ) + 
  ylab( "Elevation optima variability (CV)" ) + 
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
localrich_elevvar


elev_elevvar <- ggplot( var.div.elev, aes(y = mean.elev.var,x = meanelev, shape = Site, fill = Zone)) + 
  geom_smooth( data = cover_quads_all_pre, aes( y = elev.var, x = Shore_height_cm, group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE) +
  geom_point( data = cover_quads_all_pre, aes( y = elev.var, x = Shore_height_cm), show.legend = FALSE, alpha = 1, size = 1.5 ) +
  # geom_smooth( data = cover_quads_all_post, aes( y = elev.var, x = Shore_height_cm, group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE, lty=2) +
  # geom_point( data = cover_quads_all_post, aes( y = elev.var, x = Shore_height_cm), show.legend = FALSE, alpha = 0.5, size = 0.75 ) +
  scale_shape_manual( values=c(21,22,24)) +
  # geom_smooth( aes( group = 1 ), method="lm", se=F, col='black', lwd=0.5, show.legend = FALSE, lty = 2) +
  # geom_point( show.legend = FALSE, size = 4 ) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab( "Elevation (cm)" ) + 
  ylab( "Elevation optima variability (CV)" ) + 
  # coord_cartesian( ylim = c(0.005,50), xlim = c(-2,10)) +
  # scale_y_log10() +
  theme_classic() + theme( legend.position = c(0.01,0.01), 
                           legend.justification = c(0,0),
                           legend.background = element_blank()) 
elev_elevvar

# slight erosion of elevation variability before and after heatwave, but not significant perhaps due to little change in richness
ggplot( cover_quads_all, aes( x = prepost, y =  elev.var) ) + geom_boxplot(notch = T)


cowplot::plot_grid( elev_elevvar, localrich_elevvar, 
                    localrich_localstab, localrich_synchrony,
                    ncol = 2, labels = "auto", align = 'hv', vjust = 1.05)
ggsave( "R/Figs/elevation_optima_richness_stability.svg", width = 6, height = 6)

