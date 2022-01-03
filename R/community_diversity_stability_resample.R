# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen 
# created 17 Dec 2021

# Using Wang et al. 2019 framework for stability and synchrony
# https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.04290

# adding a resampling piece to try and remove the effect of richness
# grab n species at random, then calculate elevation variability and stability components



# create the array first or sample first?
# 
library(tidyverse)
library(codyn)
library(str2str)
library(quantreg)
## var.partition function
source( "R/var.partition.R" )



## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

# hmsc predicted elevation shifts
elev_abun_shifts <- read_csv( "R/output/shifts_predicted.csv" )


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



# just get algae, and just those species looked at in HMSC to get elevation variability
d.comm.algae <- d.select %>%
  filter( funct_2021 != "animal" ) %>% 
  filter( taxon %in% elev_abun_shifts$taxon ) %>% # just taxa modeled in HMSC
  select(-funct_2021) %>% 
  spread( taxon, Abundance, fill=0 )





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
dmt.split <- split( ddmean.transect, f = transect )
# abind::abind( dmt.split, along = 3 )
# 


# compare variability of elevation optima within transects to stability and synchrony
# replace occurrences with initial elevation peaks
# variability as a mean-variance relationship?
initial_peaks_taxon <- elev_abun_shifts %>% 
  mutate( elev.final.med = elev.init.med + elev.shifts.med ) %>% 
  select( taxon, elev.init.med, elev.final.med ) 



# compare variability using composition data from first two years
cover_quads_all <- ddmean.transect 

cover_quads_all_comm <- select( cover_quads_all, Alaria.marginata:Ulva )
cover_quads_all_pa <- ifelse( cover_quads_all_comm > 0, 1, 0)

# only select columns for which we have results from HMSC
pa_hmsc <- cover_quads_all_pa[, colnames(cover_quads_all_pa) %in% initial_peaks_taxon$taxon ]
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
  for( j in 1:nrow(pa_hmsc) ){
  pa_hmsc[j,i] <-  if( (j %% 2) == 0 ) {
    ifelse( pa_hmsc[j,i] == 1, initial_peaks_taxon$elev.init.med[i], NA )
  } else {
    ifelse( pa_hmsc[j,i] == 1, initial_peaks_taxon$elev.final.med[i], NA )
  }
}}


# check basic stats on entire assemblage before resampling (sanity check)
peaks.mean <- apply( pa_hmsc, 1, mean, na.rm = T)
peaks.sd <- apply( pa_hmsc, 1, sd, na.rm = T)
boxplot(peaks.sd/peaks.mean)

elev.var.check <- data.frame( select(cover_quads_all, transect, meter.point, prepost), peaks.mean, peaks.sd, elev.var = peaks.sd/peaks.mean)
elev.var.check <- elev.var.check %>% 
  group_by( prepost )%>% 
  mutate( peaks.sd = ifelse( is.na(peaks.sd), 0, peaks.sd),
          elev.var = ifelse( is.na(elev.var), 0, elev.var)) 

elev.var.check %>% 
  summarize( mean.sd = mean(peaks.sd) ,
             mean.mean = mean(peaks.mean) ,
             mean.elev.var = mean(elev.var) ) 

ggplot( elev.var.check, aes(y = elev.var, x = prepost)) + geom_boxplot()

# resampling scheme
# choose the taxa numbers at random for each quadrat. Use this for both the var.partiion call and the elevation variability call
# only select taxa that occurred?
# if not, then stability will be based only on those present
# start by pulling column numbers for each row that are non-zero on not NA
# only use pre-values for this
present.columns.pre <- apply( pa_hmsc[seq( 2,nrow(pa_hmsc), by = 2),], 1, function(z) which( !is.na(z) ) )
sort( unlist(lapply( present.columns.pre, length )) )
present.columns.pre

####----------------------------------------------

#### sampling through a for loop
elev.stab <- list()
for( i in 1:100 ){


picks <- lapply( present.columns.pre, sample, size = 2, replace = F )
# repeat list elements twice to pick the same taxa for pre and post heatwave
picks <- rep( picks, each = 2 )
picksdf <- as.data.frame( do.call(rbind, picks) )
colnames(picksdf) <- c("x1","x2")
picks.split <- split( picksdf, f = transect )
  
elev_list <- split(pa_hmsc, seq(nrow(pa_hmsc)))
elev_picks <- mapply( function(x,y) x[y], elev_list, picks )
ddmean.transect$elev.var.pick <- apply( elev_picks, 2, function(z){
  sd(z,na.rm=T)/mean(z,na.rm=T) 
})
# if NA produced, it means one of the two species was absent, so make that zero for no variability
ddmean.transect$elev.var.pick[ is.na(ddmean.transect$elev.var.pick) ] <- 0
elev.var.transect.pre <- ddmean.transect %>% 
  filter( prepost == "pre" ) %>% 
  group_by( transect) %>% 
  summarize( elev.var.pick.pre = mean(elev.var.pick)) %>% 
  separate( transect, c("Site","Blah","Zone"), remove = F) %>% 
  unite( "Site", Site:Blah, sep = " ")
elev.var.transect.post <- ddmean.transect %>% 
  filter( prepost == "post" ) %>% 
  group_by( transect) %>% 
  summarize( elev.var.pick.post = mean(elev.var.pick)) %>%  # lots of NAs bc taxa were seleced based on "pre" conditions, but identities did change in some cases
  select( elev.var.pick.post )
  
# 
array.list <- lapply( dmt.split, function(z){
  tmp <- z %>% 
    unite( "quadrat", transect, meter.point)
  test.comm <- tmp[, -c(1,2)]
  comm.picks <- test.comm[, ]
  test.split <- split( test.comm, f= factor(tmp$quadrat) )
  test.array <- ld2a( test.split )
  return( test.array )
} )

array.list <- mapply( function(x,y){
  tmp <- x %>% 
    unite( "quadrat", transect, meter.point)
  test.comm <- tmp[, -c(1,2)]
  xsplit <- split(test.comm, seq(nrow(test.comm)))
  ysplit <- split(y, seq(nrow(y)))
  comm.picks <- t(mapply( function(z,w) z[ c(w$x1,w$x2)], xsplit, ysplit ))
  test.split <- split( as.data.frame(comm.picks), f = factor(tmp$quadrat) )
  test.array <- ld2a( test.split )
  return( test.array )
}, dmt.split, picks.split  )


# var.partition accepts the array so
# provide arrays to var.partition in another lapply call
var.synch.quads.in.transects <- lapply( array.list, function(z) {
  var.partition(z)
})

var.synch.quads.in.transects <- do.call( rbind, var.synch.quads.in.transects )

elev.stab[[i]] <- cbind( elev.var.transect.pre, elev.var.transect.post, var.synch.quads.in.transects)

}


elev.stab

elev.stab.bind <- data.table::rbindlist( elev.stab, idcol = T )
elev.stab.bind$Zone <- factor( elev.stab.bind$Zone, levels = c("LOW","MID","HIGH"))
elev.stab.bind$Site <- factor( elev.stab.bind$Site, levels = c("Foggy Cove", "Fifth Beach", "North Beach" ))

mmean <- mclean %>% 
  group_by( transect ) %>% 
  summarize( mean.elev = mean(Shore_height_cm) )
  
elev.stab.join <- left_join( elev.stab.bind, mmean )

ggplot( elev.stab.bind, aes(x = elev.var.pick.pre)) + 
  geom_histogram() +
  facet_wrap(  Site ~ Zone )
ggplot( elev.stab.bind, aes(x = elev.var.pick.post)) + 
  geom_histogram() +
  facet_wrap(  Site ~ Zone )

elev.var.all <- elev.stab.join %>% 
  select( transect, Site, Zone, mean.elev, elev.var.pick.pre, elev.var.pick.post ) %>% 
  pivot_longer( elev.var.pick.pre:elev.var.pick.post, names_to = "prepost", values_to = "elev.var" ) %>% 
  mutate( prepost = ifelse( prepost == "elev.var.pick.pre", "pre", "post")) %>% 
  mutate( prepost = factor(prepost, levels = c('pre','post')))

ggplot( elev.var.all, aes(x = elev.var, fill = prepost)) + 
  geom_histogram( alpha = 0.75, position = "identity" ) +
  facet_wrap(  Site ~ Zone )

# ddmean.transect <- left_join(ddmean.transect, select(mclean, transect, Shore_height_cm) )
# ggplot( ddmean.transect, aes(x = Shore_height_cm, y = elev.var.pick, col = prepost)) + 
#   geom_point(alpha = 0.5) + geom_smooth( method = 'lm')

qs <- seq(0.05, 0.95, by = 0.1)
ggplot( elev.var.all, aes(x = mean.elev, y = elev.var, shape = Site, fill = Zone ) ) +
  facet_wrap( ~prepost ) +
  geom_point( alpha = 0.5 ) + 
  geom_quantile(aes(group = 1), quantiles = qs) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  xlab( "Mean transect elevation (cm)" ) + 
  ylab( "Variability in elevation optima" ) + 
    theme_classic()
ggsave( "R/Figs/elev_var_resample_prepost.svg", width = 6, height = 3.5)

rqfit1 <- rq( elev.var.pick.pre ~ mean.elev, tau = qs, data = elev.stab.join)
summary(rqfit1)
ggplot( filter(elev.var.all, prepost == 'pre'), aes(x = mean.elev, y = elev.var, shape = Site, fill = Zone ) ) +
  geom_point( alpha = 0.5 ) + 
  geom_quantile(aes(group = 1), quantiles = qs) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  xlab( "Mean transect elevation (cm)" ) + 
  ylab( "Variability in elevation optima" ) + 
    theme_classic()
ggsave( "R/Figs/elev_var_resample.svg", width = 4, height = 3.5)

rqfit2 <- rq( phi_S2C_L ~ mean.elev, tau = qs, data = elev.stab.join)
summary(rqfit2)
ggplot( elev.stab.join, aes(x = mean.elev, y = phi_S2C_L, shape = Site, fill = Zone ) ) +
  # facet_wrap( ~prepost ) +
  geom_point( alpha = 0.5 ) + 
  # geom_smooth( aes(group = 1), show.legend = FALSE, col = 'black', se = F ) +
  geom_quantile(aes(group = 1), quantiles = qs) +
  scale_shape_manual( values=c(21,22,24)) +
  scale_fill_manual( values=c("black","gray50","whitesmoke") ) +
  xlab( "Mean transect elevation (cm)" ) + 
  ylab( "Local-scale species synchrony" ) + 
  theme_classic()
ggsave( "R/Figs/synchrony_elevation_resample.svg", width = 4, height = 3.5)



ggplot( elev.stab.bind, aes(x = CV_S_L)) + 
  geom_histogram() +
  facet_wrap(  Site ~ Zone )
ggplot( elev.stab.bind, aes(x = phi_S2C_L)) + 
  geom_histogram() +
  facet_wrap(  Site ~ Zone )

ggplot( elev.stab.bind, aes(x = Zone, y = phi_S2C_L, col = Site)) + 
  geom_boxplot( notch = F, fill = 'whitesmoke' ) +
  ylab( expression(atop( "Local-scale species",paste("synchrony (",italic(phi[S %->% "C,L"]),")" ))) ) +
  theme_classic()
ggsave( "R/Figs/synchrony_local_resample_boxplot.svg", width = 4.33, height = 3 )


ggplot( elev.stab.bind, aes(x = phi_S_L2R)) + 
  geom_histogram() +
  facet_wrap(  Site ~ Zone )


ggplot( elev.stab.join, aes(x = mean.elev, y = CV_S_L ) ) +
  geom_point() + geom_smooth( method = 'lm')
ggplot( elev.stab.join, aes(x = mean.elev, y = CV_C_L ) ) +
  geom_point() + geom_smooth()
ggplot( elev.stab.join, aes(x = mean.elev, y = CV_S_R ) ) +
  geom_point() + geom_smooth()
ggplot( elev.stab.join, aes(x = mean.elev, y = CV_C_R ) ) +
  geom_point() + geom_smooth()
ggplot( elev.stab.join, aes(x = mean.elev, y = CV_C_R ) ) +
  geom_point() + geom_smooth()
ggplot( elev.stab.join, aes(x = mean.elev, y = phi_S_L2R ) ) +
  geom_point() + geom_smooth()
ggplot( elev.stab.join, aes(x = mean.elev, y = phi_C_L2R ) ) +
  geom_point() + geom_smooth()
ggplot( elev.stab.join, aes(x = mean.elev, y = phi_S2C_L ) ) +
  # geom_smooth( aes(col=Site), method = 'lm' ) +
  geom_smooth( method = 'lm' ) +
  geom_point() 
ggplot( elev.stab.join, aes(x = mean.elev, y = phi_S2C_R ) ) +
  geom_point() + geom_smooth()
psych::pairs.panels( select( elev.stab.join, CV_S_L, CV_C_L, CV_S_R, CV_C_R) )
summary( lm( CV_S_L ~ elev.var.pick.pre, elev.stab.join ) )
a <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = CV_S_L ) ) +
  geom_point(alpha = 0.5) + #geom_smooth( method='lm' ) +
  geom_quantile(aes(group = 1), quantiles = qs) +
  ylab( expression(atop( "Local-scale species",paste("variability (",italic(CV["S,L"]),")" ))) ) +
  xlab( "Elevation optima variability (CV)" ) +
  theme_classic()
summary( lm( CV_C_L ~ elev.var.pick.pre, elev.stab.join ) )
b <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = CV_C_L ) ) +
  geom_point(alpha = 0.5) + #geom_smooth( method='lm' )+
  geom_quantile(aes(group = 1), quantiles = qs) +
  ylab( expression(atop( "Local-scale community",paste("variability (",italic(CV["C,L"]),")" ))) ) +
  xlab( "Elevation optima variability (CV)" ) +
  theme_classic()
summary( lm( CV_S_R ~ elev.var.pick.pre, elev.stab.join ) )
c <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = CV_S_R ) ) +
  geom_point(alpha = 0.5) + #geom_smooth( method='lm' )+
  geom_quantile(aes(group = 1), quantiles = qs) +
  ylab( expression(atop( "Metapopulation",paste("variability (",italic(CV["S,R"]),")" ))) ) +
  xlab( "Elevation optima variability (CV)" ) +
  theme_classic()
summary( lm( CV_C_R ~ elev.var.pick.pre, elev.stab.join ) )
d <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = CV_C_R ) ) +
  geom_point(alpha = 0.5) + #geom_smooth( method='lm' )+
  geom_quantile(aes(group = 1), quantiles = qs) +
  ylab( expression(atop( "Metacommunity",paste("variability (",italic(CV["C,R"]),")" ))) ) +
  xlab( "Elevation optima variability (CV)" ) +
  theme_classic()
cowplot::plot_grid( a,b,c,d, nrow = 2)
ggsave( "R/Figs/stability_elevvar_resample.svg", width = 6, height = 5)

a <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = phi_S_L2R ) ) +
  geom_point() + geom_smooth( method='lm' )
b <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = phi_C_L2R ) ) +
  geom_point() + geom_smooth( method='lm' )
c <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = phi_S2C_L ) ) +
  geom_point() + geom_smooth( method='lm' )
d <- ggplot( elev.stab.join, aes(x = elev.var.pick.pre, y = phi_S2C_R ) ) +
  geom_point() + geom_smooth( method='lm' )
cowplot::plot_grid( a,b,c,d, nrow = 2)

ggplot( elev.stab.join, aes(x = phi_S2C_L, y = CV_C_L ) ) +
  geom_point() + geom_smooth( method='lm' )
