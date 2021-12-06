# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen 
# created 2 Dec 2021

# Using Hammond et al. 2020 framework for stability and synchrony
# https://onlinelibrary.wiley.com/doi/abs/10.1002/ecs2.3078

# 

# 
library(tidyverse)



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
# d.select <- d.select %>% filter( Year %in% c(2012,2013,2018,2019))



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
  select(-UID, -Year, -Site, -Zone, -Quadrat, -Meter.point, -post, -pre, -both, -Shore_height_cm) %>% 
  group_by(quadrat, prepost) %>% 
  summarize_all( .funs = mean)


# remove taxa that do not appear
colSums(ddmean[,-c(1,2)])


# pivot_longer to get taxa and cover columns
ddmeanlong <- ddmean %>% 
  pivot_longer( cols = Acrosiphonia:Unknown.crust, names_to = "taxon", values_to = "cover" ) 


### Quadrat - level
### calculations  necessary for each measure

# only take data from before and after the heatwaves started
# d.select <- d.select %>% filter( Year %in% c(2012,2013,2018,2019))


# temporal mean, temporal standard deviation of individual populations
ddcalcs.pop <- ddmeanlong %>% 
  group_by( quadrat, taxon ) %>% 
  mutate( sd.pop = sd(cover), mean.pop = mean(cover)) %>% 
  mutate( deviation = cover-mean.pop )

# test covariance
x1 <- ddcalcs.pop$cover[ ddcalcs.pop$quadrat == "Fifth Beach HIGH 10" & ddcalcs.pop$taxon == "Acrosiphonia"]
x2 <- ddcalcs.pop$cover[ ddcalcs.pop$quadrat == "North Beach MID 8" & ddcalcs.pop$taxon == "Ulva"]

cov(x1,x2)

# with two time points can simplify
# just select pre or post, the use deviations among 

# covariances
# first filter out species that never occurred (i.e. temporal mean cover is 0)?
# or maybe not?

## FOR LOOP
# cov.pop <- array( dim = 4)
# for(i in 1:length(unique(ddmeanlong$taxon)) ){
#   for(j in 1:length(unique(ddmeanlong$taxon))){
#     for(k in 1:length(unique(ddmeanlong$quadrat))){
#       for(l in 1:length(unique(ddmeanlong$quadrat))){
#         cov.pop[i][j][k][l] = cov( ddmeanlong$cover[ddmeanlong$taxon==unique(ddmeanlong$taxon)[i] & ddmeanlong$quadrat==unique(ddmeanlong$quadrat)[k]],
#                             ddmeanlong$cover[ddmeanlong$taxon==unique(ddmeanlong$taxon)[j] & ddmeanlong$quadrat==unique(ddmeanlong$quadrat)[l]] )
#   }
#     }}}
# this for loop takes forever

# could create a matrix of community by taxon, where each cell contains all of the numbers needed


# try to make a matrix of populations by 


# temporal mean of metacommunity cover (total across all quads)
ddmeanlong %>% 
  ungroup() %>% 
  summarize( mean.meta = mean(cover) ) 











### measures
cv1 <- ddcalcs %>% 
  ungroup() %>% 
  summarize( sum.sd.pop = sum(sd.pop), sum.mean.pop = mean(mean.pop) ) %>% 
  mutate( cv1 = (sum.sd.pop/sum.mean.pop)^2 )








# Transect-level approaches
# calculate things for each transect, or across the metacommunity using transect as the focal level
# split data.frame by transect, then look at variability and asynchony

# Start with sites as metacommunities

# for covariance matrix -- make a species by time point matrix, then carry on with calculations
# 


# pivot_longer to get taxa and cover columns
ddlong <- dd %>% 
  pivot_longer( cols = Acrosiphonia:Unknown.crust, names_to = "taxon", values_to = "cover" ) 

# sum all covers across transects
ddtransect_sum <- ddlong %>% 
  group_by(transect, Year, taxon) %>% 
  summarize( cover.sum = sum(cover) )

### -- M = temporal mean metacommunity cover (total across all quadrats, transects, and sites)
M.transect <- ddtransect_sum %>% 
  group_by( Year ) %>% 
  summarize( meta.cover = sum( cover.sum) ) %>% 
  ungroup() %>% 
  summarize( M = mean(meta.cover) )

### -- temporal mean, temporal standard deviation of individual populations
ddcalcs.pop <- ddtransect_sum %>% 
  group_by( transect, taxon ) %>% 
  summarize( sd.pop = sd(cover.sum), mean.pop = mean(cover.sum))

cv1.transect <- ddcalcs.pop %>% 
  ungroup() %>% 
  summarize( cv1 = (sum(sd.pop)/mean(mean.pop))^2 )

# pivot wider
ddtswide <- ddtransect_sum %>% 
  pivot_wider( names_from = "taxon", values_from = cover.sum)


# separate by transect
test <- ddtswide %>% filter( transect == "Fifth Beach_HIGH")
testcom <- test[,-c(1:2)]


sum.var.cov.diff <- function(){
  
}
  
# variance-covariance matrix
testcov <- cov(testcom)
# variances in the diagonal
tmp <- diag(testcov)
var.prod <- combn( tmp, m =2, FUN = prod, simply = T )
var.prod.mat <- matrix( NA, nrow = length(tmp), ncol = length(tmp) )
var.prod.mat[ upper.tri(var.prod.mat) ] <- var.prod

# covariances - how are these arranged?
# cov.mat <- testcov[upper.tri(testcov)]

# calculate the numerator as the sum of differences in products of variances minus covariances
# apply()
tmp.mat <- matrix( NA, nrow = length(tmp), ncol = length(tmp) )
for( i in 1:(length(tmp)-1) ){
  for( j in i+1:(length(tmp)-i) ){
    tmp.mat[i,j] <- var.prod.mat[i,j] - testcov[i,j]
  }
}
sum( tmp.mat, na.rm = T )


###  -- del = stabilization from asynchrony among species in local communities (type II)



