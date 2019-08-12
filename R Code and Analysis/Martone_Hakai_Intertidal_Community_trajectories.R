# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses community trajectory analysis to Martone seaweed community data

# load libraries
library( tidyverse )
library( vegclust )   # CTA, length, angles, directionality, projection, distances between segments/trajectories
library( adespatial ) # dynamic-based beta diversity, local contributions to beta diversity
library( smacof )     # stress minimization using majorization (smacof) -- new approaches for MDS
library( vegan )
library( RColorBrewer )

# To specify community dynamics, we need three data items:
#   
# A set of community states (i.e. coordinates in a space Ω), described using a distance matrix d;
# A vector specifying the site (i.e. sampling unit) corresponding to each community state;
# A vector specifying the survey (i.e. time point) corresponding to the sampling of each community state.
# 

## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read.csv( "Data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv", stringsAsFactors = FALSE )
# all metadata
am <- read.csv( "Data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv", stringsAsFactors = TRUE )





## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
muse <- am[ am$Year != "2011", ]
# remove Meay Channel
## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
muse <- muse[ muse$Site != "Meay Channel", ]
# Only use Mid-shore transects for now
muse <- muse[ muse$Zone == "MID", ]
# muse <- muse[ muse$Site == "North Beach", ]
muse <- droplevels(muse)


# restrict rows of d to ones with UID's in the metadata
duse <- ad[ ad$UID %in% muse$UID, ]
dm <- left_join( duse, muse )


# for now, restrict community analysis to algae only
d <- dm %>% 
  filter( motile_sessile != "motile" )
  # filter( non.alga.flag =="Algae" )


# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which will reduce the size of the dataset a bit
# restrict this to seaweeds and sessile invertebrates
d.simple <- d %>%
  group_by( UID, Year, Site, Zone, taxon_lumped2 ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))

# average cover per transect
dmean <- d.simple %>% 
  spread( taxon_lumped2, Abundance, fill=0 ) %>%
  gather( taxon_lumped2, Abundance, -UID, -Year, -Site, -Zone ) %>%
  group_by( Year, Site, Zone, taxon_lumped2 ) %>%
  summarise( Abundance=mean(Abundance) )

# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- dmean %>%
  spread( taxon_lumped2, Abundance, fill=0 )

# order by site and zone
d.comm <- d.comm %>%
  arrange( Site, factor(Zone,levels=c("LOW","MID","HIGH")) )
d.comm$Zone <- factor( d.comm$Zone,levels=c("LOW","MID","HIGH") )


# isolate the community, site, and sample data
comm <- d.comm[,-c(1:3)]
site <- paste( d.comm$Site, d.comm$Zone, sep="." )
tran <- apply( d.comm[,c(2,3)],1,paste, collapse="." )
year <- d.comm$Year
year2 <- as.character(year)
year2[year2!=2016] <- ""
yearmod <- factor(year,ordered=T)
zone <- d.comm$Zone

# what is the distribution of values like?
k <- sample(nrow(comm),1)
cs <- ceiling(comm[k,])
(x<- fisherfit( cs ))
(x<- prestonfit( unlist(cs) ))
plot(x)
histogram(unlist(cs))
x<-radfit(cs)
plot(x)

# which species are exceedingly rare?
sort(colSums(comm))
hist(colSums(comm),breaks = seq(0,1200,2))
# which species are most abundant through the time series?
dominance <- apply( comm, 1, function(z) names(z)[order(z,decreasing = T)]  )
dominance[1:3,]

# quick beta diversity by group
z <- betadiver(comm,"z")
mod <- betadisper(z, site)
plot(mod)
plot(mod, axes=c(3,1) )
boxplot(mod)
anova(mod)
(mod3B <- betadisper(z, site, type = "median", bias.adjust=TRUE))
anova(mod3B)
permutest(mod3B, permutations = 99)
plot(mod3B)
boxplot(mod3B)

adonis2( comm~site*year, by='margin'  )

# Trajectory analysis
# calculate a distance 
D_man <- vegdist( comm, method="manhattan", transform = function(x){log(x+1)})
D_bray <- vegdist( comm, method="bray" )
D_use <- D_bray

# custom color scheme (three shades of three colors)
lighten <- function(color, factor=1.5){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}
darken <- function(color, factor=2){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}
base <- c('#1b9e77', '#7570b3', '#d95f02')
cols = as.vector(matrix( c(sapply( base, darken ), sapply( base, lighten ), base), ncol=3, byrow = TRUE ))

# display trjectories in PCoA
par(mar=c(5,4,1,1)+0.1)
x <- trajectoryPCoA( D_use,  as.numeric(factor(site)), year,
                     traj.colors = cols, 
                     axes=c(1,2), length=0.1, lwd=2 )
text( x$points[,1:2],labels = year2, pos = 1, col=rep(cols,each=length(unique(year))) )
legend("topleft", bty="n", legend = unique(site), col = cols, lwd=3 )

# centered trajectories
x<-trajectoryPCoA( centerTrajectories(D_use, as.numeric(factor(site)) ),  as.numeric(factor(site)), year,
                traj.colors = cols, 
                axes=c(1,2), length=0.1, lwd=2 )
text( x$points[,1:2],labels = year2, pos = 1, col=rep(cols,each=length(unique(year))) )
legend("topright", bty="n", legend = unique(site), col = cols, lwd=2)

# can also use MDS to represent trajectories
# stress plot

mMDS  <- mds( D_use, ndim=8 )
mMDS
mMDS2 <- mds( centerTrajectories(D_man, as.numeric(factor(site)) ) ) 
mMDS2
use <- mMDS$conf
trajectoryPlot( use,  as.numeric(factor(site)), year,
                traj.colors = cols, 
                axes=c(1,2), length=0.1, lwd=2 )
text( use,labels = year2, pos = 1, col=rep(cols,each=length(unique(year))) )
legend("topright", bty="n", legend = unique(site), col = cols, lwd=2)




## Trajectory statistics
trajectoryLengths(        D_man, site, year ) # trajectory lengths always greater in the LOW zone for seaweeds
trajectoryAngles(         D_man, site, year )
trajectoryAngles(         D_man, site, year, all=TRUE ) # high degree of angle homogeneity
trajectoryDirectionality( D_man, site, year ) # despite longer segements in LOW, MID and HIGH often are more directional
plot( trajectoryDirectionality( D_man, site, year ), cex=3, pch=21, bg=cols ) # similar directionality among sites
trajectoryProjection(     D_man, 1,2:6 )
trajectoryConvergence(    D_man, site, year, symmetric = FALSE )

Ds = segmentDistances( D_man, site, year )$Dseg
Ds
mMDS = mds(Ds)
mMDS

xret = mMDS$conf
par(mar=c(4,4,1,1))
plot(xret, xlab="axis 1", ylab = "axis 2", asp=1, pch=21,
     bg=rep(cols, each=7), 
     xlim=c(-1.5,1), ylim=c(-1,1.5))
text(xret, labels=rep(paste0("s",1:7),9), pos=1)
legend("topleft", pt.bg=cols, pch=21, bty="n", legend=paste0("trajectory",1:9))


trajectoryDistances( D_man, site, year, distance.type = "Hausdorff")

                    