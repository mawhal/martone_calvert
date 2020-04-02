# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses community trajectory analysis on Martone seaweed community data

# load libraries
library( tidyverse )
library( betapart )
library( vegan )
library( adespatial )
library( ggrepel )
library( RColorBrewer )
#


## define color scheme for temperature anomalies
# from the plot of seasonal anomalies the range of SST anomalies is -2 to 2
anom.range <- c(-2,2)
n=9
cols <- brewer.pal(n,"RdBu")
data.frame(
  cols,
  anom=seq(anom.range[1],anom.range[2],leng=n) )
# the range of anomalies is more like -0.5 to 1.5
pal <- colorRampPalette(rev(cols)[4:8])
pal2 <- colorRampPalette(rev(cols))

## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )


## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
muse <- am[ am$Year != "2011", ]
# remove Meay Channel
## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
muse <- muse[ muse$Site != "Meay Channel", ]
# Only use Mid-shore transects for now
# muse <- muse[ muse$Zone == "MID", ]
# muse <- muse[ muse$Site == "North Beach", ]
muse <- droplevels(muse)

# restrict rows of d to ones with UID's in the metadata
duse <- ad[ ad$UID %in% muse$UID, ]
dm <- left_join( duse, muse )


# for now, restrict community analysis to algae only
d <- dm %>% 
  # filter( non.alga.flag %in% c("Algae"  )
  filter( motile_sessile == "sessile" )

# average elevation per zone
dmeanelev <- d %>% 
  group_by( Year, Site, Zone ) %>%
  summarise( Elevation=mean(Shore_height_cm, na.rm=T) )
  
# calculate mean abundance per transect in each year
dmean <- d %>% 
  group_by( Year, Site, Zone, taxon_lumped ) %>%
  summarise( Abundance=mean(Abundance) )

# spread out
d.comm.mean <- dmean %>%
  spread( taxon_lumped, Abundance, fill=0 ) %>% 
  ungroup() %>% 
  mutate( Zone = factor(Zone, levels=c("LOW","MID","HIGH")) ) %>% 
  arrange( Year, Site, Zone )

d.comm.mean <- d.comm.mean %>% filter( !is.na(Elevation) )

# remove UID column from community data
meta <- d.comm.mean[ ,1:3 ]
meta <- left_join(meta,dmeanelev)
comm <- as.matrix(d.comm.mean[,-c(1:3)])

# # define site as the particular trasect
# Zone <- factor( meta$Zone, levels=c("LOW","MID","HIGH"), labels=c("Low","Mid","High") )
# Site <- factor( meta$Site, labels=c("5","N","W") )
# site <- paste( Site, Zone, sep="." )
# site <- factor( site, levels= c("5.Low","5.Mid","5.High","N.Low","N.Mid","N.High","W.Low","W.Mid","W.High"),
#                 labels= c("5.L","5.M","5.H","N.L","N.M","N.H","W.L","W.M","W.H"))
# 
# # define year
# year <- meta$Year
# # define 2016 for plotting later
# year2 <- as.character(year)
# year2[year2!=2016] <- ""
# 

## temperature anomaly data from Pine Island
anoms <-  read_csv("Data/environmetal_data/Lighthouse Data/through May 2019_Peter Chandler/output from R/PineIsland_monthly_SST_anomaly.csv")

# calculate winter and summer temperature anomalies
library(zoo)
yq <- as.yearqtr( as.yearmon( paste(anoms$month,anoms$year,sep="/"), "%m/%Y") + 1/12)
anoms$season <- factor(format(yq, "%q"), levels = 1:4, 
                       labels = c("winter", "spring", "summer", "fall"))
anoms.season <- anoms %>% 
  group_by( year, season ) %>% 
  summarize( temp.anom=mean(temp.anom) )
ggplot(anoms.season, aes(x=year,y=temp.anom,col=season)) + geom_point()
ggplot(anoms.season, aes(x=year,y=temp.anom,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
ggplot(filter(anoms.season, year>=2010), aes(x=year,y=temp.anom,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
# figure out how to make an anomly plot with vertical lines from zero

# extract data for 2011 to 2019
as.survey <- anoms.season %>% 
  filter( year>=2010 ) %>% 
  spread( key = season, value=temp.anom )


# which species are exceedingly rare?
sort(colSums(comm))
hist(colSums(comm),breaks = seq(0,1600,2))
# which species are most abundant through the time series?
dominance <- apply( comm, 1, function(z) names(z)[order(z,decreasing = T)]  )
dominance[1:3,]

# quick beta diversity by group
z <- betadiver(comm,"z")
mod <- betadisper(z, meta$Site)
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



### Comparing transects at different dates
## betapart temporal
# need two matrices of same dimensions with data from two different years
commpa <- comm 
commpa[commpa>0] <- 1
d12 <- comm[ meta$Year == 2012, ]
d19 <- comm[ meta$Year == 2019, ]


bpt <- bind_cols(meta[ meta$Year == 2012, ], beta.temp( d12, d19 ) )
bpt %>% mutate( turn=beta.sim/beta.sor, nest=beta.sne/beta.sor )

# Legendre method (TBI)
tbia <- TBI( d12, d19, method = "%difference" )
plot(tbia)
tbip <- TBI( d12, d19, method = "%difference", pa.tr = TRUE )
plot(tbip)

# get TBI matrix for all comparisons to 2012
tbi.sums.abun <- matrix( nrow = length(unique(meta$Year))-1, ncol=5 )
tbi.sums.pa   <- matrix( nrow = length(unique(meta$Year))-1, ncol=5 )
for( i in 2:length(unique(meta$Year)) ){
  di <- comm[ meta$Year == unique(meta$Year)[i], ]
  tbi.sums.abun[i-1,] = unlist(TBI( d12, di, method = "%difference" )$BCD.summary[1:5])
  tbi.sums.pa[i-1,]   = unlist(TBI( d12, di, method = "%difference", pa.tr = TRUE )$BCD.summary[1:5])
}

# add row and column names
tbi.sums.abun <- as.data.frame(tbi.sums.abun)
names(tbi.sums.abun) <- c("losses","gains","total","B%","C%")
tbi.sums.abun$comparison <- paste("2012",unique(meta$Year)[-1],sep = "-")
tbi.sums.abun$comp <- 1:7
tbi.sums.pa <- as.data.frame(tbi.sums.pa)
names(tbi.sums.pa) <- c("losses","gains","total","B%","C%")
tbi.sums.pa$comparison <- paste("2012",unique(meta$Year)[-1],sep = "-")
tbi.sums.pa$comp <- 1:7

# gather
gath.abun <- tbi.sums.abun %>% select(losses,gains,total,comparison,comp) %>% 
  group_by(comparison,comp) %>% 
  gather(key,value,losses:total)
gath.pa <- tbi.sums.pa %>% select(losses,gains,total,comparison,comp) %>% 
  group_by(comparison,comp) %>% 
  gather(key,value,losses:total)

a <- ggplot( gath.abun, aes(x=comp,y=value,col=key,shape=key)) + 
  geom_line() + geom_point(size=2) + 
  ylim(c(0,0.75)) + scale_color_manual(values=c("red","blue","black")) +
  ylab("Dissimilarity") + xlab("Survey year pairs" ) +
  scale_x_continuous( breaks=1:7,labels=tbi.sums.abun$comparison ) +
  theme_classic()
b <- ggplot( gath.pa, aes(x=comp,y=value,col=key,shape=key)) + 
  geom_line() + geom_point(size=2) + 
  ylim(c(0,0.5)) + scale_color_manual(values=c("red","blue","black")) +
  ylab("Dissimilarity") + xlab("Survey year pairs" ) +
  scale_x_continuous( breaks=1:7,labels=tbi.sums.abun$comparison ) +
  theme_classic()

cowplot::plot_grid( a,b, ncol=1, rel_heights = c(2,1.5) )



### consider dbRDA with factor for year and site (and transect?)
meta$year <- factor(meta$Year, ordered=F)
meta$site <- factor(meta$Site, ordered=F)
db1 <- dbrda( comm~year+Elevation, data=meta )
plot(db1)
os1 <- ordisurf( db1, meta$Elevation )
ordiellipse( db1, meta$Zone, conf=0.6 )


RsquareAdj(db1)  # explains 20-25% of variation in consumption rate?
R2 <- eigenvals(db1)/sum(eigenvals(db1))
R2
summary(db1)
cummr2 <- R2
for(i in 2:length(eigenvals(db1))){
  cummr2[i] <- cummr2[i]+cummr2[i-1]  
}

# extract axes
nax <- 1:2
sites     <- data.frame(scores(db1,choices = nax, scaling = 0)$sites)
sites$zone <- factor( meta$Zone, levels=c("LOW","MID","HIGH") )
sites$year <- meta$Year
centroids <- data.frame(scores(db1,choices = nax, scaling=0)$centroids)
centroids$year <- 2012:2019

ordi.grid <- os1$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.na <- data.frame(na.omit(ordi)) #gets rid of the nas
ordi.na #looks ready for plotting!


### compare centroids to temperature anomaly data
centroids
as.survey$summer1 <- c(NA, as.survey$summer[1:length(as.survey$summer)-1])
surv.cent <- left_join( as.survey, centroids )
centroid2 <- centroids
centroid2$year <- centroid2$year-1
names(centroid2) <- c("shift1","shift2","year")
surv.cent <- left_join( surv.cent, centroid2  )
surv.cent$year2 <- paste0("'",substr(surv.cent$year+1,start=3,stop=4))

ggplot( sites, aes(x=dbRDA1,y=dbRDA2)) + 
  # stat_contour(data = ordi.na, aes(x = x, y = y, z = z, colour = rev(..level..)),
  #              binwidth = 20)+ #can change the binwidth depending on how many contours you want
  geom_smooth( aes(group=year), method='lm', se=F, size=0.5,col='gray85') +
  geom_point( data=sites, aes(shape=zone), size=1 ) +
  # geom_smooth( aes(group=1), method='lm', se=F, size=0.5,col='gray25') +
  geom_path( data=surv.cent, aes(col=summer1),size=1 ) +
  geom_point( data=surv.cent, aes(fill=summer1), size=3, shape=21 ) + 
  geom_text_repel( data=surv.cent, label=surv.cent$year, 
                   fontface="bold",
                   point.padding = 0.3, box.padding = 0.1 ) +
  scale_shape_manual(values = c(1,0,2)) +
  xlab(paste0(names(sites)[1],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(sites)[2],' (',round(R2[2],3)*100, '%)')) +
  scale_fill_gradientn(colours=pal2(100),limits=c(-2,2)) +
  scale_color_gradientn(colours=pal2(100),limits=c(-2,2)) +
  theme_bw() + theme( panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
ggsave("R Code and Analysis/Figs/rda_bwr_pal2.svg", width=4, height=3 )





windows(2,2)
ggplot( surv.cent, aes(y=dbRDA2,x=winter) ) + geom_point()
ggplot( surv.cent, aes(y=shift2,x=winter) ) + geom_point()
ggplot( surv.cent, aes(y=shift2,x=summer) ) + 
  # geom_smooth(method="lm",se=F) + 
  geom_path() +
  geom_point(size=2) + 
  geom_text_repel(label=surv.cent$year2, 
                  point.padding = 0.1,box.padding = 0.3,
                  direction="both",
                  nudge_y = c(-0.02, 0.08, -0.02, 0.01,
                              -0.01, 0.01, -0.01, 0.01)) +
  ylab("centroids | elevation") +
  xlab(expression(paste(degree,"C (summer-1)"))) +
  theme_classic()
ggsave( "R Code and Analysis/Figs/rda_summer_centroid.svg", width=2, height=2 )
ggplot( surv.cent, aes(y=dbRDA2,x=winter) ) + 
  geom_smooth(method="lm") + geom_point() + 
  geom_text_repel(label=surv.cent$year) 




# show summer versus winter temperature anomaly
as.survey.all <- anoms.season %>% 
  # filter( year>=2010 ) %>% 
  spread( key = season, value=temp.anom )

ggplot( as.survey, aes(y=summer,x=winter) ) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_smooth(se=T,size=0.5) +
  geom_point(size=2) +  geom_path(size=1) +
  geom_text_repel(label=surv.cent$year, box.padding = 0.5) +
  ylab(expression(paste("Summer SST anomaly (",degree,"C)"))) +
  xlab(expression(paste("Winter SST anomaly (",degree,"C)"))) +
  theme_bw()

ggplot( as.survey.all, aes(y=summer,x=winter) ) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_smooth(se=T,size=0.5) +
  geom_point(size=1,shape=1) +
  geom_point(data=as.survey,size=2) +  geom_path(data=as.survey,size=1) +
  geom_text_repel(data=as.survey,label=as.survey$year, box.padding = 0.5) +
  ylab(expression(paste("Summer anomaly (",degree,"C)"))) +
  xlab(expression(paste("Winter anomaly (",degree,"C)"))) +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
ggsave( "R Code and Analysis/Figs/winter_summer_SST.svg", width=2, height=2 )
