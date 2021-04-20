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


# choose which to include
d <- dm %>% 
  filter( motile_sessile == "sessile" ) #%>% 
  # filter( non.alga.flag %in% c("Algae")  ) # restrict community analysis to algae only

# take a closer look
d %>% 
  group_by(taxon_lumped2) %>% 
  summarize(total=sum(Abundance)) %>% 
  arrange( total )
# remove a few taxa that are unidentified or extremely low abundance
d %>% 
  filter( (taxon_lumped %in% c("articulated coralline","Unknown red blade","")) )
  
# Haliclona to sponge - fix these below using taxon_lumped, but likely need to update in data prep
# Mytilus trossolus combine with others - using taxon_lumped2 fixes this, but changes many other things, too. changes total taxa from 165 to 141

# average elevation per zone
dmeanelev <- d %>% 
  group_by( Year, Site, Zone ) %>%
  summarise( Elevation=mean(Shore_height_cm, na.rm=T) )
  
# calculate mean abundance per transect in each year
# first fix some taxa names
d$taxon_lumped[d$taxon_lumped=="Haliclona"] <- "Porifera"
d$taxon_lumped[d$taxon_lumped=="Mytilus trossulus"] <- "Mytilus sp."
d$taxon_lumped2[d$taxon_lumped2=="Bossiella_articulate"] <- "Bossiella articulate"
d$taxon_lumped2[d$taxon_lumped2=="Bossiella_crust"] <- "Bossiella crust"
dmean <- d %>% 
  filter( !(taxon_lumped %in% c("articulated coralline","Unknown red blade")) ) %>% 
  group_by( Year, Site, Zone, taxon_lumped2 ) %>%
  summarise( Abundance=mean(Abundance) )

### quick trial of "genus" level
d$genus <- unlist(lapply(strsplit( d$taxon_lumped2, split = " " ), function(z) z[1]))
## family level
tax <- read_csv("Data/taxa/algal_taxonomy.csv")
# animal taxonomy , just replace with "genus" for now because this is already above genus level 
dtax <- left_join( d, tax )
dtax$family[ d$non.alga.flag == "Animal" ] <- dtax$genus[ d$non.alga.flag == "Animal" ]
dtax$class[ d$non.alga.flag == "Animal" ] <- dtax$genus[ d$non.alga.flag == "Animal" ]
dtax$order[ d$non.alga.flag == "Animal" ] <- dtax$genus[ d$non.alga.flag == "Animal" ]

dtax$family[ d$genus == "Petrocelis" ] <- unique(dtax$family[ d$genus == "Mastocarpus" ])
dtax$class[ d$genus == "Petrocelis" ] <- unique(dtax$class[ d$genus == "Mastocarpus" ])
dtax$order[ d$genus == "Petrocelis" ] <- unique(dtax$order[ d$genus == "Mastocarpus" ])
dtax$family[ d$genus == "Ralfsioid" ] <- "Ralfsiaceae"
dtax$order[ d$genus == "Ralfsioid" ] <- "Ralfsiales"
dtax$class[ d$genus == "Ralfsioid" ] <- "Phaeophyceae"
dtax$family[ d$genus == "coralline" ] <- "Corallinaceae"
dtax$order[ d$genus == "coralline" ] <- "Corallinales"
dtax$class[ d$genus == "coralline" ] <- "Florideophyceae"
dtax$family[ d$genus == "Colonial" ] <- "Bacillariophyceae"
dtax$order[ d$genus == "Colonial" ] <- "Bacillariophyceae"
dtax$class[ d$genus == "Colonial" ] <- "Bacillariophyceae"
dtax$family[ d$genus == "Opuntiella" ] <- "Furcellariaceae"
dtax$order[ d$genus == "Opuntiella" ] <- "Gigartinales"
dtax$class[ d$genus == "Opuntiella" ] <- "Florideophyceae"


dmean <- dtax %>% 
  filter( !(taxon_lumped %in% c("articulated coralline","Unknown red blade")) ) %>% 
  group_by( Year, Site, Zone, family ) %>%
  summarise( Abundance=mean(Abundance) )

# spread out
d.comm.mean <- dmean %>%
  filter( !is.na(family)) %>% 
  spread( family, Abundance, fill=0 ) %>% 
  ungroup() %>% 
  mutate( Zone = factor(Zone, levels=c("LOW","MID","HIGH")) ) %>% 
  arrange( Year, Site, Zone )

d.comm.mean <- d.comm.mean %>% filter( !is.na(Elevation) )

# remove UID column from community data
meta <- d.comm.mean[ ,1:3 ]
meta <- left_join(meta,dmeanelev)
comm <- as.matrix(d.comm.mean[,-c(1:3)])

# interrogate the dataset
sort( colSums(comm), decreasing = T )
# should we remove 


# # read temperature data
# pine <- read_csv( "R Code and Analysis/output from r/PineIsland_summary.csv" )
# # merge temperature and the rest of the metadata
# M <- left_join( meta, pine )
# # define site as the particular trasect
# Zone <- factor( meta$Zone, levels=c("LOW","MID","HIGH"), labels=c("Low","Mid","High") )
# Site <- factor( meta$Site, labels=c("5","N","W") )
# site <- paste( Site, Zone, sep="." )
# site <- factor( site, levels= c("5.Low","5.Mid","5.High","N.Low","N.Mid","N.High","W.Low","W.Mid","W.High"),
#                 labels= c("5.L","5.M","5.H","N.L","N.M","N.H","W.L","W.M","W.H"))
# # define year
# year <- meta$Year
# # define 2016 for plotting later
# year2 <- as.character(year)
# year2[year2!=2016] <- ""
# 

## temperature anomaly data from Pine Island
anoms <-  read_csv("Data/R code for Data Prep/output from R/PineIsland_monthly_SST_anomaly.csv")

# calculate winter and summer temperature anomalies
library(zoo)
yq <- as.yearqtr( as.yearmon( paste(anoms$month,anoms$year,sep="/"), "%m/%Y") + 1/12)
anoms$season <- factor(format(yq, "%q"), levels = 1:4, 
                       labels = c("winter", "spring", "summer", "fall"))
anoms.season <- anoms %>% 
  group_by( year, season ) %>% 
  summarize( temp.anom=mean(temp.anom) )
ggplot(anoms.season, aes(x=year,y=temp.anom,col=season)) +  geom_line()
ggplot(anoms.season, aes(x=year,y=temp.anom,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
ggplot(filter(anoms.season, year>=2010), aes(x=year,y=temp.anom,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
# figure out how to make an anomly plot with vertical lines from zero
# pairwise correlations among seasonal anomalies
as.all <- anoms.season %>% 
  spread( key = season, value=temp.anom )
winter = ts(as.all$winter)
spring = ts(as.all$spring)
summer = ts(as.all$summer)
fall   = ts(as.all$fall)
ccf(winter, summer)
ccf(spring, summer)
ccf(summer, fall)
ccf(winter, fall)
cor(winter, summer) # 0.77
cor(spring, summer) # 0.88
cor(winter, fall)   # 0.69
cor(fall, summer)   # 0.73

# extract data for 2011 to 2019
as.survey <- anoms.season %>% 
  filter( year>=2010 ) %>% 
  spread( key = season, value=temp.anom )



test <- c("Station_1_CTD_42","Station_1_CTD_42_2")
gsub( "Station_1_CTD_42*", "Station_1", test )

# 
# extract data for 2011 to 2019
as.survey <- anoms.season %>% 
  filter( year>=2010 ) %>% 
  spread( key = season, value=temp.anom )
as.survey$summer1 <- c(NA, as.survey$summer[1:length(as.survey$summer)-1])
meta$year <- meta$Year
M <- left_join( meta, as.survey )

# show summer versus winter temperature anomaly
as.survey.all <- anoms.season %>% 
  # filter( year>=2010 ) %>% 
  spread( key = season, value=temp.anom )

ggplot( as.survey, aes(y=summer,x=winter) ) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_smooth(se=T,size=0.5) +
  geom_point(size=2) +  geom_path(size=1) +
  geom_text_repel(label=as.survey$year, box.padding = 0.5) +
  ylab(expression(paste("Summer SST anomaly (",degree,"C)"))) +
  xlab(expression(paste("Winter SST anomaly (",degree,"C)"))) +
  theme_bw()

old.anom <- as.survey.all %>% filter( winter < -1.1 | winter >1 | summer > 1 | summer < -1) %>% filter(year<2015 )
ggplot( as.survey.all, aes(y=summer,x=winter) ) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  stat_ellipse( level = 0.9 ) +
  geom_point(size=1,shape=1) +
  geom_point(data=as.survey,size=2) +  geom_path(data=as.survey,size=1) +
  geom_text_repel(data=old.anom, label=old.anom$year, box.padding = 0.5, size=5) +
  geom_text(data=filter(as.survey.all,year %in% c(2010,2019)), 
            label=c(2010,2019), nudge_x = 0.25, size=5 ) +
  ylab(expression(paste("Summer anomaly (",degree,"C)"))) +
  xlab(expression(paste("Winter anomaly (",degree,"C)"))) +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
ggsave( "R/Figs/winter_summer_SST.svg", width=4, height=4 )





# which species are exceedingly rare?
sort(colSums(comm))
hist(colSums(comm),breaks = seq(0,1600,2))
# which species are most abundant through the time series?
dominance <- apply( comm, 1, function(z) names(z)[order(z,decreasing = T)]  )
dominance[1:3,]

# quick beta diversity by group
# z <- betadiver(comm,"z")
# mod <- betadisper(z, meta$Site)
# plot(mod)
# plot(mod, axes=c(3,1) )
# boxplot(mod)
# anova(mod)
# (mod3B <- betadisper(z, site, type = "median", bias.adjust=TRUE))
# anova(mod3B)
# permutest(mod3B, permutations = 99)
# plot(mod3B)
# boxplot(mod3B)
# 
# adonis2( comm~site*year, by='margin'  )



### Comparing transects at different dates
## betapart temporal
# need two matrices of same dimensions with data from two different years
commpa <- comm 
commpa[commpa>0] <- 1
d12 <- comm[ meta$Year == 2012, ]
d19 <- comm[ meta$Year == 2019, ]
d12p <- commpa[ meta$Year == 2012, ]
d19p <- commpa[ meta$Year == 2019, ]


bpt <- bind_cols(meta[ meta$Year == 2012, ], beta.temp( d12p, d19p ) )
bpt %>% mutate( turn=beta.sim/beta.sor, nest=beta.sne/beta.sor )

# Legendre method (TBI)
tbia <- TBI( d12, d19, method = "%difference" )
plot(tbia)
tbip <- TBI( d12, d19, method = "%difference", pa.tr = TRUE )
plot(tbip)
data.frame( meta[meta$Year == 2019, ], tbia$BCD.mat )
# get TBI matrix for all comparisons to 2012
tbi.sums.abun.loop <- matrix( nrow = length(unique(meta$Year))-1, ncol=8 )
tbi.sums.pa.loop   <- matrix( nrow = length(unique(meta$Year))-1, ncol=8 )
for( i in 2:length(unique(meta$Year)) ){
  di <- comm[ meta$Year == unique(meta$Year)[i], ]
  tbi.sums.abun.loop[i-1,] = unlist(lapply(TBI( d12, di, method = "%difference" )[5:6],function(z) z[1:4]) )
  tbi.sums.pa.loop[i-1,]   = unlist(lapply(TBI( d12, di, method = "%difference", pa.tr = T )[5:6],function(z) z[1:4]) )
}

# add row and column names
tbi.sums.abun <- as.data.frame(tbi.sums.abun.loop)
names(tbi.sums.abun) <- c("losses","gains","total","B%","C-B","stat","p.val",'p.perm')
tbi.sums.abun$comparison <- paste("2012",unique(meta$Year)[-1],sep = "-")
tbi.sums.abun$comp <- 1:7
tbi.sums.abun$second <- unique(meta$Year)[-1]
tbi.sums.abun$sig <- ifelse( tbi.sums.abun$p.val < 0.05, 8, NA )
tbi.sums.pa <- as.data.frame(tbi.sums.pa.loop)
names(tbi.sums.pa) <- c("losses","gains","total","B%","C-B","stat","p.val",'p.perm')
tbi.sums.pa$comparison <- paste("2012",unique(meta$Year)[-1],sep = "-")
tbi.sums.pa$comp <- 1:7
tbi.sums.pa$sig <- ifelse( tbi.sums.pa$p.val < 0.05, 8, NA )

# gather
gath.abun <- tbi.sums.abun %>% select(losses,gains,total,comparison,comp,sig) %>% 
  group_by(comparison,comp) %>% 
  gather(key,value,losses:total)
gath.pa <- tbi.sums.pa %>% select(losses,gains,total,comparison,comp,sig) %>% 
  group_by(comparison,comp) %>% 
  gather(key,value,losses:total)

a <- ggplot( gath.abun, aes(x=comp,y=value,col=key,shape=key)) + 
  geom_line() + geom_point(size=2) + 
  geom_point(data=filter(gath.abun,!is.na(sig)),aes(y=0.7,x=comp),shape=8) +
  ylim(c(0,0.75)) + scale_color_manual(values=c("red","blue","black")) +
  ylab("Dissimilarity") + xlab("Year compared to 2012") +
  scale_x_continuous( breaks=1:7,labels=tbi.sums.abun$second ) +
  theme_bw()
b <- ggplot( gath.pa, aes(x=comp,y=value,col=key,shape=key)) + 
  geom_line() + geom_point(size=2) + 
  # geom_point(data=filter(gath.pa,!is.na(sig)),aes(y=0.7,x=comp),shape=8) +
  ylim(c(0,0.75)) + scale_color_manual(values=c("red","blue","black")) +
  ylab("Dissimilarity") + xlab("Year compared to 2012") +
  scale_x_continuous( breaks=1:7,labels=tbi.sums.abun$second ) +
  theme_bw()

cowplot::plot_grid( a,b, ncol=1, rel_heights = c(1,1) )
ggsave( "R/Figs/beta_temporal_pairs.svg", width =4, height=3 )





### consider dbRDA with factor for year and site (and transect?)
meta$year <- factor(meta$Year, ordered=F)
M$year <- factor(meta$Year, ordered=F)
meta$site <- factor(meta$Site, ordered=F)
db1 <- dbrda( comm~Year+Elevation, distance="bray", data=meta )
db1 <- dbrda( comm~summer1+Elevation, distance="bray", data=M )
db2 <- dbrda( ifelse(comm>0,1,0)~year+Elevation, distance="jaccard", data=meta )
# summary(db1)
anova( db1, by = "terms" )
# db1 <- dbrda( comm~anom.pine.sum.1+Elevation+year, distance="bray", data=M )
# db2 <- dbrda( ifelse(comm>0,1,0)~year+Elevation, distance="jaccard", data=meta )
plot(db1)
os1 <- ordisurf( db1, meta$Elevation )
os1 <- ordisurf( db1, M$summer1 )
ordiellipse( db1, meta$Zone, conf=0.6 )


RsquareAdj(db1)  # explains 25-33% of variation in community composition?
R2 <- eigenvals(db1)/sum(eigenvals(db1))
R2
cummr2 <- R2
for(i in 2:length(eigenvals(db1))){
  cummr2[i] <- cummr2[i]+cummr2[i-1]  
}

##### 
# https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/
anova(db1) # overall test of the significant of the analysis
anova(db1, by="axis", perm.max=500) # test axes for significance
anova(db1, by="terms", permu=200) # test for sign. environ. variables
scores_dbRDA=scores(db1)
site_scores=scores_dbRDA$sites # separating out the site scores, get CAP1 and CAP2 scores
species_scores=scores_dbRDA$species # separating out the species scores
site_scores_environment=cbind(site_scores,select(meta,Elevation,anom.pine.sum.1)) # merge
correlations=cor(site_scores_environment) # calculate correlations
# fix(correlations)
#####

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
write_csv( as.survey, "R/output/sst_anoms_survey.csv" )
surv.cent <- left_join( as.survey, centroids )
centroid2 <- centroids
centroid2$year <- centroid2$year-1
names(centroid2) <- c("shift1","shift2","year")
surv.cent <- left_join( surv.cent, centroid2  )
surv.cent$year2 <- paste0("'",substr(surv.cent$year+1,start=3,stop=4))

# get rough centroids for zones and years
sites.mean <- sites %>% 
  group_by( zone, year ) %>% 
  summarize( dbRDA1 = mean(dbRDA1), dbRDA2 = mean(dbRDA2) )

as.zone <- left_join( sites.mean, as.survey)
as.zone$year2 <- NA
as.zone$year2[ as.zone$zone == "MID" & as.zone$year %in% c(2012,2019)] <- c("2012","2019")

library(broom)
tidied <- sites.mean %>%
  nest(data = -zone) %>% 
  mutate(
    fit = map(data, ~ lm(dbRDA1 ~ dbRDA2, data = .x)),  # S3 list-col
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  ) %>% 
  unnest(tidied)
unnest.predict <- tidied %>% unnest(augmented)
zone.reg <- unnest.predict %>% 
  filter(term == "dbRDA2") %>% 
  select(zone, dbRDA2, dbRDA1 = .fitted )

# model difference in slopes for elevation trajectory with year
lmm1 <- lme4::lmer( dbRDA2~dbRDA1 + (dbRDA1|year), data=sites )
lmm1 <- lm( dbRDA2~dbRDA1, data=sites )
summary(lmm1)
coef(lmm1)
# windows(3.5,3.5)
par(mar=c(4,6,0,2)+0.1,pty='s',las=1)
plot( y = coef(lmm1)$year$dbRDA1, x= 2012:2019, axes=F,
      ylab="slope of\ndbRDA2~dbRDA1\n(BLUPs)", xlab="year" )
abline( h=0, lty=2 )
axis(1, at = 2012:2019, labels=F )
text(x=2012:2019, par("usr")[3] - 0.12, labels = 2012:2019, srt=-45, xpd=T, pos=3 ) #paste0("'",12:19)
axis(2)
box()
plot( dbRDA2~dbRDA1, data=sites)
predict(lmm1)
preds <- merTools::predictInterval(lmm1, newdata=sites, n.sims=999)
preds <- predict(lmm1, newdata=sites, n.sims=999)
preds.df <- bind_cols( sites, data.frame(fit=preds) )
preds.env <- left_join( preds.df, select(surv.cent, year, summer1) )
#

a <- ggplot( sites, aes(x=dbRDA1,y=dbRDA2)) + 
  # stat_contour(data = ordi.na, aes(x = x, y = y, z = z, colour = rev(..level..)),
  #              binwidth = 20)+ #can change the binwidth depending on how many contours you want
  # geom_smooth( aes(group=year), method='lm', se=F, size=0.5,col='gray85') +
  geom_point( data=sites, aes(shape=zone), size=1 )  +
  # geom_smooth( aes(group=1), method='lm', se=F, size=0.5,col='gray25') +
  geom_smooth( data=preds.env, aes(y=fit, group=year,col=summer1), method='lm', se=F, size=0.5 ) +
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
a
ggsave("R Code and Analysis/Figs/rda_bwr_pal2_temp.svg", width=4, height=3 )


# flip axes
sites$zone <-  factor( sites$zone, levels = c("HIGH","MID","LOW"))
library(ggnewscale)
ggplot( sites, aes(x=dbRDA1,y=dbRDA2)) + 
  # stat_contour(data = ordi.na, aes(x = x, y = y, z = z, colour = rev(..level..)),
               # binwidth = 20)+ #can change the binwidth depending on how many contours you want
  # geom_smooth( aes(group=year), method='lm', se=F, size=0.5,col='gray85') +
  # geom_point( aes(shape=zone), size=1, col='darkslategrey' )  +
  geom_point( aes(fill = zone), size=1, shape = 21, alpha=0.75,col = 'slategrey' )  +
  scale_fill_manual(values = c("white","grey","black")) +
  labs(fill = "Zone", col = "Summer\nSST anomaly", lty = "Zone" ) +
  # geom_point( data = filter(sites, zone == "HIGH"), size=1, shape = 21, fill="white", alpha=0.5 )  +
  # geom_point( data = filter(sites, zone == "MID"), size=1, shape = 21, fill="gray", alpha=0.5 )  +
  # geom_point( data = filter(sites, zone == "LOW"), size=1, shape = 21, fill="black", alpha=0.5 )  +
  # geom_path( data = as.zone, aes(lty=zone)) +
  stat_ellipse( aes(lty=zone), level = 0.6, col='slategrey', lwd = 0.33 ) +
  # geom_smooth( aes(group=1), method='lm', se=F, size=0.5,col='gray25') +
  # geom_smooth( data=preds.env, aes(y=fit, group=year,col=summer1), method='lm', se=F, size=0.5 ) +
  # geom_path( data=surv.cent, aes(col=summer1),size=1 ) +
  # geom_point( data=surv.cent, aes(fill=summer1), size=3, shape=21 ) +
  geom_path( data=as.zone, aes(col=summer1,group=zone),size = 1, ) +
  geom_path( data=as.zone, aes(group=zone),size = 0.25, col='black' ) +
  ggnewscale::new_scale_fill() +
  geom_point( data=as.zone, aes(fill=summer1, group=zone), size=3, shape=21 ) +
  # geom_path( data = zone.reg, aes(group=zone), col='black') +
  annotate("text", label = "2012", y = -0.29, x = 0.01 ) +
  annotate("text", label = "2019", y = 0.29, x = 0.01 ) +
  xlab(paste0(names(sites)[1],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(sites)[2],' (',round(R2[2],3)*100, '%)')) +
  scale_fill_gradientn(colours=pal2(100),limits=c(-2,2)) +
  scale_color_gradientn(colours=pal2(100),limits=c(-2,2)) +
  scale_shape_manual(values = c(2,0,1)) +
  # scale_linetype_manual(values = c(3,2,1)) +
  theme_bw() + theme( panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank()) +
  theme( legend.title = element_text(size=8),
         legend.text = element_text(size=8),
         legend.key.size = unit(0.59, "cm")) +
  labs(fill = "Summer\nSST anomaly", col = "Summer\nSST anomaly", shape = "Zone", lty = "Zone" ) +
  coord_flip()
ggsave("R/Figs/rda_family_bwr_pal2_temp_flip.svg", width=4, height=3 )
#










windows(2,2)
ggplot( surv.cent, aes(y=dbRDA2,x=winter) ) + geom_point()
ggplot( surv.cent, aes(y=shift2,x=winter) ) + geom_point()
ggplot( surv.cent, aes(y=shift2,x=summer) ) + 
  geom_smooth(method="lm",se=T) + 
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



