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
  group_by(taxon_lumped3) %>% 
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
d$taxon_lumped[d$taxon_lumped=="Mytilus trossulus"] <- "Mytilus sp."
d$taxon_lumped2[d$taxon_lumped2=="Haliclona"] <- "Tunicata/Porifera"
d$taxon_lumped2[d$taxon_lumped2=="Halichondria"] <- "Tunicata/Porifera"

d <- d %>% 
  filter( !(taxon_lumped %in% c("articulated coralline","Unknown red blade",
                                "Colonial Diatoms", "Pleonosporium vancouverianum",
                                "Acrochaetium sp.", "Fauchea") ) )



## taxa not in final list
# remove "Chiharaea bodegensis", "Desmarestia herbacea", "Ectocarpus siliculosus"
# remove "Gloiocladia sp.", "Kornmannia leptoderma", "Petalonia fascia", "Tiffaniella snyderae", "Rhodymenia sp."
# "Phycodrys sp." = Polyneura latissima
d <- d %>% 
  filter( !(taxon_lumped3 %in% c("Chiharaea bodegensis", "Desmarestia herbacea", "Ectocarpus siliculosus",
                                 "Ulothrix-Urospora sp.") ) )
d <- d %>% 
  filter( !(taxon_lumped3 %in% c("Gloiocladia sp.", "Kornmannia leptoderma", "Petalonia fascia", "Tiffaniella snyderae", "Rhodymenia sp.") ) )
d$taxon_lumped3[ d$taxon_lumped3 == "Phycodrys sp."] <- "Polyneura latissima"
sort(unique(d$taxon_lumped3))

# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which will reduce the size of the dataset a bit
# restrict this to seaweeds and sessile invertebrates
d.simple <- d %>%
  mutate( taxon = gsub(" ",".",taxon_lumped3) ) %>% 
  group_by( UID, Year, Site, Zone, Quadrat, taxon, funct_2021 ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T)) 

# write to disk so we are using the same dataset here and in HMSC
write_csv(d.simple, "R/output/data_select_rda_HMSC.csv")

# pivot wider
dwide <-  d.simple %>%
  group_by( Year, Site, Zone, Quadrat, taxon ) %>%
  summarise( Abundance=mean(Abundance) ) %>% 
  pivot_wider( names_from = taxon, values_from = Abundance, values_fill=0 ) %>% 
  ungroup()
# Remove taxa with less than one percent 
# remove taxa with less than two occurrences (singletons)? some are good
occurrences <- colSums(ifelse( dwide[,-c(1:4)] == 0, 0, 1 ) )
total_abundances <- colSums(dwide[,-c(1:4)])
sort(occurrences, decreasing = T)
sort(total_abundances, decreasing = T)
# taxa_keep <- names(total_abundances[total_abundances>=1])

dmean <- d %>% 
  group_by( Year, Site, Zone, taxon = taxon_lumped3 ) %>%
  # filter( taxon %in% taxa_keep ) %>% 
  summarise( Abundance=mean(Abundance) )


# spread out
d.comm.mean <- dmean %>%
  pivot_wider( names_from = taxon, values_from = Abundance, values_fill=0 ) %>%
  ungroup() %>% 
  mutate( Zone = factor(Zone, levels=c("LOW","MID","HIGH")) ) %>% 
  arrange( Year, Site, Zone )

# d.comm.mean <- d.comm.mean %>% filter( !is.na(Elevation) )

# remove UID column from community data
meta <- d.comm.mean[ ,1:3 ]
meta <- left_join(meta,dmeanelev)
comm <- as.matrix(d.comm.mean[,-c(1:3)])
dim(comm)
# interrogate the dataset
sort( colSums(comm), decreasing = T )


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
anoms <-  read_csv("Data/R code for Data Prep/output from R/Lightstation_monthly_anomaly.csv")
plot(anoms$sal[anoms$site == "mccinnis" & anoms$year >= 2000], type='l')
filter(anoms, site == "mccinnis") %>%  summarize(mean.sal = mean(sal,na.rm=T))
# collect temperature, sal, precip and use these in a PCA
temp_wide <- anoms %>% 
  select( site, year, month, temp.anom ) %>% 
  pivot_wider( names_from = site, values_from = temp.anom, names_prefix = "temp_" ) 
sal_wide <- anoms %>% 
  filter( site != "addenbroke" ) %>% 
  select( site, year, month, sal.anom ) %>% 
  pivot_wider( names_from = site, values_from = sal.anom, names_prefix = "sal_" ) 
precip <- anoms %>% 
  filter( site == "addenbroke" ) %>% 
  select( year, month, precip.anom ) %>% 
  mutate( precip.anom = -precip.anom )

alld <- full_join(full_join(temp_wide, sal_wide),precip)
# alld <- full_join(temp_wide,precip)
alld <- alld %>% 
  unite(date, year,month,remove = F) %>% 
  mutate( date = lubridate::ym(date) ) %>% 
  select( -precip.anom )
# PCA
dna <- alld[ !(apply( alld, 1, function(z) any(is.na(z)) )), ]
pca1 <- princomp( select(dna,temp_pine:sal_mccinnis ) )
summary(pca1)
plot(pca1)
biplot(pca1,scale = 0)
biplot(pca1,scale = 0, choice = c(2,3))
# missing and imputation ---- see http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html -----
library(missMDA)
# ignore the earlier years of the dataset
filtd <- filter( alld,year > 1977 & year < 2020 )
nb <- estim_ncpPCA( select(filtd,temp_pine:sal_mccinnis ), method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
#(available methods include GCV to approximate CV)
nb$ncp
res.comp <- imputePCA( select(filtd,temp_pine:sal_mccinnis ), ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set
imp <- res.comp$completeObs
library(FactoMineR)
png(file="R/Figs/BC_Lightstation_PCA.png", res = 600, width = 4.5, height = 4.5, units = "in")
# windows(4.5,6)
par(mfrow = c(2,1),las = 3, mar = c(3,4,3,1)+0.1 )
res.pca <- PCA(imp, ncp = nb$ncp, graph = TRUE)
dev.off()
res.pca$var
summary(res.pca)
# windows(2,2)
par( mar = c(2.25,2.25,1,1)+0.1, las = 0, cex=0.8, pty = "s" )
plot( res.pca, choix = "var", cex = 0.8, title = "", 
     col.var = c("darkslateblue","darkslateblue", "darkslateblue","firebrick4","firebrick4"),
     label = "none", graph.type = "classic" )
plot(res.pca, choix = "var", axes = 2:3, cex = 0.8)
plot(res.pca, choix = "var", axes = 3:4, cex = 0.8)
dimdesc( res.pca )
pcscores <- data.frame( res.pca$ind$coord )
names(pcscores) <- paste0("pca",1:ncol(pcscores))

# biplot(pca1,scale = 0, choice = c(3,4))
# ccf(pcscores$pca1, pcscores$pca2 )
ccf(pcscores$pca1, pcscores$pca2 )

# dates for events
# 1997-1998 El Nino
warm_times <- lubridate::ymd(c("1997-06-01","1998-06-01","2013-12-01","2016-01-01"))
lubridate::ymd(c("2013-12-01","2016-01-01"))
png(file="R/Figs/BC_Lightstation_PCA_timeseries.png", res = 600, width = 5, height = 5, units = "in")
par( mar=c(2,4,0,1)+0.1, mfrow=c(3,1), las = 1 )
plot(x = filtd$date, y = pcscores$pca1, type = 'l', col = 'slateblue', ylab = "PCA1" ); abline(h = 0)
lines( lowess(x = filtd$date, y = pcscores$pca1, f = 1/20, iter = 10), col = "darkslateblue", lwd=2)
abline( v = warm_times, lty = 4, col = "slategrey" )
axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
plot(x = filtd$date, y = pcscores$pca2, type = 'l', col = 'firebrick', ylab = "PCA2" ); abline(h = 0)
lines( lowess(x = filtd$date, y = pcscores$pca2, f = 1/20, iter = 10), col = "firebrick4", lwd=2)
abline( v = warm_times, lty = 4, col = "slategrey" )
axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
plot(x = filtd$date, y = pcscores$pca3, type = 'l', col = 'darkorange', ylab = "PCA3" ); abline(h = 0)
lines( lowess(x = filtd$date, y = pcscores$pca3, f = 1/20, iter = 10), col = "darkorange4", lwd=2)
abline( v = warm_times, lty = 4, col = "slategrey" )
axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
dev.off()
dna <- bind_cols( filtd, pcscores )

# calculate winter and summer temperature anomalies
library(zoo)
yq <- as.yearqtr( as.yearmon( paste(dna$month,dna$year,sep="/"), "%m/%Y") + 1/12)
dna$season <- factor(format(yq, "%q"), levels = 1:4,
                       labels = c("winter", "spring", "summer", "fall"))
# add year groupings - for instance, seaweeds in summer 2016 would be influenced by conditions over the previous year, 
#                      so, count previous summer, fall, and current winter and spring towards a give year
dna$survey.year <- dna$year 
dna$survey.year[ dna$month %in% 6:12 ] <- dna$survey.year[ dna$month %in% 6:12 ] + 1
dna %>% select(date, season, survey.year, temp_pine ) %>% filter( survey.year > 2010)
anoms.season <- dna %>%
  group_by( survey.year, season ) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
ggplot(anoms.season, aes(x=survey.year,y=pca2,col=season)) + geom_line()
ggplot(anoms.season, aes(x=survey.year,y=pca1,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
ggplot(filter(anoms.season, survey.year>=2010), aes(x=survey.year,y=pca1,col=season)) + geom_path() + geom_point()

# include previous year as summer to summer rather than winter to winter
# anoms.season$yeargroup <-  anoms.season$year + 1
# anoms.season$yeargroup <-  ifelse(anoms.season$season %in% c("winter","spring"), anoms.season$yeargroup - 1, anoms.season$yeargroup)
# compare with means across years
anoms.annual <- anoms.season %>%  # or dna if just using annual means
  group_by( survey.year ) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)
ggplot(filter(anoms.annual, year>=2010), aes(x=survey.year,y=pca1)) + geom_path() + geom_point()
# ggplot(filter(anoms.annual, year>=2010), aes(x=survey.year,y=pca2)) + geom_path() + geom_point()
# ggplot(filter(anoms.annual, year>=2010), aes(x=survey.year,y=pca3)) + geom_path() + geom_point()
# ggplot(filter(anoms.annual), aes(x=survey.year,y=sal_mccinnis)) + geom_hline(yintercept = 0) +  geom_path(col = 'slateblue') + geom_point(col = 'slateblue') + ylim(c(-1.75,1.75))
# ggplot(filter(anoms.annual), aes(x=survey.year,y=sal_pine)) + geom_hline(yintercept = 0) +  geom_path(col = 'slateblue') + geom_point(col = 'slateblue') + ylim(c(-1.75,1.75))
# ggplot(filter(anoms.annual, year>=2010), aes(x=year,y=temp_pine)) + geom_path() + geom_point()

# # figure out how to make an anomaly plot with vertical lines from zero
# # pairwise correlations among seasonal anomalies
# as.all <- anoms.season %>% 
#   spread( key = season, value=temp.anom )
# winter = ts(as.all$winter)
# spring = ts(as.all$spring)
# summer = ts(as.all$summer)
# fall   = ts(as.all$fall)
# ccf(winter, summer)
# ccf(spring, summer)
# ccf(summer, fall)
# ccf(winter, fall)
# cor(winter, summer) # 0.77
# cor(spring, summer) # 0.88
# cor(winter, fall)   # 0.69
# cor(fall, summer)   # 0.73

# extract data for 2011 to 2019
as.survey <- anoms.annual %>% 
  filter( survey.year>=2010 ) 
# as.survey <- anoms.season %>% 
  # filter( year>=2010 ) # %>% spread(season, Comp.2)

# read in PCA data using raw and then imputed values
as2 <- read_csv( "Data/R code for Data Prep/Output from R/Lightstation_raw_PCA_impute.csv" )

# pick
as.use <- as2 # as2 or as.survey

# add lag into pca axes (this is already done above)
# # as.survey$summer1 <- c(NA, as.survey$summer[1:length(as.survey$summer)-1])
# as.survey$pca11 <- c(NA, as.survey$pca1[1:length(as.survey$pca1)-1])
# as.survey$pca21 <- c(NA, as.survey$pca2[1:length(as.survey$pca2)-1])
# as.survey$pca31 <- c(NA, as.survey$pca3[1:length(as.survey$pca3)-1])
meta$year <- meta$Year
M <- left_join( meta, as.use, by = c("year" = "survey.year") )

# # show summer versus winter temperature anomaly
# as.survey.all <- anoms.season %>% 
#   # filter( year>=2010 ) %>% 
#   spread( key = season, value=temp.anom )
# 
# ggplot( as.survey, aes(y=summer,x=winter) ) + 
#   geom_hline(yintercept=0) +
#   geom_vline(xintercept=0) +
#   geom_smooth(se=T,size=0.5) +
#   geom_point(size=2) +  geom_path(size=1) +
#   geom_text_repel(label=as.survey$year, box.padding = 0.5) +
#   ylab(expression(paste("Summer SST anomaly (",degree,"C)"))) +
#   xlab(expression(paste("Winter SST anomaly (",degree,"C)"))) +
#   theme_bw()
# 
# old.anom <- as.survey.all %>% filter( winter < -1.1 | winter >1 | summer > 1 | summer < -1) %>% filter(year<2015 )
# ggplot( as.survey.all, aes(y=summer,x=winter) ) + 
#   geom_hline(yintercept=0) +
#   geom_vline(xintercept=0) +
#   stat_ellipse( level = 0.9 ) +
#   geom_point(size=1,shape=1) +
#   geom_point(data=as.survey,size=2) +  geom_path(data=as.survey,size=1) +
#   geom_text_repel(data=old.anom, label=old.anom$year, box.padding = 0.5, size=5) +
#   geom_text(data=filter(as.survey.all,year %in% c(2010,2019)), 
#             label=c(2010,2019), nudge_x = 0.25, size=5 ) +
#   ylab(expression(paste("Summer anomaly (",degree,"C)"))) +
#   xlab(expression(paste("Winter anomaly (",degree,"C)"))) +
#   theme_bw() +
#   theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
# ggsave( "R/Figs/winter_summer_SST.svg", width=4, height=4 )





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

cowplot::plot_grid( a,b, ncol=1, rel_heights = c(1,1), labels = "auto", vjust = 0.8 )
ggsave( "R/Figs/beta_temporal_pairs.svg", width =4, height=3 )







###########   distance-based redundancy discrimanant analysis



### consider dbRDA with factor for year and site (and transect?)
with( distinct(select(M,pca1,pca2)), ccf(pca1,pca2) )
with( distinct(select(M,pca1,Year)), ccf(pca1,Year) )
with( distinct(select(M,pca1,pdo)), ccf(pca1,pdo) )
with( distinct(select(M,Elevation,Year)), ccf(Elevation,Year) )
with( distinct(select(M,Elevation,pca2)), ccf(Elevation,pca2) )
with(M, ccf(pca1,Year) )

meta$year <- factor(meta$Year, ordered=F)
M$year <- factor(meta$Year, ordered=F)
meta$site <- factor(meta$Site, ordered=F)
# # compare transformations, note sure what to make of this
# source("box.cox.chord.R")
# Y <- comm
# dbraw      <- dbrda( Y~Elevation+pca1, distance="bray", data=M )
# Y <- ifelse(comm>0,1,0)
# dbpa <- dbrda( Y~Elevation+pca1, distance="bray", data=M )
# Y <- sqrt(comm)
# hist(Y[Y>0])
# dbroot     <- dbrda( Y~Elevation+pca1, distance="bray", data=M )
# Y <- sqrt(sqrt(comm))
# dbrootroot <- dbrda( Y~Elevation+pca1, distance="bray", data=M )
# # find a box-cox transform
# picks <- seq(0,1,by=0.01)
# shap <- NA
# for(i in 1:length(picks)){
#   tempY <- box.cox.chord( comm, bc.exp = picks[i] )
#   tempdb <- dbrda( tempY~Elevation+pca1+pca2, distance="bray", data=M )
#   shap[i] <- shapiro.test(resid(tempdb))$p.value
# }
# plot(shap)
# picks[which(shap==max(shap))]
# # pretty clear winner is 0.77, but this seems to depend a lot on model formulation
# Y <- box.cox.chord( comm, bc.exp = 0.77 )
# hist(Y[Y>0])
# dbbox <- dbrda( Y~Elevation+pca1+pca2, distance="bray", data=M )

###
###
#
Y <- comm
db0 <- dbrda( Y~Elevation, distance="bray", data=M )
dbpca <- dbrda( Y~Elevation+pca1+pca2+pca3+pca4, distance="bray", data=M )
anova( db0, dbpca )
db1 <- dbrda( Y~Elevation+pca1+pca2, distance="bray", data=M )
anova( db1, dbpca)

# presence-absence
db2 <- dbrda( ifelse(comm>0,1,0)~Elevation+pca1+pca2, distance="jaccard", data=M )

#
#
db <- db1
plot(db)
os1 <- ordisurf( db, meta$Elevation )
os1 <- ordisurf( db, M$pca1 )
os1 <- ordisurf( db, M$pca2 )
ordiellipse( db, meta$Zone, conf=0.6 )


RsquareAdj(db)  # explains 25-33% of variation in community composition?
R2 <- eigenvals(db)/sum(eigenvals(db))
R2
cummr2 <- R2
for(i in 2:length(eigenvals(db))){
  cummr2[i] <- cummr2[i]+cummr2[i-1]  
}
cummr2[1:7]

##### 
# https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/
anova(db) # overall test of the significant of the analysis
anova(db, by="axis", perm.max=500) # test axes for significance
anova(db, by="terms", permu=999) # test for sign. environ. variables
scores_dbRDA=scores(db)
site_scores=scores_dbRDA$sites # separating out the site scores, get CAP1 and CAP2 scores
species_scores=scores_dbRDA$species # separating out the species scores
site_scores_environment=cbind(site_scores,select(M,Elevation,pca21,pca31)) # merge
correlations=cor(site_scores_environment) # calculate correlations
correlations
#####

# extract axes
nax <- 1:2
sites     <- data.frame(scores(db,choices = nax, scaling = 0)$sites)
sites$zone <- factor( meta$Zone, levels=c("LOW","MID","HIGH") )
sites$year <- meta$Year
sites$transect <-  factor(unite( select(meta,Site,Zone), "transect", Site, Zone )$transect)
centroids <- data.frame(scores(db,choices = nax, scaling=0)$centroids)
centroids$year <- 2012:2019

# make the values of dbRDA2 go from smallest to largest over time
is.increasing <- sites %>% 
  filter(year %in% c(2012,2019)) %>% 
  select(year,dbRDA2) %>% 
  group_by(year) %>% 
  summarize( dbRDA2 = mean(dbRDA2) ) %>% 
  summarize( diff = diff(dbRDA2) ) 
if(is.increasing$diff < 0) {
  sites$dbRDA2 <- -sites$dbRDA2
}

ordi.grid <- os1$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.na <- data.frame(na.omit(ordi)) #gets rid of the nas
ordi.na #looks ready for plotting!


### compare centroids to temperature anomaly data
centroids
# as.survey$summer1 <- c(NA, as.survey$summer[1:length(as.survey$summer)-1])
write_csv( as.use, "R/output/sst_anoms_survey.csv" )
surv.cent <- left_join( as.use, centroids )
centroid2 <- centroids
centroid2$year <- centroid2$year-1
names(centroid2) <- c("shift1","shift2","year")
surv.cent <- left_join( surv.cent, centroid2  )
surv.cent$year2 <- paste0("'",substr(surv.cent$year+1,start=3,stop=4))

# get rough centroids for zones and years
sites.mean <- sites %>% 
  group_by( zone, year ) %>% 
  summarize( dbRDA1 = mean(dbRDA1), dbRDA2 = mean(dbRDA2) )

as.zone <- left_join( sites.mean, as.survey, by = c("year" = "survey.year"))
as.zone$year2 <- NA
as.zone$year2[ as.zone$zone == "MID" & as.zone$year %in% c(2012,2019)] <- c("2012","2019")
as.zone$pca1_scale <- scale(as.zone$pca1)
as.zone$pca2_scale <- scale(as.zone$pca2)
as.zone$pca3_scale <- scale(as.zone$pca3)
as.zone$pca11_scale <- scale(as.zone$pca11)
as.zone$pca21_scale <- scale(as.zone$pca21)
as.zone$pca31_scale <- scale(as.zone$pca31)

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

#

a <- ggplot( sites, aes(x=dbRDA1,y=dbRDA2)) + 
  # stat_contour(data = ordi.na, aes(x = x, y = y, z = z, colour = rev(..level..)),
  #              binwidth = 20)+ #can change the binwidth depending on how many contours you want
  # geom_smooth( aes(group=year), method='lm', se=F, size=0.5,col='gray85') +
  geom_point( data=sites, aes(shape=zone), size=1 )  +
  # geom_smooth( aes(group=1), method='lm', se=F, size=0.5,col='gray25') +
  # geom_smooth( data=preds.env, aes(y=fit, group=year,col=summer1), method='lm', se=F, size=0.5 ) +
  # geom_path( data=surv.cent, aes(col=summer1),size=1 ) +
  # geom_point( data=surv.cent, aes(fill=summer1), size=3, shape=21 ) + 
  # geom_text_repel( data=surv.cent, label=surv.cent$year, 
  #                  fontface="bold",
  #                  point.padding = 0.3, box.padding = 0.1 ) +
  scale_shape_manual(values = c(1,0,2)) +
  xlab(paste0(names(sites)[1],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(sites)[2],' (',round(R2[2],3)*100, '%)')) +
  scale_fill_gradientn(colours=pal2(100),limits=c(-2,2)) +
  scale_color_gradientn(colours=pal2(100),limits=c(-2,2)) +
  theme_bw() + theme( panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
a
# ggsave("R Code and Analysis/Figs/rda_bwr_pal2_temp.svg", width=4, height=3 )


# flip axes
sites$zone <-  factor( sites$zone, levels = c("HIGH","MID","LOW"))
library(ggnewscale)
summary(sites)

# adjust pca scores for paths so that we show environment "leading up to" each survey
as.zone$pca11 <- c(as.zone$pca1[2:length(as.zone$pca1)],NA)
xrange <- range(sites$dbRDA1)*1.01
yrange <- range(sites$dbRDA2)*1.01
zone_y <- 0.35
as.zone %>% filter(year %in% c(2012,2019), zone == "MID") %>% 
  select(year, dbRDA1, dbRDA2)
rda1 <- ggplot( sites, aes(x=dbRDA1,y=dbRDA2)) + 
  # geom_point( aes(fill = zone), size=1, shape = 21, alpha=1,col = 'slategrey' )  +
  # scale_fill_manual(values = c("white","grey","black")) +
  labs(fill = "Zone", col = "Temperature\nanomaly", lty = "Zone" ) +
  # stat_ellipse( aes(lty=zone), level = 0.6, col='slategrey', lwd = 0.33 ) +
  geom_path( data=as.zone, aes(col=pca11,group=zone),size = 1 ) +
  ggnewscale::new_scale_fill() +
  geom_point( data=as.zone, aes(fill=pca1, group=zone), size=3, shape=21 ) +
  annotate("text", label = "2012", y = -0.16, x = 0.023, hjust = 1 ) +
  annotate("text", label = "2019", y = 0.19, x = 0.055, hjust = 0 ) +
  annotate("text", label = "HIGH", y = zone_y, x = 0.2, hjust = 1 ) +
  annotate("text", label = "MID", y = zone_y, x = 0, hjust = 1  ) +
  annotate("text", label = "LOW", y = zone_y, x = -0.2, hjust = 1 ) +
  xlab(paste0(names(sites)[1],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(sites)[2],' (',round(R2[2],3)*100, '%)')) +
  scale_fill_gradientn(colours=pal2(100),limits=c(-2.84,2.84)) +
  scale_color_gradientn(colours=pal2(100),limits=c(-2.84,2.84)) +
  scale_shape_manual(values = c(2,0,1)) +
  theme_bw() + theme( panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank()) +
  theme( legend.title = element_text(size=8),
         legend.text = element_text(size=8),
         legend.key.size = unit(0.59, "cm")) +
  guides(lty = FALSE) +
  labs(fill = "Environment\nPC1", col = "Environment\nPC1", shape = "Zone", lty = "Zone" ) +
  coord_flip(xlim = xrange, ylim = yrange, clip = "off")
rda1
ggsave("R/Figs/rda_bwr_pal2_temp_flip.svg", width=4, height=3 )

# show trajectories
sites$end1 <- sites$dbRDA1
sites$end1[!(sites$year %in% c(2012,2019))] <- NA
sites$end2 <- sites$dbRDA2
sites$end2[!(sites$year %in% c(2012,2019))] <- NA
# point size for beginning and end
sites$endsize <- sites$end1
sites$endsize[sites$year == 2012] <- 1
sites$endsize[sites$year == 2019] <- 2

sites <- sites %>% 
  separate(transect, c("Site", "Zone"), sep = "_", remove = FALSE) %>% 
  mutate( Site = factor(Site, levels = c("Foggy Cove","Fifth Beach","North Beach")))
rda2 <- ggplot( sites, aes(x=dbRDA1,y=dbRDA2)) + 
  geom_path( aes(group = transect, lty=Site), lwd=0.33, alpha=0.5,col = 'slategrey') +
  geom_point( aes(fill = zone, x = end1, y = end2,size = endsize), shape = 21, alpha=1,col = 'black' )  +
  scale_fill_manual(values = c("white","grey","black")) +
  scale_size(range = c(1.5,3)) +
  annotate("text", label = "HIGH", y = zone_y, x = 0.2, hjust = 1 ) +
  annotate("text", label = "MID", y = zone_y, x = 0, hjust = 1  ) +
  annotate("text", label = "LOW", y = zone_y, x = -0.2, hjust = 1 ) +
  labs(fill = "Zone", col = "Temperature\nanomaly", lty = "Site" ) +
  xlab(paste0(names(sites)[1],' (',round(R2[1],3)*100, '%)')) +
  ylab(paste0(names(sites)[2],' (',round(R2[2],3)*100, '%)')) +
  theme_bw() + theme( panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank()) +
  theme( legend.title = element_text(size=8),
         legend.text = element_text(size=8),
         legend.key.size = unit(0.59, "cm")) +
  guides(size = FALSE) +
  coord_flip(xlim = xrange, ylim = yrange)
rda2
ggsave("R/Figs/rda_bwr_pal2_temp_flip_transect.svg", width=4, height=3 )
cowplot::plot_grid(rda1,rda2,labels = "auto", axis="rlbt", align="hv")
ggsave("R/Figs/rda_bwr_pal2_temp_flip_2panel.svg", width=8, height=2.8 )
#














#










#







# take a closer look at trajectories relative to dbRDA model
library(vegclust)
library(ecotraj)
## Trajectory statistics
D_use <- comm
# define site as the particular trasect
Zone <- factor( meta$Zone, levels=c("LOW","MID","HIGH"), labels=c("Low","Mid","High") )
Site <- factor( meta$Site, labels=c("5","F","N") )
transect <- paste( Site, Zone, sep="." )
transect <- factor( site, levels= c("5.Low","5.Mid","5.High","F.Low","F.Mid","F.High","N.Low","N.Mid","N.High"),
                labels= c("5.L","5.M","5.H","F.L","F.M","F.H","N.L","N.M","N.H"))
# define year
year <- meta$Year
trajectoryLengths(        D_use, transect, year ) 
trajectoryAngles(         D_use, transect, year )
trajectoryAngles(         D_use, transect, year, all=TRUE ) # high degree of angle homogeneity
trajectoryDirectionality( D_use, transect, year ) # despite longer segements in LOW, MID and HIGH often are more directional
plot( trajectoryDirectionality( D_use, transect, year ), cex=3, pch=21, bg=cols ) # similar directionality among transects
trajectoryProjection(     D_use, 1,2:6 )
trajectoryConvergence(    D_use, transect, year, symmetric = FALSE )
trajectoryPlot( x, transect, year, axes=1:2  )

