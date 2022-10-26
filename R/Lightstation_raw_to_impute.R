## investigate raw data from lightstations instead of monthly averages or anomalies
library(tidyverse)
d <- read_csv("Data/R code for Data Prep/Output from R/Lightstation_raw.csv") # see Lightstation_spectral_anomaly.R


# arrange with columns for each measurement at each site
names(d)
head(d)

# take one piece of data at a time and then merge them
dtemp <- d %>% select(date,temp,site) %>% group_by(date) %>% 
  pivot_wider( names_from = site, values_from = temp, names_prefix = "temp_")
dsal <- d %>% select(date,sal,site) %>% group_by(date) %>% 
  pivot_wider( names_from = site, values_from = sal, names_prefix = "sal_")

# merge and clean
dwide <- left_join(dtemp,dsal)# %>% select(-sal_addenbroke)

drecent <- dwide %>% filter( date >= "1978-01-01" & date <= "2019-10-31") %>% ungroup()


# check the data for outliers, etc.
summary(drecent)


# use missing data algorithms
# missing and imputation ---- see http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html -----
library(missMDA)
# # ignore the earlier years of the dataset
# nb <- estim_ncpPCA( select(drecent,-date), method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
# remove salinity
nb <- estim_ncpPCA( select(drecent,temp_pine, temp_mcinnes, temp_addenbroke), method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
#(available methods include GCV to approximate CV)
nb$ncp
# res.comp <- imputePCA( select(drecent,-date), ncp = nb$ncp) # iterativePCA algorithm
res.comp <- imputePCA( select(drecent,temp_pine, temp_mcinnes, temp_addenbroke), ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set
imp <- res.comp$completeObs
library(FactoMineR)
res.pca <- PCA(imp, ncp = nb$ncp, graph = TRUE)
summary(res.pca)
res.pca$eig
res.pca$var
plot(res.pca, choix = "var", cex = 0.8)
plot(res.pca, choix = "var", axes = 2:3, cex = 0.8)
png(file="R/Figs/BC_Lightstation_daily_loadings.png", res = 600, width = 3, height = 3, units = "in")
par(mar = c(5,4,2,2)+0.1, pty="s" )
plot( res.pca, choix = "var", cex = 0.8, title = "", 
      col.var = c("darkslateblue","darkslateblue", "darkslateblue","firebrick4","firebrick4"),
       graph.type = "classic", label="none" ) #label = "none",
dev.off()
# dimdesc( res.pca )
pcscores <- data.frame( res.pca$ind$coord )
names(pcscores) <- paste0("pca",1:ncol(pcscores))

# merge with dates and write to disk
pca.meta <- bind_cols(drecent,pcscores)
write_csv( pca.meta, "Data/R code for Data Prep/Output from R/Lightstation_raw_PCA_impute_daily.csv" )


ccf(pcscores$pca1, pcscores$pca2 )

# dates for events
# 1997-1998 El Nino
warm_times <- lubridate::ymd(c("1997-06-01","1998-06-01","201408-01","2016-12-31"))
lubridate::ymd(c("2013-12-01","2016-01-01"))
png(file="R/Figs/BC_Lightstation_PCA_timeseries.png", res = 600, width = 7, height = 5, units = "in")
par( mar=c(2,4,0,1)+0.1, mfrow=c(ncol(pcscores),1), las = 1 )
plot(x = drecent$date, y = pcscores$pca1, type = 'l', col = 'slateblue', ylab = "PCA1" ); abline(h = 0)
lines( lowess(x = drecent$date, y = pcscores$pca1, f = 1/20, iter = 10), col = "darkslateblue", lwd=2)
abline( v = warm_times, lty = 4, col = "slategrey" )
axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
plot(x = drecent$date, y = pcscores$pca2, type = 'l', col = 'firebrick', ylab = "PCA2" ); abline(h = 0)
lines( lowess(x = drecent$date, y = pcscores$pca2, f = 1/20, iter = 10), col = "firebrick4", lwd=2)
abline( v = warm_times, lty = 4, col = "slategrey" )
axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
# plot(x = drecent$date, y = pcscores$pca3, type = 'l', col = 'darkorange', ylab = "PCA3" ); abline(h = 0)
# lines( lowess(x = drecent$date, y = pcscores$pca3, f = 1/20, iter = 10), col = "darkorange4", lwd=2)
# abline( v = warm_times, lty = 4, col = "slategrey" )
# axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
dev.off()

png(file="R/Figs/BC_Lightstation_PCA1_daily_timeseries.png", res = 600, width = 6, height = 1.5, units = "in")
par( mar=c(2,4,0,1)+0.1, mfrow=c(1,1), las = 1 )
plot(x = drecent$date, y = pcscores$pca1, type = 'l', col = 'slateblue', ylab = "PCA1" ); abline(h = 0)
lines( lowess(x = drecent$date, y = pcscores$pca1, f = 1/20, iter = 10), col = "darkslateblue", lwd=2)
abline( v = warm_times, lty = 4, col = "slategrey" )
axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
dev.off()

par( mar=c(2,4,0,1)+0.1, mfrow=c(1,1), las = 1 )
plot(x = drecent$date, y = pcscores$pca1, type = 'n', col = 'slateblue', ylab = "PCA1" ); abline(h = 0)
lines( lowess(x = drecent$date, y = pcscores$pca1, f = 1/20, iter = 10), col = "darkslateblue", lwd=2)
lines( lowess(x = drecent$date, y = pcscores$pca2, f = 1/20, iter = 10), col = "darkslateblue", lwd=2)
x <- lowess(x = drecent$date, y = pcscores$pca1, f = 1/20, iter = 10)
x <- do.call(cbind,x)
x <- as.data.frame(as.matrix(x))
names(x) <- c("x","pca1")
y <- lowess(x = drecent$date, y = pcscores$pca2, f = 1/20, iter = 10)
y <- do.call(cbind,y)
y <- as.data.frame(as.matrix(y))
names(y) <- c("x","pca2")
xy <- left_join(x, y)

#####
# PCA2 with temperature only:
# does this represent a more de-trended analysis that 
# represents anomalies within years (anomalous anomalies?)
library(viridis)
library(lubridate)
library(pals)
xy$date <- drecent$date
xy$doy <- yday(xy$date)
xy$year <- year(xy$date)
# ad survey dates
am <- read_csv("Data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )
surveys <- am %>% 
  filter(Site != "Meay Channel") %>% 
  filter(Year %in% 2012:2019) %>% 
  mutate(Date = ymd(Date)) %>% 
  select(Year,Date) %>% 
  distinct()
surveys <- xy[xy$date %in% surveys$Date,]
surveys <- surveys %>% 
  group_by(year) %>% 
  summarize(date=date[1], pca1=pca1[1], pca2 = pca2[2]) %>% 
  mutate(doy=yday(date))

ggplot( xy, aes(x = pca1, y = pca2, col=doy)) +
  geom_point(size = 0.5) +
  geom_point(data = surveys, aes(col=doy),size = 3, shape = 19 ) +
  geom_point(data = surveys, size = 3, shape = 21, col = 'black' ) +
  scale_color_gradientn(colours = pals::kovesi.cyclic_mygbm_30_95_c78_s25(365)) 

# PCA2 alone
dsurvey <- drecent
dsurvey$date <- ymd( dsurvey$date )
dsurvey$year <- year( dsurvey$date )
dsurvey$doy <- yday( dsurvey$date )
dsurvey <- data.frame(dsurvey, pcscores)
dsurvey <- dsurvey %>% 
  filter(year %in% 2011:2020)
ggplot(data = dsurvey, aes(x = date, y = pca2, col = doy) ) + 
  geom_hline(yintercept = 0 ) +
  geom_line() +
  geom_smooth() +
  geom_point(data = surveys,size = 3, shape = 19, col = "black" ) +
  scale_color_gradientn(colours = pals::kovesi.cyclic_mygbm_30_95_c78_s25(365)) 

#####
























dna <- bind_cols( drecent, pcscores )
dna$month = lubridate::month(dna$date)
dna$year = lubridate::year(dna$date)

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
ggplot(filter(anoms.annual, year>=2010 & survey.year <2020), aes(x=survey.year,y=pca1)) + geom_path() + geom_point()
ggplot(filter(anoms.annual, year>=2010), aes(x=survey.year,y=pca2)) + geom_path() + geom_point()
ggplot(filter(anoms.annual, year>=2010), aes(x=survey.year,y=pca3)) + geom_path() + geom_point()
ggplot(filter(anoms.annual), aes(x=survey.year,y=sal_mccinnis)) + geom_hline(yintercept = 0) +  geom_path(col = 'slateblue') + geom_point(col = 'slateblue') + ylim(c(-1.75,1.75))
ggplot(filter(anoms.annual), aes(x=survey.year,y=sal_pine)) + geom_hline(yintercept = 0) +  geom_path(col = 'slateblue') + geom_point(col = 'slateblue') + ylim(c(-1.75,1.75))
ggplot(filter(anoms.annual, year>=2010 &year<2020), aes(x=year,y=temp_pine)) + geom_path() + geom_point()


# extract data for 2011 to 2019
as.survey <- anoms.annual %>% 
  filter( survey.year>=2010 ) 


# write to disk
write_csv( as.survey, "Data/R code for Data Prep/Output from R/Lightstation_raw_PCA_impute.csv")
