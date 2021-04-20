###### BC Lighthouse Data on Sea Surface Temperature and Salinity
### NOTE: 14 day artifact as a result of sampling methods
### https://www.dfo-mpo.gc.ca/science/data-donnees/lightstations-phares/index-eng.html

## The purpose of this script is to identify the 14-day artifact, 
# to see if it can be removed through averaging or some other strategy
# AND to save monthly SST anomalies for downstream analyses

# by Matt Whalen
# created 16 July 2019
# updated 22 July 2020

library(tidyverse)
library(lubridate)
library(fpp)
library(imputeTS)

##### read data ####
# SST - daily from lightstation keepers
pineraw <- read_csv( "Data/environmetal_data/Lighthouse Data/2020 update/DATA_-_Active_Sites/Pine_Island/Pine_Island_-_Daily_Sea_Surface_Temperature_and_Salinity_1937-2019.csv",
                  skip=1 )
names(pineraw) <- c( 'date','sal','temp','latitude', 'longitude' )
pineraw$site <- 'pine'
mcinnraw <- read_csv( "Data/environmetal_data/Lighthouse Data/2020 update/DATA_-_Active_Sites/McInnes_Island/McInnes_Island_-_Daily_Sea_Surface_Temperature_and_Salinity_1954-2018.csv",
                  skip=1 )
names(mcinnraw) <- c( 'date','sal','temp','latitude', 'longitude' )
mcinnraw$site <- 'mccinnis'
# merge sst
sst <- bind_rows( pineraw, mcinnraw ) %>% 
  filter( !is.na(date) ) %>% 
  mutate( date = mdy(date) ) %>% 
  mutate( temp=replace(temp, temp==99.9, NA)) %>%
  mutate( sal=replace(sal, sal==99.9, NA)) 
# air temperature - Addenbrooke air temperature - hourly
addenraw <- read_csv( "data/environmetal_data/Addenbroke Air Temperature/EC/1060080.ascii", skip=1 )  # data from https://data.pacificclimate.org/portal/pcds/map/

names(addenraw) <- c( 'precip','rain','temp','snow','time','snow_ground','temp_max' )
addenraw$site <- 'addenbrooke'
adden <- addenraw %>% 
  mutate(date = ymd(time) ) %>% 
  select(-time)
# merge all data
d <- bind_rows(sst, adden)

ggplot( data = d, aes( x = date, y = temp)) +
  facet_wrap(~ site, ncol =1, scales = "free_y" ) +
  geom_line() + geom_smooth() 


## average by month and decompose the signal
dm <- d %>%
  mutate( year=year(date), month=month(date) ) %>% 
  group_by(year,month,site) %>%
  summarize( temp=mean(temp,na.rm=T),sal=mean(sal,na.rm=T),
             temp_max = mean(temp_max,na.rm=T),
             precip = mean(precip, na.rm=T) )

# # impute missing months
# dm$temp.na <- na_ma(dm$temp,k=12, weighting = "exponential")
# dm$sal.na <- na_ma(dm$sal,k=12, weighting = "exponential")
# dm[is.na(dm$temp),]
# dm[is.na(dm$sal),]

# calculate anomalies as deviations from expected (average) monthly mean tempearture
dm <- dm %>%
  group_by(month, site) %>%
  mutate( month.mean.temp = mean(temp,na.rm=T), month.mean.sal = mean(sal,na.rm=T),
          month.mean.temp.max = mean(temp_max,na.rm=T), month.mean.precip = mean(precip,na.rm=T)) %>%
  ungroup() %>%
  mutate( temp.anom = temp-month.mean.temp, sal.anom = sal - month.mean.sal,
          temp.max.anom = temp_max - month.mean.temp.max, precip.anom = precip - month.mean.precip )
## write tempeartures and anomalies to disk

write_csv( dm, "Data/R code for Data Prep/output from R/Lightstation_monthly_anomaly.csv" )





# 
# 
# ## spectral decomposition after seasons and trend removed
# dmts <- ts( dm$temp.na, frequency=12, start=c(1937,1) )
# stl_temp = stl(dmts, "periodic" )
# plot( stl_temp, col="dodgerblue" )
# # remove the seasons and trend
# dts.rem <- stl_temp$time.series[,3]
# # create the spectrum
# x.spec <- spectrum( dts.rem, log="no", plot=T )
# spx <- x.spec$freq
# spy <- 2*x.spec$spec
# plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(0,6))
# 
# # add in colors for sign of anomalies
# dm <- dm %>%
#   mutate( temp.col=ifelse(temp.anom>0,"red","blue") )
# 
# # make a ts object
# dmts <- ts( dm$temp.anom, frequency=12, start=c(1937,1) )
# dmts.col <- ts( dm$temp.col, frequency=12, start=c(1937,1) )
# windows(8,3)
# par(mar=c(5,4,2,2)+0.1)
# plot(dmts,type="n")
# points(dmts,col=dmts.col,pch=20)
# segments(x0=time(dmts),y0=0,x1=time(dmts),y1=dmts,col=dmts.col)
# lines(ma(dmts, 12),lwd=2,lty=1, col='lightblue')
# 
# 
# 
# # time period of Patrick and Sandra's dataset
# calvert <- time(dmts)>=2010 
# par(mar=c(5,4,1,2)+0.1)
# plot( subset(dmts, subset=calvert), type='n', axes=F,
#       ylab=expression('Sea Surface Temperature ('*~degree*C*')'), xlab="Date",
#       col='dodgerblue' )
# abline( h=0,lty=2 )
# axis( 2, las=1 )
# xtick <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
#                   seq.Date(as.Date("2010-01-01"),as.Date("2019-01-01"), by="year" ) )
# axis( 1,at=xtick,labels=2010:2019, srt=45,las=3 )
# # get uniuue dates from surveys
# meta <- read_csv( "../../../Martone_Hakai_metadata.csv")
# mean.date <- meta %>%
#   mutate( year=year(Date) ) %>%
#   group_by( year ) %>%
#   summarise( date=mean(Date) )
# surveys <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
#                     ymd(mean.date$date) )
# abline( v=surveys, lty=2 )
# # moving average
# x.ma <- ma(dts.na,order=410, centre=T)
# ma.calvert <- subset( x.ma, subset=calvert )
# lines( ma.calvert, col='red' )
# 






# create a time series object
dts <- lapply( split(d,d$site), function(z)
  ts(z$temp, frequency = 365, start = c(year(min(z$date)),1) ) )
dts2 <- lapply( split(d,d$site), function(z)
  ts(z$sal, frequency = 365, start = c(year(min(z$date)),1) ) )
# add max air temp, precip, rain?
dts[[4]] <- ts(adden$temp_max, frequency = 365, start = c(year(min(adden$date)),1) )
dts[[5]] <- ts(adden$precip, frequency = 365, start = c(year(min(adden$date)),1) )
dts[[6]] <- dts2[[2]]
dts[[7]] <- dts2[[3]]
par( mfrow = c(5,1), mar = c(3,4,0,2)+0.1 )
plot( dts[[1]], type='l', col="dodgerblue" )
plot( dts[[2]], type='l', col="royalblue" )
plot( dts[[3]], type='l', col="midnightblue" )
plot( dts[[4]], type='l', col="black" )
plot( dts[[5]], type='l', col="orange" )


# impute data
dts.na <- lapply( dts, function(z) na_ma( z, k=4 ) )

# moving average -- does not work well with missing values
line_cols <- c("dodgerblue","royalblue","midnightblue","black","orange","pink","magenta")
trend_temp = lapply( dts.na, function(z) ma( z, order = 1000, centre = T) )

par( mfrow = c(length(dts),1), mar = c(3,4,0,2)+0.1 )
for(i in 1:length(dts)){
  plot(dts[[i]], col = line_cols[i])
  lines(trend_temp[[i]], col="red", lw=2.5)
}
# decompose seasons
decompose_temp = lapply( dts.na, function(z) decompose(z, "additive") )
plot( decompose_temp[[1]] )
plot( decompose_temp[[7]] )

# # stl - LOESS version
# stl_temp = stl(dts.na, "periodic" )
# plot( stl_temp, col="dodgerblue" )

# examine the stereotypical seasonal pattern
season <- time(dts.na)>1940 & time(dts.na)<1945
plot( decompose_temp$seasonal )
plot( subset(decompose_temp$seasonal, subset=season), type='l' )
calendar <- time(dts.na)>=1940 & time(dts.na)<1941

windows(6,4)
par(mar=c(6,4,1,2)+0.1)
plot( subset(decompose_temp$seasonal, subset=calendar), type='l', axes=F,
      ylab=expression('SST Deviation ('*~degree*C*')'), xlab="",
      col='dodgerblue' )
abline( h=0,lty=2 )
axis( 2, las=1 )
xtick <- which( seq.Date(as.Date("1940-01-01"),as.Date("1941-12-31"), by="day" ) %in% 
  seq.Date(as.Date("1940-01-01"),as.Date("1940-12-31"), by="month" ) )
axis( 1,at=xtick,labels=month.name, srt=45,las=3 )
mtext("Date",1,4)

# time period of Patrick and Sandra's dataset
calvert <- time(dts.na)>=2010 
par(mar=c(5,4,1,2)+0.1)
plot( subset(dts, subset=calvert), type='l', axes=F,
      ylab=expression('Sea Surface Temperature ('*~degree*C*')'), xlab="Date",
      col='dodgerblue' )
abline( h=0,lty=2 )
axis( 2, las=1 )
xtick <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
                  seq.Date(as.Date("2010-01-01"),as.Date("2019-01-01"), by="year" ) )
axis( 1,at=xtick,labels=2010:2019, srt=45,las=3 )
# get unique dates from surveys
meta <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv")
mean.date <- meta %>%
  mutate( year=year(Date) ) %>%
  group_by( year ) %>%
  summarise( date=mean(Date) )
surveys <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
                    ymd(mean.date$date) )
abline( v=surveys, lty=2 )
# moving average
x.ma <- ma(dts.na,order=410, centre=T)
ma.calvert <- subset( x.ma, subset=calvert )
lines( ma.calvert, col='red' )
# mean of entire time series
abline( h=mean(dts,na.rm=T),col='blue' )
abline( h=mean(subset(dts, subset=calvert),na.rm=T),col='goldenrod')
box()

# remove the seasons and trend
dts.rem <- stl_temp$time.series[,3]


# simple Fourier transform
x.k <- fft(dts.rem)


# create the spectrum
x.spec <- spectrum( dts.rem, log="no", plot=T )
spx <- x.spec$freq
spy <- 2*x.spec$spec
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(0,100))
# there should be an ~14 day signal here 
# frequency here defined per year 
spectrum(dts.na,log="no",xlim=c(0,10))
# so, 14 days should be at frequency 26 (365/14)
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(20,30))
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(25.7,25.75))  # 25.72 == 14.19129 days
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(23.7,23.75))  # 23.72 == 15.38786 days
# yes. we see signal every 14 to 16 days, but this is relatively small (peaks of 0.1 compared to 0.35 for yearly signal)
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(13,15))       # 13.25 == 27.5 days, 14.25 == 25.6 days
# these additional peaks in the spectrum may have to do with spring tides
# offsets due to leap years?
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(0,26))
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l",xlim=c(0,5))


