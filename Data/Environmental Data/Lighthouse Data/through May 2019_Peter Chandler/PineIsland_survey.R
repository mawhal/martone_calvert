###### Calvert Island thermal environment data

## The purpose of this script is to examine temporal trends in temperature data
# this includes Pine Island Lighthouse SST, iButton data and tidbit data from Calvert Island, and 
# subtidal temperature data from Calvert

# by Matt Whalen
# created 16 July 2019

library(tidyverse)
library(lubridate)
library(fpp)
library(zoo)

# read data
draw <- read_table( file = "PineDailySalTemp.txt", skip=3 )

# renames columns, make date columns, replace 999.9 with NA
d <- draw %>%
  select( year=Year, month=Month, day=Day, temp=`Temperature(C)`, sal=`Salinity(psu)` ) %>%
  mutate( date=ymd(paste(year,month,day))) %>%
  mutate( temp=replace(temp, temp==999.9, NA)) %>%
  mutate( sal=replace(sal, sal==999.9, NA)) 

d[ is.na(d$temp),]

plot(temp~date,d,type='l')


# create a time series object
dts <- ts( d$temp, frequency=365, start=c(1937,1) )
plot( dts, type='l', col="dodgerblue" )

# impute data
dts.na <- na_ma( dts, k=4 )

# moving average -- does not work well with missing values
trend_temp = ma( dts.na, order = 2000, centre = T)
plot(dts)
lines(trend_temp, col="red", lw=3)
plot(as.ts(trend_temp))

# decompose seasons
decompose_temp = decompose(dts.na, "additive")




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
      col='dodgerblue', ylim=c(3,20.5) )
abline( h=0,lty=2 )
axis( 2, las=1 )
xtick <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
                  seq.Date(as.Date("2010-01-01"),as.Date("2019-01-01"), by="year" ) )
axis( 1,at=xtick,labels=2010:2019, srt=45,las=3 )
# get uniuue dates from surveys
meta <- read_csv( "../../../Martone_Hakai_metadata.csv")
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



## get data from iButtons all in one place
# for now, all from West Beach High Shore
# 2012+2013
ib13 <- read_csv( "../../iButtons/ibuttons/12_WBB_High_BlueFalcon.csv", skip=14 )
ib13$`Date/Time` <- as.POSIXlt(ib13$`Date/Time`, format="%d/%m/%Y %H:%M")
ib13$`Date/Time` <- ymd_hms(ib13$`Date/Time`)
# 2015
ib15 <- read_csv( "../../iButtons/ibuttons/2015 ibuttons/WBB_High_1A.csv", skip=14)
ib15$`Date/Time` <- as.POSIXlt(ib15$`Date/Time`, format="%m/%d/%y %I:%M:%S %p")
ib15$`Date/Time` <- ymd_hms(ib15$`Date/Time`)
# 2016
ib16 <- read_csv( "../../iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_High_N_BlackPVCCapW.csv",
                  skip=14 )
ib16$`Date/Time` <- as.POSIXlt(ib16$`Date/Time`, format="%m/%d/%y %I:%M:%S %p")
ib16$`Date/Time` <- ymd_hms(ib16$`Date/Time`)

# merge all of these
ib <- full_join( full_join( ib13,ib15),ib16 )
# posix
ib <- ib %>%
  select( date.time=`Date/Time`, value=Value ) 
# split date and time
split.df <- data.frame(do.call(rbind,strsplit(as.character(ib$date.time),split=" ")))
names(split.df) <- c('date','time')
ib <- data.frame( ib, split.df )
# summarize
ib.sum <- ib %>%
  mutate(date=ymd(date)) %>%
  group_by(date) %>%
  summarize(max=max(value),mean=mean(value))

#  plot it
ibs <- which( seq.Date(as.Date("2010-01-01"),as.Date("2014-01-01"), by="day" ) %in%
                           ib.sum$date )
ibs2 <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
                ib.sum$date[ib.sum$date>ymd("2014-01-01")] )

# lines( x=ibs, y=ib.sum$max, pch=20, col=rgb(0,0,0,0.2), cex=0.1  )
points( x=ibs, y=ib.sum$mean[ib.sum$date<ymd("2014-01-01")],lwd=0.5, col=rgb(0,0,0,0.5),pch=20,cex=0.2  )
points( x=ibs2, y=ib.sum$mean[ib.sum$date>ymd("2014-01-01")],lwd=0.5, col=rgb(0,0,0,0.5),pch=20,cex=0.2  )
  