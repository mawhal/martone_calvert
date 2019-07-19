# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Sam Starko and Matt Whalen

# This script is to summarize and visualize environmental data from Calvert and nearby stations

#load relavent packages
library(tidyverse)
library(lubridate)
# library(pracma)
# library(DataCombine)
library(imputeTS)
library(fpp)
library(zoo)


#### NOTE ABOUT THE DATA FROM LIGHTHOUSES
# :::WARNING:::
# The daily sampling strategy at the BC Lighthouse Stations was designed long ago 
# by Dr. John P. Tully. We have chosen not to change the strategy in the interests 
# of a homogeneous data set. Sampling occurs at or near the daytime high tide. 
# This means, for example, that if an observer starts sampling one day at 6 a.m.
# (local time), and continues to sample at the daytime hightide, as instructed, 
# then on the 2nd day he/she will take samples at about 06:50 the next day, 
# 07:40 the day after etc. When the daytime high tide gets close to 6 p.m. 
# then it snaps back to 6 a.m. and the cycle starts again. Since there is a 
# diurnal signal in sea-surface temperature the sampling creates a 14-day signal 
# as an artifact
##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Do we need to think about how to remove this artifact? Probably not since 
#### the temporal scale of interest here is months and years
#### See scripts in "Environental Data/Lighthouse Data/through May 2019 Peter Chandler",
####  one of which shows that the spectral signal of the 14-day-ish artifact is relatively small
####  it should go away when averaged across months



#################################################
############Temperature Data####################
###############################################

### CTD data is too sparse to say anything meaningful about water temperature,
### but we also have iButton and tidbit data that can help

### BC Lighthouse Data -- the best we have for a climatological signal
### cross-correlate with sparser intertidal temperature loggers and CTD casts

# read data
draw <- read_table( file = "../Lighthouse Data/through May 2019_Peter Chandler/PineDailySalTemp.txt", skip=3 )

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


###Average by month
temp.pine.sum<-temp.pine2[temp.pine2$Temperature!="NA",] %>%
  group_by(Year, Month) %>%
  summarize(Temperature = mean(Temperature))

temp.pine.byMonth <- temp.pine.sum %>%
  group_by(Month) %>%
  summarize(avgTemperature = mean(Temperature))

temp.pine.byMonth<-temp.pine.byMonth[1:12,]
temp.pine.byMonth

##Calculate Anomalies
temp.pine.comb<-left_join(temp.pine.sum, temp.pine.byMonth, by="Month")
temp.pine.comb$Anomaly <- temp.pine.comb$Temperature - temp.pine.comb$avgTemperature
temp.pine.comb<-temp.pine.comb[complete.cases(temp.pine.comb$Anomaly),]
temp.pine.RockyTime<-temp.pine.comb[temp.pine.comb$Year>2008,]

#Add in months that have no data
temp.pine.RockyTime<-InsertRow(temp.pine.RockyTime, NewRow=c(2013,4,NA,8.50, NA), RowNum=28)
temp.pine.RockyTime<-InsertRow(temp.pine.RockyTime, NewRow=c(2015,1,NA,7.87, NA), RowNum=49)

#Interpolate those months

temp.pine.RockyTime$Anomaly<-na.ma(temp.pine.RockyTime$Anomaly,k=12, weighting = "exponential")

#Assign colour to temperature anomaly direction
for (i in 1:length(temp.pine.RockyTime$Anomaly)) {
if (temp.pine.RockyTime$Anomaly[i] > 0) {
  temp.pine.RockyTime$An.col[i]<-"red"
} else {
  temp.pine.RockyTime$An.col[i]<-"blue"
}
}

temp.pine.RockyTime$Date<-as.yearmon(paste(temp.pine.RockyTime$Year, temp.pine.RockyTime$Month), "%Y %m")
temp.pine.RockyTime$Date<-as.Date(temp.pine.RockyTime$Date)

##Plot SST through time
par(mar=c(4,4,2,1))
par(mfrow=c(1,1))
barplot(temp.pine.RockyTime$Anomaly, col=temp.pine.RockyTime$An.col, ylim=c(-2.5,2.5), space=0, border=FALSE, las=1, width=1)
axis(1, at=c(0,12,24,36,48,12*5,12*6,12*7, 12*8, 12*9))
lines(movavg(temp.pine.RockyTime$Anomaly, 12, type="s"),lwd=2,lty=1)


# write data to disk
write.csv( temp.pine.RockyTime, "Data/Environmental Data/PineIsland_anomaly.csv", row.names=FALSE )


#######SST data from McInnes Island############

####Note that data from this station is very patchy and several months worth of data are missing
temp.mc<-read.csv("Data/Environmental Data/McInnesIsland_SST.csv")
#remove data from before 1985
temp.mc2<-temp.mc[temp.mc$Year>1984,]
#view temperature data
temp.mc2$Temperature

###Average by month
temp.mc.sum<-temp.mc2[temp.mc2$Temperature!="NA",] %>%
  group_by(Year, Month) %>%
  summarize(Temperature = mean(Temperature))

temp.mc.byMonth <- temp.mc.sum %>%
  group_by(Month) %>%
  summarize(avgTemperature = mean(Temperature))

temp.mc.byMonth<-temp.mc.byMonth[1:12,]
temp.mc.byMonth

##Calculate Anomalies
temp.mc.comb<-left_join(temp.mc.sum, temp.mc.byMonth, by="Month")
temp.mc.comb$Anomaly <- temp.mc.comb$Temperature - temp.mc.comb$avgTemperature
temp.mc.comb<-temp.mc.comb[complete.cases(temp.mc.comb$Anomaly),]
temp.mc.RockyTime<-temp.mc.comb[temp.mc.comb$Year>2009,]

#Interpolate those months

temp.mc.RockyTime$Anomaly<-na.ma(temp.mc.RockyTime$Anomaly,k=12, weighting = "exponential")

#Assign colour to temperature anomaly direction
for (i in 1:length(temp.mc.RockyTime$Anomaly)) {
  if (temp.mc.RockyTime$Anomaly[i] > 0) {
    temp.mc.RockyTime$An.col[i]<-"red"
  } else {
    temp.mc.RockyTime$An.col[i]<-"blue"
  }
}

temp.mc.RockyTime$Date<-as.yearmon(paste(temp.mc.RockyTime$Year, temp.mc.RockyTime$Month), "%Y %m")
temp.mc.RockyTime$Date<-as.Date(temp.mc.RockyTime$Date)

temp.mc.narrow<-data.frame(temp.mc.RockyTime$Date, temp.mc.RockyTime$Temperature, temp.mc.RockyTime$Anomaly,temp.mc.RockyTime$An.col )
colnames(temp.mc.narrow)<-c("Date", "Temperature.mc", "Anomaly.mc", "An.col.mc")
temp.compare<-left_join(temp.pine.RockyTime,temp.mc.narrow, by="Date")


##Plot SST through time (note lots of missing data)
par(mar=c(4,4,2,1))
par(mfrow=c(1,1))
barplot(temp.compare$Anomaly.mc, col=c("blue", "red", "black")[temp.compare$An.col.mc], ylim=c(-3.5,3.5), space=0, border=FALSE, las=1, width=1)
axis(1, at=c(0,12,24,36,48,12*5,12*6,12*7, 12*8, 12*9))


##Average anomalies from both light stations 
temp.compare$Average <- (temp.compare$Temperature+temp.compare$Temperature.mc, na.rm = TRUE)


####################################
######Addenbroke Island Air Temp#####
######################################

#Import Data (Each year is in a separate .csv file)
ad.1985<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1985.csv", skip=24, stringsAsFactors = FALSE)
ad.1986<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1986.csv", skip=24, stringsAsFactors = FALSE)
ad.1987<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1987.csv", skip=24, stringsAsFactors = FALSE)
ad.1988<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1988.csv", skip=24, stringsAsFactors = FALSE)
ad.1989<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1989.csv", skip=24, stringsAsFactors = FALSE)
ad.1990<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1990.csv", skip=24, stringsAsFactors = FALSE)
ad.1991<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1991.csv", skip=24, stringsAsFactors = FALSE)
ad.1992<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1992.csv", skip=24, stringsAsFactors = FALSE)
ad.1993<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1993.csv", skip=24, stringsAsFactors = FALSE)
ad.1994<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1994.csv", skip=24, stringsAsFactors = FALSE)
ad.1995<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1995.csv", skip=24, stringsAsFactors = FALSE)
ad.1996<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1996.csv", skip=24, stringsAsFactors = FALSE)
ad.1997<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1997.csv", skip=24, stringsAsFactors = FALSE)
ad.1998<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1998.csv", skip=24, stringsAsFactors = FALSE)
ad.1999<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/1999.csv", skip=24, stringsAsFactors = FALSE)
ad.2000<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2000.csv", skip=24, stringsAsFactors = FALSE)
ad.2001<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2001.csv", skip=24, stringsAsFactors = FALSE)
ad.2002<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2002.csv", skip=24, stringsAsFactors = FALSE)
ad.2003<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2003.csv", skip=24, stringsAsFactors = FALSE)
ad.2004<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2004.csv", skip=24, stringsAsFactors = FALSE)
ad.2005<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2005.csv", skip=24, stringsAsFactors = FALSE)
ad.2006<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2006.csv", skip=24, stringsAsFactors = FALSE)
ad.2007<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2007.csv", skip=24, stringsAsFactors = FALSE)
ad.2008<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2008.csv", skip=24, stringsAsFactors = FALSE)
ad.2009<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2009.csv", skip=24, stringsAsFactors = FALSE)
ad.2010<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2010.csv", skip=24, stringsAsFactors = FALSE)
ad.2011<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2011.csv", skip=24, stringsAsFactors = FALSE)
ad.2012<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2012.csv", skip=24, stringsAsFactors = FALSE)
ad.2013<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2013.csv", skip=24, stringsAsFactors = FALSE)
ad.2014<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2014.csv", skip=24, stringsAsFactors = FALSE)
ad.2015<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2015.csv", skip=24, stringsAsFactors = FALSE)
ad.2016<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2016.csv", skip=24, stringsAsFactors = FALSE)
ad.2017<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2017.csv", skip=24, stringsAsFactors = FALSE)
ad.2018<-read.csv("Data/Environmental Data/Addenbroke Air Temperature/2018.csv", skip=24, stringsAsFactors = FALSE)

#take only the first 10 rows, which contain temperature and date info
ad.1985<-ad.1985[,1:10]
ad.1986<-ad.1986[,1:10]
ad.1987<-ad.1987[,1:10]
ad.1988<-ad.1988[,1:10]
ad.1989<-ad.1989[,1:10]
ad.1990<-ad.1990[,1:10]
ad.1991<-ad.1991[,1:10]
ad.1992<-ad.1992[,1:10]
ad.1993<-ad.1993[,1:10]
ad.1994<-ad.1994[,1:10]
ad.1995<-ad.1995[,1:10]
ad.1996<-ad.1996[,1:10]
ad.1997<-ad.1997[,1:10]
ad.1998<-ad.1998[,1:10]
ad.1999<-ad.1999[,1:10]
ad.2000<-ad.2000[,1:10]
ad.2001<-ad.2001[,1:10]
ad.2002<-ad.2002[,1:10]
ad.2003<-ad.2003[,1:10]
ad.2004<-ad.2004[,1:10]
ad.2005<-ad.2005[,1:10]
ad.2006<-ad.2006[,1:10]
ad.2007<-ad.2007[,1:10]
ad.2008<-ad.2008[,1:10]
ad.2009<-ad.2009[,1:10]
ad.2010<-ad.2010[,1:10]
ad.2011<-ad.2011[,1:10]
ad.2012<-ad.2012[,1:10]
ad.2013<-ad.2013[,1:10]
ad.2014<-ad.2014[,1:10]
ad.2015<-ad.2015[,1:10]
ad.2016<-ad.2016[,1:10]
ad.2017<-ad.2017[,1:10]
ad.2018<-ad.2018[,1:10]


#Combine data from all of the years into a single data frame
ad.comb<-bind_rows(ad.1985, ad.1986, ad.1987, ad.1988, ad.1989, ad.1990, ad.1991,
          ad.1992, ad.1993, ad.1994, ad.1995, ad.1996, ad.1997, ad.1998,
          ad.1999, ad.2000, ad.2001, ad.2002, ad.2003, ad.2004,ad.2005,
          ad.2006, ad.2007, ad.2008, ad.2009, ad.2010, ad.2011, ad.2012,
          ad.2013, ad.2014, ad.2015, ad.2016, ad.2017, ad.2018)

#rename the columns to make them easier to call
names(ad.comb)<-c("Date", "Year", "Month", "Day", "Data.Quality","MaxTemp", "MaxTempFlag", "MinTemp", "MinTempFlag", "MeanTemp")

#Average MeanTemp data by year and by month
ad.summ<-ad.comb %>%
  group_by(Year, Month) %>%
  summarize(AvgTemp=mean(MeanTemp, na.rm=TRUE))

#Average MeanTemp data by month only
ad.avg<-ad.comb %>%
  group_by(Month) %>%
  summarize(AvgTemp=mean(MeanTemp, na.rm=TRUE))

#Combine dataframes with left_join and calculate anomalies
ad.anon<-left_join(ad.summ, ad.avg, by="Month")
names(ad.anon)<-c("Year","Month", "Temp", "BaselineTemp")
ad.anon$Anomaly <- ad.anon$Temp - ad.anon$BaselineTemp

#Interpolate missing values for three months
ad.anon$Anomaly<-na.ma(ad.anon$Anomaly,k=12, weighting = "exponential")

#Create dummy variable that determines colour
for (i in 1:length(ad.anon$Anomaly)) {
  if (ad.anon$Anomaly[i] > 0) {
    ad.anon$An.col[i]<-"red"
  } else {
    ad.anon$An.col[i]<-"blue"
  }
}



#Plot data
##Plot SST through time
par(mar=c(4,3,0,1))
par(mfrow=c(2,1))
barplot(temp.pine.RockyTime$Anomaly[temp.pine.RockyTime$Year>2008][1:115], col=temp.pine.RockyTime$An.col[temp.pine.RockyTime$Year>2008][1:115], ylim=c(-2.5,2.5), space=0, border=FALSE, las=1, width=1)
axis(1, at=c(0,12,24,36,48,12*5,12*6,12*7, 12*8, 12*9))
lines(movavg(temp.pine.RockyTime$Anomaly[temp.pine.RockyTime$Year>2008][1:115], 12, type="s"),lwd=2,lty=1)

#Air Temp
barplot(ad.anon$Anomaly[ad.anon$Year>2008][1:115], col=ad.anon$An.col[ad.anon$Year>2008][1:115], ylim=c(-3.5,3.5), space=0, border=FALSE, las=1, width=1)
lines(movavg(ad.anon$Anomaly[ad.anon$Year>2008][1:115], 12, type="s"),lwd=2,lty=1)
axis(1, at=c(0,12,24,36,48,12*5,12*6,12*7, 12*8, 12*9))

