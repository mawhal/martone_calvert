# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Sam Starko 
# updated 8 Mar 2019

# This script is to summarize and visualize environmental data from Calvert and nearby stations

#load relavent packages
library(tidyverse)
library(lubridate)
library(pracma)
library(DataCombine)
library(imputeTS)
library(zoo)

# set y axis limits for anomaly plots
ylimits <- c(-3.5,3.5)

#################################################
########Nutrient Data from QCS01#################
################################################

#Summarize data by date and calculate average nitrogen for sites that have data, remove data that does not equal 5 m depth
nut<-read.csv("Data/environmental_data/CTD_nutrients_Hakai.csv")
nut.avg<-nut[nut$NO2.NO3..uM.!="NA"&nut$Line.Out.Depth==5,] %>% 
  group_by(Date, Line.Out.Depth) %>%
  summarize(avg.N = mean(NO2.NO3..uM.))
nut.avg<-nut.avg[complete.cases(nut.avg$Date),]

####Plot CTD Nitrogen through time
nut.avg$Date<-as.Date(nut.avg$Date)
nut.avg$Month<-month(as.POSIXlt(nut.avg$Date, format="%Y/%m/%d"))
nut.avg$Year<-year(as.POSIXlt(nut.avg$Date, format="%Y/%m/%d"))
plot(avg.N~Date, data=nut.avg, pch=19, cex=0.8, col="#296BA8", ylab="Total Nitrogen", las=1)
lines(nut.avg$Date,nut.avg$avg.N, lwd=3)

####Summarize by month and year
nut.avg2<-nut.avg %>%
  group_by(Month, Year) %>%
  summarize(avg.N = mean(avg.N))

plot(avg.N~Year, data=nut.avg2[nut.avg2$Month==9,])


#################################################
########Nutrient Data from Pruth#################
################################################

#Summarize data by date and calculate average nitrogen for sites that have data, remove data that does not equal 5 m depth
nut<-read.csv("Data/environmental_data/PRUTHCTD_2012-2018_HakaiData_nutrients.csv")
nut.avg<-nut[nut$NO2.NO3..uM.!="NA"&nut$Line.Out.Depth==5,] %>% 
  group_by(Date, Line.Out.Depth) %>%
  summarize(avg.N = mean(NO2.NO3..uM.))
nut.avg<-nut.avg[complete.cases(nut.avg$Date),]

####Plot CTD Nitrogen through time
nut.avg$Date<-as.Date(nut.avg$Date)
nut.avg$Month<-month(as.POSIXlt(nut.avg$Date, format="%Y/%m/%d"))
nut.avg$Year<-year(as.POSIXlt(nut.avg$Date, format="%Y/%m/%d"))
plot(avg.N~Date, data=nut.avg, pch=19, cex=0.8, col="#296BA8", ylab="Total Nitrogen", las=1)
lines(nut.avg$Date,nut.avg$avg.N, lwd=3)

####Summarize by month and year
nut.avg2<-nut.avg %>%
  group_by(Month, Year) %>%
  summarize(avg.N = mean(avg.N, na.rm=TRUE))

####Summarize by month
nut.avg3<-nut.avg2 %>%
  group_by(Month) %>%
  summarize(baseline.N = mean(avg.N, na.rm=TRUE))

###Combine using left_join and calculate anomalies
n.anom<-left_join(nut.avg2, nut.avg3, by="Month")
n.anom$Anomaly<-n.anom$avg.N-n.anom$baseline.N
n.anom$Date<-as.yearmon(paste(n.anom$Year, n.anom$Month), "%Y %m")
n.anom$Date<-as.Date(n.anom$Date)

#Plot nitrogen anomalies
plot(Anomaly~Date, data=n.anom[order(n.anom$Date),], pch=19, cex=0.8, col="#296BA8", ylab="Nitrogen Anomaly", las=1)
lines(n.anom$Anomaly[order(n.anom$Date)]~n.anom$Date[order(n.anom$Date)], lwd=3)


#################################################
############Temperature Data####################
###############################################

#CTD data is too sparse to say anything meaningful about water temperature
## try with data we already have 
temp.pine <- read_csv( "Data/environmetal_data/Lighthouse Data/2021_update/DATA_-_Active_Sites/Pine_Island/Pine_Island_-_Daily_Sea_Surface_Temperature_and_Salinity_1937-2021.csv",
                    skip = 1 )
temp.mc <- read_csv( "Data/environmetal_data/Lighthouse Data/2021_update/DATA_-_Active_Sites/McInnes_Island/McInnes_Island_-_Daily_Sea_Surface_Temperature_and_Salinity_1954-2021.csv",
                     skip = 1 )

# air temperature - Addenbrooke air temperature
adden <- read_csv( "data/environmetal_data/Addenbroke Air Temperature/EC/1060080.ascii", skip=1 )  # data from https://data.pacificclimate.org/portal/pcds/map/



#######SST data from Pine Island############
#temp.pine<-read.csv("./Data/environmetal_data/Lighthouse Data/through May 2019_Peter Chandler/PineDailySalTemp.txt", skip=2)
temp.pine 
#remove data from before 1985
temp.pine$temperature <-gsub("999.9","NA", temp.pine$`TEMPERATURE ( C )`)
temp.pine2<-temp.pine
#view temperature data
temp.pine2$temperature <- as.numeric(temp.pine2$temperature)
temp.pine2 <- temp.pine2 %>% select( data = `DATE (YYYY-MM-DD)`, temperature )

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
temp.pine.RockyTime<-temp.pine.comb[temp.pine.comb$Year>2010,]

#Add in months that have no data
temp.pine.RockyTime<-InsertRow(temp.pine.RockyTime, NewRow=c(2013,4,NA,8.50, NA), RowNum=28)
temp.pine.RockyTime<-InsertRow(temp.pine.RockyTime, NewRow=c(2015,1,NA,7.87, NA), RowNum=49)

#Interpolate those months

# temp.pine.RockyTime$Anomaly<-na.ma(temp.pine.RockyTime$Anomaly,k=12, weighting = "exponential")

#Assign colour to temperature anomaly direction
for (i in 1:length(temp.pine.RockyTime$Anomaly)) {
if ( all(c(!is.na(temp.pine.RockyTime$Anomaly[i]), temp.pine.RockyTime$Anomaly[i] > 0),na.rm=T) ) {
  temp.pine.RockyTime$An.col[i]<-"#B51D2C"
} else {
  temp.pine.RockyTime$An.col[i]<-"#296BA8"
}
}

temp.pine.RockyTime$Date<-as.yearmon(paste(temp.pine.RockyTime$Year, temp.pine.RockyTime$Month), "%Y %m")
temp.pine.RockyTime$Date<-as.Date(temp.pine.RockyTime$Date)

##Plot SST through time
par(mar=c(4,4,2,1))
par(mfrow=c(1,1))
barplot(temp.pine.RockyTime$Anomaly, col=temp.pine.RockyTime$An.col, ylim=c(-3,3), space=0, border=TRUE, las=1, width=1, xlim=c(0,12*8))
title(ylab=expression(paste("Temperature (", degree, "C)")), line=2, cex.lab=1.2)
axis(1, at=seq(0,12*8,12), labels = seq(2011,2019,1), srt=45, cex=2)
#lines(movavg(temp.pine.RockyTime$Anomaly, 4, type="s"),lwd=2,lty=1)
box()
####Export 3 x 5.5

# set x limits
xlimits <- range(temp.pine.RockyTime$Date)

yprange <- temp.pine.RockyTime %>% filter(Year > 2010 ) %>% ungroup() %>% select(Anomaly) %>% range(na.rm=T)
yp <- round(yprange)
diffp <- diff(c(yp[1],ylimits[1]))
png(filename="Data/environmetal_data/R Code and Analysis/figs/anomaly_pine.png", units="in", res=300, width=5, height=3, pointsize=15)
par(mar=c(3,5,1,1)+0.1,las=1, bg=NA, lend=2 )
plot(Anomaly~Date, data=temp.pine.RockyTime, col=temp.pine.RockyTime$An.col, type="n", axes=F,
     ylim=ylimits, xlim=xlimits,
     ylab = expression(paste("SST anomaly (", degree, "C)")), xlab="" )
# rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
points( Anomaly~Date, data=temp.pine.RockyTime, col=temp.pine.RockyTime$An.col, type="h", lwd=4 )
axis( 2, at = yp[1]:yp[2] )
axis.Date( x = ad.anon$Date, side = 1, line = 0+diffp, at = seq.Date(xlimits[1],xlimits[2],by="year") )
dev.off()

# write data to disk
write.csv( temp.pine.RockyTime, "Data/Environmental Data/PineIsland_anomaly.csv", row.names=FALSE )


#######SST data from McInnes Island############

####Note that data from this station is very patchy and several months worth of data are missing
#temp.mc<-read.csv("Data/Environmental Data/McInnesIsland_SST.csv")
#remove data from before 1985
temp.mc$`Temperature(C)`<-gsub("999.9","NA", temp.mc$`Temperature(C)`)
temp.mc2<-temp.mc

#view temperature data
temp.mc2$`Temperature(C)`<-as.numeric(temp.mc2$`Temperature(C)`)
colnames(temp.mc2)<-c("Year","Month","Day","Temperature","Salinity")


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
temp.mc.RockyTime<-temp.mc.comb[temp.mc.comb$Year>2010,]

#Interpolate those months

#temp.mc.RockyTime$Anomaly<-na.ma(temp.mc.RockyTime$Anomaly,k=12, weighting = "exponential")



#Assign colour to temperature anomaly direction
for (i in 1:length(temp.mc.RockyTime$Anomaly)) {
  if (temp.mc.RockyTime$Anomaly[i] > 0) {
    temp.mc.RockyTime$An.col[i]<-"#B51D2C"
  } else {
    temp.mc.RockyTime$An.col[i]<-"#296BA8"
  }
}

temp.mc.RockyTime$Date<-as.yearmon(paste(temp.mc.RockyTime$Year, temp.mc.RockyTime$Month), "%Y %m")
temp.mc.RockyTime$Date<-as.Date(temp.mc.RockyTime$Date)

temp.mc.narrow<-data.frame(temp.mc.RockyTime$Date, temp.mc.RockyTime$Temperature, temp.mc.RockyTime$Anomaly,temp.mc.RockyTime$An.col )
colnames(temp.mc.narrow)<-c("Date", "Temperature.mc", "Anomaly.mc", "An.col.mc")
temp.mc.narrow$An.col.mc <- as.character( temp.mc.narrow$An.col.mc )
temp.compare<-left_join(temp.pine.RockyTime,temp.mc.narrow, by="Date")


##Plot SST through time (note lots of missing data)
par(mar=c(4,4,2,1))
par(mfrow=c(1,1))
barplot(temp.compare$Anomaly.mc, col=c("#296BA8", "#B51D2C", "black")[temp.compare$An.col.mc], ylim=c(-3.5,3.5), space=0, border=TRUE, las=1, width=1, xlim=c(0,12*8))
title(ylab=expression(paste("Temperature (", degree, "C)")), line=2, cex.lab=1.2)
axis(1, at=c(0,12,24,36,48,12*5,12*6,12*7,12*8), labels = seq(2011,2019,1), srt=45)
#lines(movavg(temp.pine.RockyTime$Anomaly, 4, type="s"),lwd=2,lty=1)
box()
####Export 3 x 5.5

ymrange <- temp.compare %>% filter(Year > 2010 ) %>% ungroup() %>% select(Anomaly.mc) %>% range(na.rm=T)
ym <- round(ymrange)
diffm <- diff(c(ym[1],ylimits[1]))
png(filename="Data/environmetal_data/R Code and Analysis/figs/anomaly_mcinnis.png", units="in", res=300, width=5, height=3, pointsize=15)
par(mar=c(3,5,1,1)+0.1,las=1, bg=NA, lend=2 )
plot(Anomaly.mc~Date, data=temp.compare, col=temp.compare$An.col.mc, type="n", axes=F,
     ylim=ylimits, xlim=xlimits,
     ylab = expression(paste("SST anomaly (", degree, "C)")), xlab="" )
# rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
points( Anomaly.mc~Date, data=temp.compare, col=temp.compare$An.col.mc, type="h", lwd=4 )
axis( 2, at = ym[1]:ym[2] )
axis.Date( x = ad.anon$Date, side = 1, line = 0.5+diffm, at = seq.Date(xlimits[1],xlimits[2],by="year") )
dev.off()

##Average anomalies from both light stations 
temp.compare$Average <- (temp.compare$Temperature+temp.compare$Temperature.mc, na.rm = TRUE)


####################################
######Addenbroke Island Air Temp#####
######################################

#Import Data (Each year is in a separate .csv file)
ad.1985<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1985.csv", skip=24, stringsAsFactors = FALSE)
ad.1986<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1986.csv", skip=24, stringsAsFactors = FALSE)
ad.1987<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1987.csv", skip=24, stringsAsFactors = FALSE)
ad.1988<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1988.csv", skip=24, stringsAsFactors = FALSE)
ad.1989<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1989.csv", skip=24, stringsAsFactors = FALSE)
ad.1990<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1990.csv", skip=24, stringsAsFactors = FALSE)
ad.1991<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1991.csv", skip=24, stringsAsFactors = FALSE)
ad.1992<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1992.csv", skip=24, stringsAsFactors = FALSE)
ad.1993<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1993.csv", skip=24, stringsAsFactors = FALSE)
ad.1994<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1994.csv", skip=24, stringsAsFactors = FALSE)
ad.1995<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1995.csv", skip=24, stringsAsFactors = FALSE)
ad.1996<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1996.csv", skip=24, stringsAsFactors = FALSE)
ad.1997<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1997.csv", skip=24, stringsAsFactors = FALSE)
ad.1998<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1998.csv", skip=24, stringsAsFactors = FALSE)
ad.1999<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/1999.csv", skip=24, stringsAsFactors = FALSE)
ad.2000<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2000.csv", skip=24, stringsAsFactors = FALSE)
ad.2001<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2001.csv", skip=24, stringsAsFactors = FALSE)
ad.2002<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2002.csv", skip=24, stringsAsFactors = FALSE)
ad.2003<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2003.csv", skip=24, stringsAsFactors = FALSE)
ad.2004<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2004.csv", skip=24, stringsAsFactors = FALSE)
ad.2005<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2005.csv", skip=24, stringsAsFactors = FALSE)
ad.2006<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2006.csv", skip=24, stringsAsFactors = FALSE)
ad.2007<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2007.csv", skip=24, stringsAsFactors = FALSE)
ad.2008<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2008.csv", skip=24, stringsAsFactors = FALSE)
ad.2009<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2009.csv", skip=24, stringsAsFactors = FALSE)
ad.2010<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2010.csv", skip=24, stringsAsFactors = FALSE)
ad.2011<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2011.csv", skip=24, stringsAsFactors = FALSE)
ad.2012<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2012.csv", skip=24, stringsAsFactors = FALSE)
ad.2013<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2013.csv", skip=24, stringsAsFactors = FALSE)
ad.2014<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2014.csv", skip=24, stringsAsFactors = FALSE)
ad.2015<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2015.csv", skip=24, stringsAsFactors = FALSE)
ad.2016<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2016.csv", skip=24, stringsAsFactors = FALSE)
ad.2017<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2017.csv", skip=24, stringsAsFactors = FALSE)
ad.2018<-read.csv("Data/environmetal_data/Addenbroke Air Temperature/separate_years/2018.csv", skip=24, stringsAsFactors = FALSE)

#take only the first 10 columns, which contain temperature and date info
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
# ad.anon$Anomaly<-na.ma(ad.anon$Anomaly,k=12, weighting = "exponential")

#Create dummy variable that determines colour
for (i in 1:length(ad.anon$Anomaly)) {
  if ( all(c(!is.na(ad.anon$Anomaly[i]), ad.anon$Anomaly[i] > 0),na.rm=T) ) {
    ad.anon$An.col[i]<-"#B51D2C"
  } else {
    ad.anon$An.col[i]<-"#296BA8"
  }
}

ad.anon$Date<-as.yearmon(paste(ad.anon$Year, ad.anon$Month), "%Y %m")
ad.anon$Date<-as.Date(ad.anon$Date)


#Air Temp
barplot(ad.anon$Anomaly[ad.anon$Year>2010][1:115], col=ad.anon$An.col[ad.anon$Year>2010][1:115], ylim=c(-3.5,3.5), space=0, border=TRUE, las=1, width=1, xlim=c(0,12*8))
#lines(movavg(ad.anon$Anomaly[ad.anon$Year>2010][1:115], 12, type="s"),lwd=2,lty=1)
title(ylab=expression(paste("Temperature (", degree, "C)")), line=2, cex.lab=1.2)
axis(1, at=c(0,12,24,36,48,12*5,12*6,12*7,12*8), labels = seq(2011,2019,1), srt=45)
#lines(movavg(temp.pine.RockyTime$Anomaly, 4, type="s"),lwd=2,lty=1)
box()

# set custom locations for plotting axes
yarange <- ad.anon %>% filter(Year > 2010 ) %>% ungroup() %>% select(Anomaly) %>% range(na.rm=T)
ya <- round(yarange)
diffa <- diff(c(ya[1],ylimits[1]))
png(filename="Data/environmetal_data/R Code and Analysis/figs/anomaly_addenbrooke.png", units="in", res=300, width=5, height=3, pointsize=15)
par(mar=c(3,5,1,1)+0.1,las=1, bg=NA, lend=2, bty="n" )
plot(Anomaly~Date, data=ad.anon, col=ad.anon$An.col, type="n", 
     ylim=ylimits, xlim=xlimits,
     ylab = expression(paste("AT anomaly (", degree, "C)")), xlab="", axes=F )
# rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
points( Anomaly~Date, data=ad.anon, col=ad.anon$An.col, type="h", lwd=4 )
axis( 2, at = ya[1]:ya[2] )
axis.Date( x = ad.anon$Date, side = 1, line = 0+diffa, at = seq.Date(xlimits[1],xlimits[2],by="year") )
dev.off()

