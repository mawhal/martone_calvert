# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Sam Starko 
# updated 8 Mar 2019

# This script is to summarize and visualize environmental data from Calvert and nearby stations

#load relavent packages
library(tidyverse)
library(lubridate)
library(pracma)

#################################################
########Nutrient Data from QCS01#################
################################################

#Summarize data by date and calculate average nitrogen for sites that have data, remove data that does not equal 5 m depth
nut<-read.csv("./Data/Environmental Data/CTD_nutrients_Hakai.csv")
nut.avg<-nut[nut$NO2.NO3..uM.!="NA"&nut$Line.Out.Depth==5,] %>% 
  group_by(Date, Line.Out.Depth) %>%
  summarize(avg.N = mean(NO2.NO3..uM.))
nut.avg<-nut.avg[complete.cases(nut.avg$Date),]

####Plot CTD Nitrogen through time
nut.avg$Date<-as.Date(nut.avg$Date)
nut.avg$Month<-month(as.POSIXlt(nut.avg$Date, format="%Y/%m/%d"))
nut.avg$Year<-year(as.POSIXlt(nut.avg$Date, format="%Y/%m/%d"))
plot(avg.N~Date, data=nut.avg, pch=19, cex=0.8, col="Blue", ylab="Total Nitrogen", las=1)
lines(nut.avg$Date,nut.avg$avg.N, lwd=3)

####Summarize by month and year
nut.avg2<-nut.avg %>%
  group_by(Month, Year) %>%
  summarize(avg.N = mean(avg.N))

plot(avg.N~Year, data=nut.avg2[nut.avg2$Month==9,])


#################################################
############Temperature Data####################
###############################################

#CTD data is too sparse to say anything meaningful about water temperature

#######SST data from Pine Island############
temp.pine<-read.csv("./Data/Environmental Data/PineIsland_SST.csv")
#remove data from before 1985
temp.pine2<-temp.pine[temp.pine$Year>1984,]
#view temperature data
temp.pine2$Temperature

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

#Assign colour to temperature anomaly direction
for (i in 1:length(temp.pine.comb$Anomaly)) {
if (temp.pine.comb$Anomaly[i] > 0) {
  temp.pine.comb$An.col[i]<-"red"
} else {
  temp.pine.comb$An.col[i]<-"blue"
}
}

##Plot SST through time
par(mar=c(4,4,2,1))
par(mfrow=c(1,1))
barplot(temp.pine.comb$Anomaly[temp.pine.comb$Year > 1989], col=temp.pine.comb$An.col[temp.pine.comb$Year > 1989], ylim=c(-3.2,3.2), space=0, border=FALSE, las=1, width=1)
axis(1, at=c(0,48,48*2,48*3,48*4,48*5,48*6,48*7))
lines(movavg(temp.pine.comb$Anomaly[temp.pine.comb$Year > 1989], 12, type="s"),lwd=2,lty=1)





