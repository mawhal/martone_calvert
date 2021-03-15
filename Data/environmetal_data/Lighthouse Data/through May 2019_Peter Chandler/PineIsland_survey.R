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
library(imputeTS)
library(cowplot)

# read data
draw <- read_table( file = "Data/environmetal_data/Lighthouse Data/through May 2019_Peter Chandler/PineDailySalTemp.txt", skip=3 )

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
calendar2 <- time(dts.na)>=1980 & time(dts.na)<2000

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
      col='dodgerblue', ylim=c(0,30) )
abline( h=0,lty=2 )
axis( 2, las=1 )
xtick <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
                  seq.Date(as.Date("2010-01-01"),as.Date("2019-01-01"), by="year" ) )
axis( 1,at=xtick,labels=2010:2019, srt=45,las=3 )
# get uniuue dates from surveys
meta <- read_csv( "R/output/sampling_date_range.csv")
mean.date <- meta %>%
  pivot_longer( cols = c(start, end)) %>% 
  group_by( Year ) %>%
  summarise( date=mean(value) )
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
ib13 <- read_csv( "Data/environmetal_data/iButtons/ibuttons/12_WBB_High_BlueFalcon.csv", skip=14 )
ib13$`Date/Time` <- as.POSIXlt(ib13$`Date/Time`, format="%d/%m/%Y %H:%M")
ib13$`Date/Time` <- ymd_hms(ib13$`Date/Time`)
# 2015
ib15 <- read_csv( "Data/environmetal_data/iButtons/ibuttons/2015 ibuttons/WBB_High_1A.csv", skip=14)
ib15$`Date/Time` <- as.POSIXlt(ib15$`Date/Time`, format="%m/%d/%y %I:%M:%S %p")
ib15$`Date/Time` <- ymd_hms(ib15$`Date/Time`)
# 2016
ib16 <- read_csv( "Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_High_N_BlackPVCCapW.csv",
                  skip=14 )
ib16$`Date/Time` <- as.POSIXlt(ib16$`Date/Time`, format="%m/%d/%y %I:%M:%S %p")
ib16$`Date/Time` <- ymd_hms(ib16$`Date/Time`)

# merge all of these
ib <- full_join( full_join( ib13,ib15),ib16 )
write_csv( ib, "Data/environmetal_data/Lighthouse Data/through May 2019_Peter Chandler/output from R/FoggyCove_HIGH_iButton.csv")
# posix
ib <- ib %>%
  select( date.time=`Date/Time`, value=Value ) 
# split date and time
split.df <- data.frame(do.call(rbind,strsplit(as.character(ib$date.time),split=" ")))
names(split.df) <- c('date','time')
ib <- data.frame( ib, split.df )
ib <- ib %>% 
  mutate(date=ymd(date))
ib$lube.date <- ymd_hms(ib$date.time)
# summarize
ib.sum <- ib %>%
  group_by(date) %>%
  summarize(max=max(value),mean=mean(value))

#  plot it
dateseq <- data.frame( date=seq.Date(as.Date("2010-01-01"),as.Date("2020-01-01"), by="day") )
dateseq <- dateseq %>% mutate( rank=1:nrow(dateseq))
ib <- left_join( ib, dateseq )

# lines( x=ibs, y=ib.sum$max, pch=20, col=rgb(0,0,0,0.2), cex=0.1  )
points( x=ib$rank, y=ib$value,lwd=0.5, col=rgb(0,0,0,0.5),pch=20,cex=0.2  )
points( x=ibs, y=ib.sum$mean[ib.sum$date<ymd("2014-01-01")],lwd=0.5, col=rgb(0,0,0,0.5),pch=20,cex=0.2  )
points( x=ibs2, y=ib.sum$mean[ib.sum$date>ymd("2014-01-01")],lwd=0.5, col=rgb(0,0,0,0.5),pch=20,cex=0.2  )
  

svg(filename="Data/environmetal_data/Lighthouse Data/through May 2019_Peter Chandler/Figs/Pine_Foggy_compare.svg", 
    width=5, 
    height=3, 
    pointsize=12)

par(mar=c(5,4,1,2)+0.1)
plot( subset(dts, subset=calvert), type='n', axes=F,
      ylab=expression('Temperature ('*~degree*C*')'), xlab="Date", ylim=c(-5,40) )
points( x=ib$rank, y=ib$value,lwd=0.5, col=rgb(0,0,0,0.2),pch=20,cex=0.2  )
lines( subset(dts, subset=calvert),col='red' )
abline( h=0,lty=1 )
axis( 2, las=1 )
xtick <- which( seq.Date(as.Date("2010-01-01"),as.Date("2019-05-31"), by="day" ) %in%
                  seq.Date(as.Date("2010-01-01"),as.Date("2019-01-01"), by="year" ) )
axis( 1,at=xtick,labels=2010:2019, srt=45,las=3 )

abline( v=surveys, lty=2 )
# moving average
x.ma <- ma(dts.na,order=410, centre=T)
ma.calvert <- subset( x.ma, subset=calvert )
# lines( ma.calvert, col='red' )
# mean of entire time series
# abline( h=mean(dts,na.rm=T),col='blue' )
# abline( h=mean(subset(dts, subset=calvert),na.rm=T),col='goldenrod')
box()
dev.off()


# make date/time for pine island data
# ad a time each day, say noon
d$date.time <- ymd_hms(paste( d$date, "12:00:00"))
pt.size <- 0.6
fcib <- ggplot( data = filter(d,year >= 2010), aes(x=date.time,y=temp)) + 
  geom_point(data = ib, aes(x=lube.date, y=value), alpha=0.2, size=pt.size) +
  geom_line(col='red') +
  ylim(c(-3,37)) +
  ylab(expression(paste("Temperature (",degree,"C)"))) +
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

# just look at the second deployment 
ib2 <- ib %>% 
  filter( date > as.Date("2014-01-01"))
ib.month <- ib2 %>% 
  mutate( month=month(date), year=year(date) ) %>%
  mutate( year = ifelse(month==12,year+1,year)) 
ib.win <- ib.month %>% 
  filter( month %in% c(12,1,2) ) 
ib.sum <- ib.month %>% 
  filter( month %in% c(7,8,9) )

a <- ggplot(ib2, aes(x=date,y=value)) + geom_point(alpha=0.2) + geom_smooth() +
  xlab("Date") + ylab(expression(paste("Temperature (",degree,"C)"))) + ylim(c(-3,37))

b <- ggplot( ib.win, aes(group=year, x= factor(year), y=value)) + 
  geom_boxplot(notch=T, width=0.5, outlier.alpha = 0.2, outlier.size = pt.size, lwd=0.3) + 
    xlab("") + 
  ylab(expression(paste("Winter temperature (",degree,"C)"))) + 
  ylim(c(-3,37)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust=1))

c <- ggplot( ib.sum, aes(group=year, x= factor(year), y=value)) + 
  geom_boxplot(notch=T, width=0.5, outlier.alpha = 0.2, outlier.size = pt.size, lwd=0.3) + 
  xlab("") + 
  ylab(expression(paste("Summer temperature (",degree,"C)"))) + 
  ylim(c(-3,37)) + theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))


cowplot::plot_grid( fcib,c,b,ncol=3,rel_widths = c(4,1.25,1.25), labels = "AUTO" )
ggsave( "Data/environmetal_data/Lighthouse Data/through May 2019_Peter Chandler/Figs/ibutton_2014-2016.svg", width=6,height=3 )
