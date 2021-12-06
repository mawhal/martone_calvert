### Hakai Institute  weather station data
# data from the Lookout on Calvert Island
# Whalen

library(tidyverse)
library(lubridate)

# read data
draw <- read_csv( "Data/environmental_data/calvert_lookout/2020-04-21.1hourSamples.csv" )
draw <- draw[-c(1:3),]
names(draw)
draw2 <- read_csv( "Data/environmental_data/calvert_lookout/2020-05-04.1hourSamples.csv" )
draw2 <- draw2[-c(1:3),]
names(draw2)

dmerge <- full_join( draw, draw2 )

# select columns
dsel <- dmerge %>% select( date_time = `Measurement time`, Year, Month, WaterYear, 
                 temp1q=TAirLookout1_Q_level,temp1mean=TAirLookout1_Avg,temp2q=TAirLookout2_Q_level,temp2mean=TAirLookout2_Avg,
                 rhq=RHLookoutQ_level,rhmean=RHLookoutAvg, 
                 SolarRadQ_level, SolarMax=SolarRadMax,
                 PAR_Q_level, PARmax=PAR_Max )
unique(dsel$temp1q)
unique(dsel$temp2q)
unique(dsel$SolarRadQ_level)
unique(dsel$PAR_Q_level)


# dates, numeric
d <- dsel %>% mutate( date_time=ymd_hms(date_time), date=date(date_time) ) %>% 
  mutate( temp1mean=as.numeric(temp1mean),
          temp2mean=as.numeric(temp2mean),
          rhmean=as.numeric(rhmean),
          SolarMax=as.numeric(SolarMax),
          PARmax=as.numeric(PARmax))

# daily data
dday <- d %>% 
  group_by( date_time, Year, date ) %>% 
  mutate( temp = mean(c(temp1mean,temp2mean),na.rm=T)) %>% 
  group_by( Year, date ) %>% 
  summarize( temp_mean=mean(temp,na.rm=T),temp_max = max(temp,na.rm=T), rhmean=mean(rhmean,na.rm=T),
             SolarMax=mean(SolarMax,na.rm=T), PARmax=mean(PARmax,na.rm=T)  ) %>% 
  arrange( date )

# use heatwaveR to get seasonal trend
library(heatwaveR)
# select tempearture only for now
d <- dday %>%
  ungroup() %>% 
  select( t=date, temp=temp_mean ) 
ts <- ts2clm(d, climatologyPeriod = range(d$t) )
ts10 <- ts2clm(d, climatologyPeriod = range(d$t), pctile=10)
ts_max <- ts2clm( select(ungroup(dday), t=date,temp=temp_max), climatologyPeriod = range(d$t))

a <- ggplot( ts, aes(x=t, y=temp, col=temp-seas) ) + 

  geom_point( alpha=1 ) +
  geom_line( aes(y=seas),col='gray50' ) +
  geom_line( aes(y=thresh),col='gray75' ) +
  geom_line( data=ts10, aes(y=thresh),col='gray75' ) +
  geom_line( data=ts_max, aes( x=date, y=temp_max))  +
  guides(col=guide_colorbar(title="seasonal\nanomaly")) +
  scale_color_gradient2(low = "#296BA8",mid="whitesmoke", high = "#B51D2C") + 
  ylab( expression(paste("Temperature (",degree,"C)")) ) + xlab("")
ggsave( "Data/environmetal_data/calvert_lookout_airtemp_humidity/lookout_daily.svg", width=6,height=2 )


# humidity
d <- dday %>%
  ungroup() %>% 
  select( t=date, temp=rhmean ) 
ts <- ts2clm(d, climatologyPeriod = range(d$t) )
ts10 <- ts2clm(d, climatologyPeriod = range(d$t), pctile=10)

b <- ggplot( ts, aes(x=t, y=temp, col=temp-seas) ) + 
  geom_point( alpha=1 ) +
  geom_line( aes(y=seas),col='gray50' ) +
  geom_line( aes(y=thresh),col='gray75' ) +
  geom_line( data=ts10, aes(y=thresh),col='gray75' ) +
  guides(col=guide_colorbar(title="seasonal\nanomaly")) +
  scale_color_gradient2(low = "#296BA8",mid="whitesmoke", high = "#B51D2C") + 
  ylab( expression(paste("Relative Humidity (%)")) ) + xlab("")
cowplot::plot_grid(a,b, ncol=1)
ggsave( "Data/environmetal_data/calvert_lookout_airtemp_humidity/lookout_daily.svg", width=6,height=4 )



# PAR
d <- dday %>%
  ungroup() %>% 
  select( t=date, temp=PARmax ) 
ts <- ts2clm(d, climatologyPeriod = range(d$t) )
ts10 <- ts2clm(d, climatologyPeriod = range(d$t), pctile=10)

c <- ggplot( ts, aes(x=t, y=temp, col=temp-seas) ) + 
  geom_point( alpha=1 ) +
  geom_line( aes(y=seas),col='gray50' ) +
  geom_line( aes(y=thresh),col='gray75' ) +
  geom_line( data=ts10, aes(y=thresh),col='gray75' ) +
  guides(col=guide_colorbar(title="seasonal\nanomaly")) +
  scale_color_gradient2(low = "#296BA8",mid="whitesmoke", high = "#B51D2C") + 
  ylab( expression(paste("PARmax")) ) + xlab("")
cowplot::plot_grid(a,c, ncol=1,align = "hv")
ggsave( "Data/environmetal_data/calvert_lookout_airtemp_humidity/lookout_daily.svg", width=6,height=4 )



# Cross-correlation
dcf <- dday %>%
  ungroup() %>% 
  select( t=date, temp, PARmax )
dcf <- na.omit(dcf)
x <- ccf( dcf$PARmax, dcf$temp, ylab="cross-correlation", main="PARmax & Temp", lag.max = 60)
