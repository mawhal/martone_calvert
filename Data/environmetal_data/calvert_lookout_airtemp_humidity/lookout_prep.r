### Hakai Institute  weather station data
# data from the Lookout on Calvert Island
# Whalen

library(tidyverse)
library(lubridate)

# read data
draw <- read_csv( "Data/environmetal_data/calvert_lookout_airtemp_humidity/2020-04-21.1hourSamples.csv" )
draw <- draw[-c(1:3),]
names(draw)

# select columns
dsel <- draw %>% select( date_time = `Measurement time`, Year, Month, WaterYear, 
                 temp1q=TAirLookout1_Q_level,temp1mean=TAirLookout1_Avg,temp2q=TAirLookout2_Q_level,temp2mean=TAirLookout2_Avg,
                 rhq=RHLookoutQ_level,rhmean=RHLookoutAvg)

# dates, numeric
d <- dsel %>% mutate( date_time=ymd_hms(date_time), date=date(date_time) ) %>% 
  mutate( temp1mean=as.numeric(temp1mean),
          temp2mean=as.numeric(temp2mean),
          rhmean=as.numeric(rhmean) )

# daily data
dday <- d %>% 
  group_by( date_time, Year, Month, date ) %>% 
  mutate( temp = mean(c(temp1mean,temp2mean),na.rm=T)) %>% 
  group_by( Year, Month, date ) %>% 
  summarize( temp=mean(temp),rhmean=mean(rhmean) )

# use heatwaveR to get seasonal trend
library(heatwaveR)
# select tempearture only for now
d <- dday %>%
  ungroup() %>% 
  select( t=date, temp ) 
ts <- ts2clm(d, climatologyPeriod = range(d$t) )
ts10 <- ts2clm(d, climatologyPeriod = range(d$t), pctile=10)

a <- ggplot( ts, aes(x=t, y=temp, col=temp-seas) ) + 
  geom_point( alpha=1 ) +
  geom_line( aes(y=seas),col='gray50' ) +
  geom_line( aes(y=thresh),col='gray75' ) +
  geom_line( data=ts10, aes(y=thresh),col='gray75' ) +
  guides(col=guide_colorbar(title="seasonal\nanomaly")) +
  scale_color_gradient2(low = "#296BA8",mid="whitesmoke", high = "#B51D2C") + 
  ylab( expression(paste("Temperature (",degree,"C)")) ) + xlab("")
ggsave( "Data/environmetal_data/calvert_lookout_airtemp_humidity/lookout_daily.svg", width=6,height=2 )


# select tempearture only for now
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
