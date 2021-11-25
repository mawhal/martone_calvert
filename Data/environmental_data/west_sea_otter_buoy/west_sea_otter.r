# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen


# This script organizes and plots wave data from the West Sea Otter Buoy, BC


# load libraries
library( tidyverse )
library( lubridate )

# read data
d <- read_csv( "c46204.csv" )

# filter out flagged 
d <- d %>% 
  filter( Q_FLAG==1 )

# dates
d <- d %>% 
  mutate( DATETIME = mdy_hm(DATE),
          date = date(DATETIME),
          month = month((DATETIME)),
          year = year(DATETIME) )

# average by day
dmean <- d %>% 
  group_by( date ) %>% 
  summarize( VCAR=mean(VCAR), VCMX=mean(VCMX) ) 
dmean2 <- d %>% 
  group_by( year,month ) %>% 
  summarize( VCAR=mean(VCAR), VCMX=mean(VCMX) ) %>% 
  unite( month, year, "month", sep="-" ) %>% 
  mutate( date=ymd(month, truncated=1) )

# plot the raw wave heights
names(d)

plot( VCAR ~ DATE, d )

ggplot( dmean, aes( x=date, y=VCAR )) + geom_point() +
  geom_smooth(se=F)
ggplot( dmean2, aes( x=date, y=VCMX )) + geom_point() + geom_line() + 
  geom_smooth(se=F) +
  geom_vline( xintercept=ymd( c("2014-03-01","2017-07-01") ), col='red' )
ggplot( dmean2, aes( x=date, y=VCAR )) + geom_point() + geom_line() + 
  geom_smooth(se=F) +
  geom_vline( xintercept=ymd( c("2014-03-01","2017-07-01") ), col='red' )


# use heatwave package
library(heatwaveR)

dsimple <- dmean %>% select( t=date, y=VCAR )
ts <- ts2clm( dsimple, y=y, climatologyPeriod = date(range(dsimple$t)) )

mhw <- detect_event(ts, y=y )

# biggest
mhw$event %>% arrange(-intensity_mean,-duration) %>% 
  select( date_start,date_end,duration,
          intensity_mean,intensity_max )

# during the heatwave
mhw$event %>% filter(date_start>ymd("2014-01-01")) %>% 
  select( date_start,date_end,duration,
          intensity_mean,intensity_max, )

mhw$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  head(6)

event_line(mhw, y=y, metric = "intensity_mean")
event_line(mhw, y=y, metric = "intensity_max")
event_line(mhw, y=y, metric = "intensity_max", #ylab="wave hieght (m)",
           start_date = "2011-01-01", end_date = "2019-12-01")
event_line(mhw, y=y, metric = "intensity_max", #ylab="wave hieght (m)",
           start_date = "2014-01-01", end_date = "2016-12-01")


lolli_plot(mhw,  metric = "intensity_mean")


mhw_clim <- mhw$climatology %>% 
  filter( t >= "2011-05-01")

windows(12,3)
ggplot(mhw_clim, aes(x = t, y = y, y2 = thresh)) + geom_line( alpha=0.2) +
  geom_line( aes(y=zoo::rollmean(y, 30, na.pad=TRUE)), alpha=0.5, col='red' ) +
  geom_flame(n=0,n_gap = 0,col='black') + 
  ylab("Characteristic significant\nwave height (m)") + xlab("") +
  theme_minimal()
