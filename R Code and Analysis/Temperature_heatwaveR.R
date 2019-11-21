## Martone Hakai Rocky Shore Seaweed Surveys
# 
# Script to classify hot and cold spells 
#
# started by Matt Whalen 7 Octoboer 2019
#
# 

# library( rerddap )
library( tidyverse )
library( lubridate )
library( heatwaveR )
# # The two packages we will need
# # NB: The packages only need to be installed from GitHub once
# devtools::install_github("tidyverse/tidyverse")
# devtools::install_github("ropensci/rerddap")


## try with data we already have 
pine <- read_table( "Data/Environmental Data/Lighthouse Data/through May 2019_Peter Chandler/PineDailySalTemp.txt",
                  skip=3 )
mcin <- read_table( "Data/Environmental Data/Lighthouse Data/through May 2019_Peter Chandler/McInnesDailySalTemp.txt",
                            skip=3 )


# clean up the data
# renames columns, make date columns, replace 999.9 with NA
d <- pine %>%
  select( year=Year, month=Month, day=Day, temp=`Temperature(C)`, sal=`Salinity(psu)` ) %>%
  # mutate( date=ymd(paste(year,month,day))) %>%
  unite( date, year,month,day, sep="-" ) %>%
  mutate( date= ymd(date) ) %>%
  mutate( temp=replace(temp, temp==999.9, NA)) %>%
  mutate( sal=replace(sal, sal==999.9, NA)) 
#
d2 <- mcin %>%
  select( year=Year, month=Month, day=Day, temp=`Temperature(C)`, sal=`Salinity(psu)` ) %>%
  # mutate( date=ymd(paste(year,month,day))) %>%
  unite( date, year,month,day, sep="-" ) %>%
  mutate( date= ymd(date) ) %>%
  mutate( temp=replace(temp, temp==999.9, NA)) %>%
  mutate( sal=replace(sal, sal==999.9, NA)) 
#

# select tempearture only for now
d <- d %>%
  select( t=date, temp ) 
d2 <- d2 %>%
  select( t=date, temp )


plot(d, type='l', ylim=c(0,20))
lines(d2, type='l',col="blue")


## heatwaveR
# Detect the events in a time series
ts <- ts2clm(d, climatologyPeriod = c("1954-01-01", "2019-05-01"))
ts10 <- ts2clm(d, climatologyPeriod = c("1954-01-01", "2019-05-01"), pctile=10)
mhw <- detect_event(ts)
mcw <- detect_event(ts10,coldSpells = TRUE)
ts2 <- ts2clm(d2, climatologyPeriod = c("1954-01-01", "2019-05-01"))
ts210 <- ts2clm(d2, climatologyPeriod = c("1954-01-01", "2019-05-01"), pctile=10)
mhw2 <- detect_event(ts2)
mcw2 <- detect_event(ts210,coldSpells = TRUE)

# View just a few metrics
mhw$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  head(6)

mhw2$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  head(6)


event_line(mhw, min_duration = 30, spread = 180, metric = "intensity_max", # big El Nino of 1997
           start_date = "1937-01-01", end_date = "2019-05-01")
event_line(mcw, spread = 180, metric = "intensity_max", # big El Nino of 1997
           start_date = "2011-01-01", end_date = "2019-05-01")
event_line(mhw, min_duration = 10, spread = 180, metric = "intensity_max", # acute event in 2013
           start_date = "2000-01-01", end_date = "2019-05-01") #strarted 2013-08-28, peaked 2013-09-01, 16 day duration
event_line(mhw, min_duration = 30, spread = 180, metric = "intensity_cumulative", 
           start_date = "2000-01-01", end_date = "2019-05-01") # most of 2015 that we have a record for
event_line(mhw2, min_duration = 30, spread = 180, metric = "intensity_cumulative", 
           start_date = "2000-01-01", end_date = "2019-05-01") # most of 2015 that we have a record for
event_line(mcw2, spread = 180, metric = "intensity_max", 
           start_date = "2011-01-01", end_date = "2019-05-01") # most of 2015 that we have a record for


lolli_plot(mhw, metric = "intensity_max")
lolli_plot(mcw, metric = "intensity_max")
lolli_plot(mhw2, metric = "intensity_max")
lolli_plot(mcw2, metric = "intensity_max")

lolli_plot(mhw, metric = "intensity_cumulative")
lolli_plot(mhw2, metric = "intensity_cumulative")
lolli_plot(mhw, metric = "duration")
lolli_plot(mhw2, metric = "duration")
lolli_plot(mhw, metric = "intensity_mean")
lolli_plot(mhw2, metric = "intensity_mean")


# look at the time period of the survey - Summer 2011 to summer 2019
# It is necessary to give geom_flame() at least one row on either side of 
# the event in order to calculate the polygon corners smoothly
mhw_clim <- mhw$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Pine Island" )
mcw_clim <- mcw$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Pine Island" )
mhw2_clim <- mhw2$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="McInnis Island")
mcw2_clim <- mcw2$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="McInnis Island")
mhw_clim_join <- full_join( mhw_clim, mhw2_clim )#, by=c("site","doy", "t") )
mcw_clim_join <- full_join( mcw_clim, mcw2_clim )#, by=c("site","doy", "t") )

ggplot(mhw_clim_join, aes(x = t, y = temp, y2 = thresh, col=site, fill=site)) +
  geom_flame(n=5,n_gap = 2, alpha=0.5) 
ggplot(mcw_clim_join, aes(x = t, y2 = temp, y = thresh, col=site, fill=site)) +
  geom_flame(n=5,n_gap = 2, alpha=0.5) 

ggplot(mhw_clim_join, aes(x = t)) +
  facet_wrap(~site, ncol=1, scales = "fixed") +
  geom_flame(aes(y = temp, y2 = thresh), 
             n=5,n_gap=2,
             show.legend = T, col="gray25", fill='red', size=0.5 ) +
  geom_flame(data=mcw_clim_join,aes(y2 = temp, y = thresh), 
             n=5,n_gap=2,
             show.legend = T, col="gray25", fill='blue', size=0.5 ) +
  # geom_flame(data = mhw_top, aes(y = temp, y2 = thresh, fill = "top"),  show.legend = T) +
  geom_line(aes(y = temp, colour = "temp"), size=0.5, alpha=0.5) +
  geom_line(aes(y = thresh, colour = "thresh"), size = 0.5) +
  geom_line(data=mcw_clim_join,aes(y = thresh, colour = "thresh2"), size = 0.5) +
  geom_line(aes(y = seas, colour = "seas"), size = 0.5) +
  scale_colour_manual(name = "Line Colour",
                      values = c("temp" = "gray25",
                                 "thresh" =  "red",
                                 "thresh2" =  "blue",
                                 "seas" = "black")) +
  # # scale_fill_manual(name = "Event Colour", 
  #                   values = c("all" = "salmon", 
  #                              "top" = "red")) +
  scale_x_date(date_labels = "%b %Y") +
  guides(colour = guide_legend(override.aes = list(fill = NA))) +
  labs(y = expression(paste("Temperature (", degree, "C)")), x = NULL) +
  theme_minimal()


ggsave( "R Code and Analysis/Figs/heatwaveR_lighthouse_2011.png",
        dpi=300, width=12, height=6 )
