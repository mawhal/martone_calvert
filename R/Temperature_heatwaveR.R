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
pine <- read_csv( "Data/environmetal_data/Lighthouse Data/2020 update/DATA_-_Active_Sites/Pine_Island/Pine_Island_-_Daily_Sea_Surface_Temperature_and_Salinity_1937-2019.csv",
skip=1 )
names(pine) <- c( 'date','sal','temp','latitude', 'longitude' )
mcin <- read_table( "data/environmetal_data/Lighthouse Data/through May 2019_Peter Chandler/McInnesDailySalTemp.txt",
                            skip=3 )

# air temperature - Addenbrooke air temperature
adden <- read_csv( "data/environmetal_data/Addenbroke Air Temperature/EC/1060080.ascii", skip=1 )  # data from https://data.pacificclimate.org/portal/pcds/map/

# clean up the data
# renames columns, make date columns, replace 999.9 with NA
d <- pine %>%
  filter( !is.na(date) ) %>% 
  mutate( date = mdy(date) ) %>% 
  # select( year=Year, month=Month, day=Day, temp=`Temperature(C)`, sal=`Salinity(psu)` ) %>%
  # mutate( date=ymd(paste(year,month,day))) %>%
  # unite( date, year,month,day, sep="-" ) %>%
  # mutate( date= ymd(date) ) %>%
  mutate( temp=replace(temp, temp==99.9, NA)) %>%
  mutate( sal=replace(sal, sal==99.9, NA))
#
d2 <- mcin %>%
  select( year=Year, month=Month, day=Day, temp=`Temperature(C)`, sal=`Salinity(psu)` ) %>%
  # mutate( date=ymd(paste(year,month,day))) %>%
  unite( date, year,month,day, sep="-" ) %>%
  mutate( date= ymd(date) ) %>%
  mutate( temp=replace(temp, temp==999.9, NA)) %>%
  mutate( sal=replace(sal, sal==999.9, NA)) 
#
d3_all <- adden %>%
  select( date=time, temp_min=MIN_TEMP, temp_max=MAX_TEMP,precip=ONE_DAY_PRECIPITATION ) %>%
  mutate( date= ymd(date), precip=log10(precip+1) ) 
#

# select tempearture only for now
d <- d %>%
  select( t=date, temp ) 
d2 <- d2 %>%
  select( t=date, temp )
d3 <- d3_all %>%
  select( t=date, temp=temp_max )
d4 <- d3_all %>%
  select( t=date, temp=temp_min )
d5 <- d3_all %>%
  select( t=date, temp=precip )


plot(d, type='l', ylim=c(0,20))
lines(d2, type='l',col="blue")


## heatwaveR
start.date <- "1937-01-01"
end.date   <- "2019-10-31"
# Detect the events in a time series
ts <- ts2clm(d, climatologyPeriod = c(start.date, end.date))
ts10 <- ts2clm(d, climatologyPeriod = c(start.date, end.date), pctile=10)
mhw <- detect_event(ts)
mcw <- detect_event(ts10,coldSpells = TRUE)
ts2 <- ts2clm(d2, climatologyPeriod = c("1954-01-01", "2019-05-01"))
ts210 <- ts2clm(d2, climatologyPeriod = c("1954-01-01", "2019-05-01"), pctile=10)
mhw2 <- detect_event(ts2)
mcw2 <- detect_event(ts210,coldSpells = TRUE)
ts3 <- ts2clm(d3, climatologyPeriod = c("1978-01-01", "2019-05-01"))
ts310 <- ts2clm(d3, climatologyPeriod = c("1978-01-01", "2019-05-01"), pctile=10)
mhw3 <- detect_event(ts3, minDuration = 1)
mcw3 <- detect_event(ts310,coldSpells = TRUE, minDuration = 1)
ts4 <- ts2clm(d4, climatologyPeriod = c("1978-01-01", "2019-05-01"))
ts410 <- ts2clm(d4, climatologyPeriod = c("1978-01-01", "2019-05-01"), pctile=10)
mhw4 <- detect_event(ts4)
mcw4 <- detect_event(ts410,coldSpells = TRUE)
ts5 <- ts2clm(d5, climatologyPeriod = c("1978-01-01", "2019-05-01"))
ts510 <- ts2clm(d5, climatologyPeriod = c("1978-01-01", "2019-05-01"), pctile=10)
mhw5 <- detect_event(ts5)
mcw5 <- detect_event(ts510,coldSpells = TRUE)

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

mhw3$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  head(6)

event_line(mhw, start_date = "2012-01-01", end_date = "2019-05-01") #, min_duration = 30, spread = 180, start_date = "1937-01-01", end_date = "2019-05-01")


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
event_line(mcw3, spread = 180, metric = "intensity_max", 
           start_date = "2011-01-01", end_date = "2019-05-01") # most of 2015 that we have a record for
event_line(mhw3, spread = 180, metric = "intensity_cumulative", 
           start_date = "2011-01-01", end_date = "2019-05-01") # most of 2015 that we have a record for


lolli_plot(mhw, metric = "intensity_max")
lolli_plot(mcw, metric = "intensity_max")
lolli_plot(mhw2, metric = "intensity_max")
lolli_plot(mcw2, metric = "intensity_max")
lolli_plot(mhw3, metric = "intensity_max")
lolli_plot(mcw3, metric = "intensity_max")

lolli_plot(mhw, metric = "intensity_cumulative")
lolli_plot(mhw2, metric = "intensity_cumulative")
lolli_plot(mhw3, metric = "intensity_cumulative")
lolli_plot(mhw, metric = "duration")
lolli_plot(mhw2, metric = "duration")
lolli_plot(mhw3, metric = "duration")
lolli_plot(mhw, metric = "intensity_mean")
lolli_plot(mhw2, metric = "intensity_mean")
lolli_plot(mhw3, metric = "intensity_mean")


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
mhw3_clim <- mhw3$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Addenbrooke Island MAX")
mcw3_clim <- mcw3$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Addenbrooke Island MAX")
mhw4_clim <- mhw4$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Addenbrooke Island MIN")
mcw4_clim <- mcw4$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Addenbrooke Island MIN")
mhw5_clim <- mhw5$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Addenbrooke Island PRECIP")
mcw5_clim <- mcw5$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( site="Addenbrooke Island PRECIP")
mhw_clim_join <- full_join(full_join(full_join(full_join( mhw_clim, mhw2_clim ),mhw3_clim),mhw4_clim),mhw5_clim)#, by=c("site","doy", "t") )
mcw_clim_join <- full_join(full_join(full_join(full_join( mcw_clim, mcw2_clim ),mcw3_clim),mcw4_clim),mcw5_clim)#, by=c("site","doy", "t") )

ggplot(mhw_clim_join, aes(x = t, y = temp, y2 = thresh, col=site, fill=site)) +
  geom_flame(n=5,n_gap = 2, alpha=0.5) 
ggplot(mcw_clim_join, aes(x = t, y2 = temp, y = thresh, col=site, fill=site)) +
  geom_flame(n=5,n_gap = 2, alpha=0.5) 

ggplot(mhw_clim_join, aes(x = t), group=site) +
  facet_wrap(~site, ncol=1, scales = "free_y") +
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
        dpi=300, width=11, height=8 )


##Separate figures
##Pine Island

ggplot(mhw_clim, aes(x = t), group=site) +
  facet_wrap(~site, ncol=1, scales = "free_y") +
  geom_flame(aes(y = temp, y2 = thresh), 
             n=5,n_gap=2,
             show.legend = FALSE, col="gray25", fill='#B51D2C', size=0.5 ) +
  geom_flame(data=mcw_clim,aes(y2 = temp, y = thresh), 
             n=5,n_gap=2,
             show.legend = FALSE, col="gray25", fill='#296BA8', size=0.5 ) +
  # geom_flame(data = mhw_top, aes(y = temp, y2 = thresh, fill = "top"),  show.legend = T) +
  geom_line(aes(y = temp, colour = "temp"), size=0.5, alpha=0.5) +
  geom_line(aes(y = thresh, colour = "thresh"), size = 0.5) +
  geom_line(data=mcw_clim,aes(y = thresh, colour = "thresh2"), size = 0.5) +
  geom_line(aes(y = seas, colour = "seas"), size = 0.5) +
  scale_colour_manual(name = "Line Colour",
                      values = c("temp" = "gray25",
                                 "thresh" =  "#B51D2C",
                                 "thresh2" =  "#296BA8",
                                 "seas" = "black")) +
  scale_x_date(date_labels = "%b %Y") + theme_classic()+
  guides(colour = guide_legend(override.aes = list(fill = NA))) +
  labs(y = expression(paste("Temperature (", degree, "C)")), x = NULL) 
  theme_minimal()
  #Save as 4 x 8

##McInnes Island
  
  ggplot(mhw2_clim, aes(x = t), group=site) +
    facet_wrap(~site, ncol=1, scales = "free_y") +
    geom_flame(aes(y = temp, y2 = thresh), 
               n=5,n_gap=2,
               show.legend = F, col="gray25", fill='#B51D2C', size=0.5 ) +
    geom_flame(data=mcw2_clim,aes(y2 = temp, y = thresh), 
               n=5,n_gap=2,
               show.legend = F, col="gray25", fill="#296BA8", size=0.5 ) +
    geom_line(aes(y = temp, colour = "temp"), size=0.5, alpha=0.5) +
    geom_line(aes(y = thresh, colour = "thresh"), size = 0.5) +
    geom_line(data=mcw2_clim,aes(y = thresh, colour = "thresh2"), size = 0.5) +
    geom_line(aes(y = seas, colour = "seas"), size = 0.5) +
    scale_colour_manual(name = "Line Colour",
                        values = c("temp" = "gray25",
                                   "thresh" =  "#B51D2C",
                                   "thresh2" =  "#296BA8",
                                   "seas" = "black")) +
    scale_x_date(date_labels = "%b %Y") + theme_classic()+
    guides(colour = guide_legend(override.aes = list(fill = NA))) +
    labs(y = expression(paste("Temperature (", degree, "C)")), x = NULL) 
  theme_minimal()
  
###Pine Island
year<-c("2011","2012","2013","2014","2015","2016","2017","2018","2019") %>% as.numeric()
heatwave_days<-c(0,0,26,92,249,191,28,48,24)
pine<-tibble(year, heatwave_days)


ggplot(pine, aes(x=year, y=heatwave_days))+
  geom_point(size=4)+
  geom_line(aes(y=heatwave_days))+
  theme_classic()+
  scale_x_continuous(name="Year", limits=c(2011, 2019), breaks=seq(2011, 2019, 1))+
  scale_y_continuous(name="Heatwave days")+
  ggtitle("Pine Island")



###McInnes Island
year<-c("2011","2012","2013","2014","2015","2016","2017","2018","2019") %>% as.numeric()
heatwave_days<-c(0,0,0,9,124,154,37,37,6)
mc<-tibble(year, heatwave_days)


ggplot(mc, aes(x=year, y=heatwave_days))+
  geom_point(size=4)+
  geom_line(aes(y=heatwave_days))+
  theme_classic()+
  scale_x_continuous(name="Year", limits=c(2011, 2019), breaks=seq(2011, 2019, 1))+
  scale_y_continuous(name="Heatwave days")+
  ggtitle("McInnes Island")


###Addenbroke
year<-c("2011","2012","2013","2014","2015","2016","2017","2018","2019") %>% as.numeric()
heatwave_days<-c(11, 18, 24, 61, 93,76,32,48,74)
ad<-tibble(year, heatwave_days)


ggplot(ad, aes(x=year, y=heatwave_days))+
  geom_point(size=4)+
  geom_line(aes(y=heatwave_days))+
  theme_classic()+
  scale_x_continuous(name="Year", limits=c(2011, 2019), breaks=seq(2011, 2019, 1))+
  scale_y_continuous(name="Heatwave days")+
  ggtitle("Addenbroke Island (Air temperature)")


