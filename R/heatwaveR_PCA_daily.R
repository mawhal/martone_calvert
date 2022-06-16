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
library( imputeTS )


## try with data we already have 
pine <- read_csv( "Data/environmental_data/Lighthouse Data/2021_update/DATA_-_Active_Sites/Pine_Island/Pine_Island_-_Daily_Sea_Surface_Temperature_and_Salinity_1937-2021.csv",
                  skip=1 )
names(pine)[1:5] <- c( 'date','sal','temp','latitude', 'longitude' )

# read daily PCA scores from Lightstation_raw_to_impute
pca.meta <-  read_csv("Data/R code for Data Prep/output from R/Lightstation_raw_PCA_impute_daily.csv")

# read study sites info
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )
survey.dates <- am %>% 
  filter( Year != "2011" ) %>% filter( Site != "Meay Channel") %>% 
  select( Date, Year) %>% distinct()


# clean up the data
# renames columns, make date columns, replace 999.9 with NA
d <- pine %>%
  filter( !is.na(date) ) %>% 
  mutate( date = ymd(date) ) %>% 
  mutate( temp=replace(temp, temp==999.9, NA)) %>%
  mutate( sal=replace(sal, sal==999.9, NA))
#
# select tempearture only for now
d <- d %>%
  select( t=date, temp ) 

## interpolate data from Pine Island
# remove the first few rows
d <- d[-c(1:13),]
dinterp <- na_interpolation(d)
dinterp <- dinterp %>% 
  filter( t < "2019-10-31")
## heatwaveR
start.date <- "1937-01-14"
end.date   <- "2019-10-30"
maxGapchoice = 2
ts <- ts2clm(dinterp, climatologyPeriod = c(start.date, end.date))
ts10 <- ts2clm(dinterp, climatologyPeriod = c(start.date, end.date), pctile=10)
mhw <- detect_event(ts, maxGap = maxGapchoice)
mcw <- detect_event(ts10,coldSpells = TRUE, maxGap = maxGapchoice)


### heatwave analysis on Pine Island

### "heatwave" analysis
# using PC scores on a daily timescale, includes salinity as well as temperature
# select PCA1 only for now
d <- pca.meta %>%
  select( t=date, temp=pca1 ) 

plot(d, type='l')


## heatwaveR
start.date <- "1978-01-01"
end.date   <- "2019-10-31"
maxGapchoice = 2

# Detect the events in a time series
ts <- ts2clm(d, climatologyPeriod = c(start.date, end.date))
ts10 <- ts2clm(d, climatologyPeriod = c(start.date, end.date), pctile=10)
mhw2 <- detect_event(ts, maxGap = maxGapchoice)
mcw2 <- detect_event(ts10,coldSpells = TRUE, maxGap = maxGapchoice)

# View just a few metrics
mhw2$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  head(6)





# look at the time period of the survey - Summer 2011 to summer 2019
# It is necessary to give geom_flame() at least one row on either side of 
# the event in order to calculate the polygon corners smoothly
mhw_clim <- mhw$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( source="Pine Island" )
mcw_clim <- mcw$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( source="Pine Island" )
mhw2_clim <- mhw2$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( source="PCA")
mcw2_clim <- mcw2$climatology %>% 
  filter( t >= "2011-05-01") %>%
  mutate( source="PCA")

mhw_clim_join <- full_join( mhw_clim, mhw2_clim )
mcw_clim_join <- full_join( mcw_clim, mcw2_clim )

ggplot(mhw_clim_join, aes(x = t, y = temp, y2 = thresh, col=source, fill=source)) +
  geom_flame(n=5,n_gap = 2, alpha=0.5) 
ggplot(mcw_clim_join, aes(x = t, y2 = temp, y = thresh, col=source, fill=source)) +
  geom_flame(n=5,n_gap = 2, alpha=0.5) 

ggplot(mhw_clim_join, aes(x = t), group=source) +
  facet_wrap(~source, ncol=1, scales = "free_y") +
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


ggsave( "R/Figs/heatwaveR_pine_pca_2011.svg",
        dpi=300, width=11, height=8 )



# just show PCA for survey times
mhw2_clim <- mhw2$climatology %>% 
  filter( t >= "2011-07-01") %>%
  mutate( source="PCA")
mcw2_clim <- mcw2$climatology %>% 
  filter( t >= "2011-07-01") %>%
  mutate( source="PCA")
ggplot(mhw2_clim, aes(x = t)) +
  geom_vline(xintercept = survey.dates$Date, lty=2, size=0.5) +
  geom_flame(aes(y = temp, y2 = thresh), 
             n=5,n_gap=2,
             show.legend = T, col="gray25", fill='red', size=0.5 ) +
  geom_flame(data=mcw2_clim,aes(y2 = temp, y = thresh), 
             n=5,n_gap=2,
             show.legend = T, col="gray25", fill='blue', size=0.5 ) +
  # geom_flame(data = mhw_top, aes(y = temp, y2 = thresh, fill = "top"),  show.legend = T) +
  geom_line(aes(y = temp, colour = "PCA1"), size=0.125) +
  geom_line(aes(y = thresh, colour = "90th percentile"), size = 0.5) +
  geom_line(data=mcw2_clim,aes(y = thresh, colour = "10th percentile"), size = 0.5) +
  geom_line(aes(y = seas, colour = "climatology"), size = 0.5) +
  scale_colour_manual(name = "Line Colour",
                      values = c("PCA1" = "black",
                                 "90th percentile" =  "red",
                                 "10th percentile" =  "blue",
                                 "climatology" = "gray25")) +
  # # scale_fill_manual(name = "Event Colour", 
  #                   values = c("all" = "salmon", 
  #                              "top" = "red")) +
  scale_x_date(date_labels = "%Y") +
  guides(colour = guide_legend(override.aes = list(fill = NA))) +
  labs(y = "PCA1", x = NULL) +
  theme_minimal() + theme(legend.justification=c(1,0), legend.position=c(1,0),legend.direction="horizontal",legend.title=element_blank(),
                          legend.background = element_rect(fill="white"),
                          legend.box.background = element_rect(colour = "black"),
                          panel.border = element_rect(colour = "black", fill=NA),
                          text = element_text(size = 20))

ggsave("R/Figs/lighthouse_heatwaveR_PCA.svg", width=10, height = 2)
  
  


# Make a table to show heatwave summary for Pine Island and PCA scores
mhw_summary <- mhw$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  mutate(source = "Pine Island")
mhw2_summary <- mhw2$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  mutate(source = "PCA")
mhw_summary_all <- bind_rows(mhw_summary, mhw2_summary)

# get yearly summary
mhw_summary_year <- mhw_summary_all %>% 
  mutate( year = year(date_start) ) %>% 
  group_by( year, source ) %>% 
  summarize( duration = sum(duration), intensity_max = max(intensity_max)) %>% 
  arrange(-duration) 
ggplot( mhw_summary_year, aes(x=year, y=duration)) +
  facet_wrap(~source)+
  geom_segment( aes(yend=duration,y=0,xend=year),col='firebrick') + geom_point( aes(size=intensity_max), alpha = 0.7, col='firebrick')
ggplot( filter(mhw_summary_year,year>=2011), aes(x=year, y=duration)) +
  facet_wrap(~source)+
  geom_segment( aes(yend=duration,y=0,xend=year),col='firebrick') + geom_point( aes(size=intensity_max), alpha = 0.7, col='firebrick')
  

lolli_plot(mhw, metric = "intensity_cumulative")
lolli_plot(mhw, metric = "intensity_max")
# write to disk
write_csv( mhw_summary_year, "R/output/heatwaveR_summary_year_pine_pca.csv" )



# get yearly summary for everything leading up to each survey
# read survey dates
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )
survey.dates <- am %>% 
  filter( Year != "2011" ) %>% filter( Site != "Meay Channel") %>% 
  select( Date, Year) %>% distinct() %>% 
  group_by(Year) %>% 
  summarize( Date = min(Date) ) %>% 
  mutate( Date_prev_year = Date - years(1) ) %>% 
  mutate( Date_prev_year5 = Date - years(5) ) 

mhw_summary_all_survey <- mhw_summary_all %>% 
  filter( date_start > "2011-01-01" ) 







# use climatology to get cumulative intensity, duration, and max intensity for PCA data
survey.dates$duration <- NA
survey.dates$intensity_cummulative <- NA
survey.dates$intensity_max <- NA
survey.dates$duration5 <- NA
survey.dates$intensity_cummulative5 <- NA
survey.dates$intensity_max5 <- NA

for(i in 1:nrow(survey.dates)){
  tmp <- mhw2$climatology[ mhw2$climatology$t %in% seq(survey.dates$Date_prev_year[i],survey.dates$Date[i],by="days"), ] %>% 
    filter(durationCriterion == T) %>% 
    mutate(anomaly = temp-seas)
  survey.dates$duration[i] <- nrow(tmp) 
  survey.dates$intensity_cummulative[i] <- sum(tmp$anomaly) 
  survey.dates$intensity_max[i] <- ifelse( nrow(tmp)>0, max(tmp$anomaly),0 )
  tmp <- mhw2$climatology[ mhw2$climatology$t %in% seq(survey.dates$Date_prev_year5[i],survey.dates$Date[i],by="days"), ] %>% 
    filter(durationCriterion == T) %>% 
    mutate(anomaly = temp-seas)
  survey.dates$duration5[i] <- nrow(tmp) 
  survey.dates$intensity_cummulative5[i] <- sum(tmp$anomaly) 
  survey.dates$intensity_max5[i] <- ifelse( nrow(tmp)>0, max(tmp$anomaly),0 )
}
psych::pairs.panels( select(survey.dates,Date,duration,intensity_max,duration5) )

write_csv( survey.dates, "R/output/heatwaveR_duration_surveyyear_pca.csv" )





# a few simple plots
# svg("R/Figs/heatwave_pca_duration_cumm_time.svg", width=2.5, height=2)
par(mar = c(3,4,1,1)+0.1, las=1, pty='s')
plot( duration5~Year, survey.dates, type='o', col="firebrick", pch=19, xaxt="n",
      ylab = "Heatwave days", xlab = "" )
lines( duration~Year, survey.dates, col='black' )
points( duration~Year, survey.dates,col='black' )
# axis(1, at = seq(2012,2018,by=2))
# axis(1, at = c(2014,2018))
axis(1, at = 2012:2019, labels=FALSE)
text( x=seq(2012,2018,by=2), 
      y=par("usr")[3] - 75, xpd = NA,
      labels = seq(2012,2018,by=2), 
      srt = 35, adj = 0.965 ) 
# dev.off()



survey.dates %>% 
  select( Date, duration, duration5 ) %>% 
  pivot_longer(cols = duration:duration5) %>% 
  mutate(name=factor(name,levels=c("duration5","duration"), labels=c("previous 5 years","previous year"))) %>% 
  ggplot( aes( x = Date, y = value, col=name, shape=name, fill=name)) + geom_line() + geom_point(size=8) +
  ylab("Heatwave\ndays") +
  scale_color_manual(values=c("firebrick","black")) +
  scale_fill_manual(values=c("firebrick","white")) +
  scale_shape_manual(values=c(20,1)) +
  coord_cartesian( xlim = c(ymd("2011-10-31"),ymd("2019-10-31"))) +
  theme_minimal() + theme(legend.justification=c(0,1), legend.position=c(0,1),legend.direction="horizontal",legend.title=element_blank(),
                          legend.background = element_rect(fill="white"),
                          legend.box.background = element_rect(colour = "black"),
                          panel.border = element_rect(colour = "black", fill=NA),
                          text = element_text(size = 20))
ggsave("R/Figs/lighthouse_heatwave_days.svg",width=10,height=1.75)




### ----------------------------------
# calculate climatology for each time series, plus pca, and make a single table to show heatwave stats
# what we want: June 1 to May 31 temperatures
