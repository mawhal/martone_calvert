###### Calvert Island thermal environment data

## The purpose of this script is to examine temporal trends in temperature data
# using iButton data from Foggy Cove to compare annual differences around the marine heatwaves

# by Matt Whalen
# created 3 November 2022

library(tidyverse)
library(lubridate)
# library(fpp)
# library(zoo)
# library(imputeTS)
library(cowplot)
library(hydroTSM)





## get data from iButtons all in one place
# for now, all from West Beach High Shore
# 2012+2013
ib13 <- read_csv( "Data/environmental_data/iButtons/ibuttons/12_WBB_High_BlueFalcon.csv", skip=14 )
ib13$`Date/Time` <- as.POSIXlt(ib13$`Date/Time`, format="%d/%m/%Y %H:%M")
ib13$`Date/Time` <- ymd_hms(ib13$`Date/Time`)
# 2015
ib15 <- read_csv( "Data/environmental_data/iButtons/ibuttons/2015 ibuttons/WBB_High_1A.csv", skip=14)
ib15$`Date/Time` <- as.POSIXlt(ib15$`Date/Time`, format="%m/%d/%y %I:%M:%S %p")
ib15$`Date/Time` <- ymd_hms(ib15$`Date/Time`)
# 2016
ib16 <- read_csv( "Data/environmental_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_High_N_BlackPVCCapW.csv",
                  skip=14 )
ib16$`Date/Time` <- as.POSIXlt(ib16$`Date/Time`, format="%m/%d/%y %I:%M:%S %p")
ib16$`Date/Time` <- ymd_hms(ib16$`Date/Time`)

# merge all of these
ib <- full_join( full_join( ib13,ib15),ib16 )
# write to disk
write_csv( ib, "R/output/FoggyCove_HIGH_iButton.csv")

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

#  add day rank - 
dateseq <- data.frame( date=seq.Date(as.Date("2010-01-01"),as.Date("2020-01-01"), by="day") )
dateseq <- dateseq %>% mutate( rank=1:nrow(dateseq))
ib <- left_join( ib, dateseq )







# check date ranges for each year
ib <- ib %>% mutate( year = year(lube.date))
ib %>% 
  group_by(year) %>% 
  summarize( range(lube.date) )



##### SEASONS
ib$season <- time2season(ib$lube.date,                # Convert dates to seasons
                         out.fmt = "seasons")
# fix spelling of autumn
ib$season[ ib$season == "autumm" ] <- "autumn"

ib %>% 
  group_by(year, season) %>% 
  summarize( min.date = min(lube.date), max.date = max(lube.date) )




# Plotting the time series
ib$season <- factor( ib$season, levels = c("winter","spring","summer","autumn") )
season.colors <- c( rgb(0,89,159, maxColorValue = 255), rgb(0,160,70, maxColorValue = 255), 
                    rgb(214,180,91, maxColorValue = 255), rgb(218,38,41, maxColorValue = 255))
ggplot( ib, aes( x = lube.date, y = value)) + 
  geom_point( aes(col=season), size = 1, alpha = 0.5 ) +
  ylab("Temperature (°C)") + xlab("Time") +
  scale_color_manual( values = c(season.colors) ) +
  theme_bw()  + theme( legend.position = "none" )
ggsave( "R/Figs/ibutton_FoggyCove_time_series.svg", height = 3, width =4 )


#### WINTERS
### can isolate full winter data for 2012/2013, 2014/2015, and 2015/2016

#### SUMMERS
# can isolate summer data for 3 July to end of summer for 2012, 
# 13 July to end summer 2014, full 2015

### SPRINGS
# can isolate start to May 21 for 2013
# full 2015, start to May 16 for 2016

### FALLS
# can isolate full fall for 2012
# full 2014
# full 2015




# remove dates that do not overlap
seasonal <- ib

## keep all winters
## keep all falls

## summers
seasonal <- seasonal %>% filter( lube.date > ymd("2012-07-13") ) %>% 
  filter( lube.date < ymd("2014-06-01") | lube.date > ymd("2014-07-13")) %>% 
  filter( lube.date < ymd("2015-06-01") | lube.date > ymd("2015-07-13"))

## springs
seasonal <- seasonal %>% filter( lube.date < ymd("2013-05-16") | lube.date > ymd("2013-05-31") )
seasonal <- seasonal %>% filter( lube.date < ymd("2015-05-16") | lube.date >= ymd("2015-06-01") ) 
seasonal <- seasonal %>% filter( lube.date < ymd("2016-05-16") | lube.date > ymd("2016-05-31") )




# group seasons by year and season, noting that winters start in one year and end in the next
seasonal <- seasonal %>% 
  mutate(season = factor(season, levels = c("spring", "summer","autumn","winter")) ) %>% 
  unite( year.season, year, season, sep = "-", remove = F)

seasonal$year.season[seasonal$year == 2012 & seasonal$season == "winter"] <- "2013-winter"
seasonal$year.season[seasonal$year == 2014 & seasonal$lube.date >= ymd("2014-12-01") ] <- "2015-winter"
seasonal$year.season[seasonal$year == 2015 & seasonal$lube.date >= ymd("2015-12-01") ] <- "2016-winter"

seasonal %>% 
  group_by(year.season) %>% 
  summarize( min.date = min(lube.date), max.date = max(lube.date) )





# plotting by season
ggplot(seasonal, aes(x = year, y = value, col= season, group = year.season )) +
  facet_wrap(~season, ncol = 4) +
  geom_boxplot( ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 



# observations above or below a threshold

# above 10, 20, 30, and 35 degrees C
over10 <- seasonal %>% 
  filter( value > 10 ) %>% 
  group_by( year.season, season ) %>% 
  summarize( over10 = n() ) %>% 
  arrange( season )
over20 <- seasonal %>% 
  filter( value > 20 ) %>% 
  group_by( year.season, season ) %>% 
  summarize( over20 = n() ) %>% 
  arrange( season )
over30 <- seasonal %>% 
  filter( value > 30 ) %>% 
  group_by( year.season, season ) %>% 
  summarize( over30 = n() ) %>% 
  arrange( season )
over35 <- seasonal %>% 
  filter( value > 35 ) %>% 
  group_by( year.season, season ) %>% 
  summarize( over35 = n() ) %>% 
  arrange( season )
under0 <- seasonal %>% 
  filter( value < 0 ) %>% 
  group_by( year.season, season ) %>% 
  summarize( under0 = n() ) %>% 
  arrange( season )
under5 <- seasonal %>% 
  filter( value < 5 ) %>% 
  group_by( year.season, season ) %>% 
  summarize( under5 = n() ) %>% 
  arrange( season )
under10 <- seasonal %>% 
  filter( value < 10 ) %>% 
  group_by( year.season, season ) %>% 
  summarize( under10 = n() ) %>% 
  arrange( season )

thresholds <- full_join( full_join( full_join( over10, over20), over30), over35 )
thresholds <- full_join(full_join(full_join(full_join( full_join( full_join( over10, over20), over30), over35), under0), under5), under10 )
# thresholds <- full_join( full_join( over20, over30), over35 )
thresholds[is.na(thresholds)] <- 0

# pivot longer
thresholds.long <- pivot_longer(thresholds, over10:under10)
# thresholds.long$year.season <- factor(thresholds.long$year.season, 
#                                       levels = c("2012-summer","2012-autumn"))
# separate (start of season) year from year.season
thresholds.long <- thresholds.long %>% 
  separate( year.season, c("year", "season"), sep = "-", remove = F) %>% 
  mutate( season = factor(season, levels = c( "summer", "autumn", "winter", "spring")))

# note exactly correct, but can covert to hours if we assume that each reading represents a certain the value over all four hours
thresholds.long$hours <- thresholds.long$value*4

# re-order the thresholds
thresholds.long <- thresholds.long %>% 
  mutate( threshold = factor(name, levels = c("under0", "under5", "under10",
                                         "over10", "over20", "over30", "over35"),
                        labels =  c("below\n0°C", "below\n5°C", "below\n10°C",
                            "above\n10°C", "above\n20°C", "above\n30°C", "above\n35°C") ) )

# label based on heatwave status
thresholds.long$status <- "during"
thresholds.long$status[ thresholds.long$year %in% c(2012,2013)  ] <- "before"
thresholds.long$status <- factor( thresholds.long$status, levels = c("before","during"))
thresholds.long$`heatwave\nstatus` <- thresholds.long$status
cols.two <- c( "#A2CDE2",  "#BF7C61" )
# plot summaries
ggplot( thresholds.long, aes(x = year, y = value, col = season, shape = `heatwave\nstatus` ) ) +  # y = hours
  facet_grid(threshold~season, scales = "free_y") + 
  geom_point( size = 2.5) + geom_line( aes(group = 1) ) +
  scale_y_continuous( breaks= scales::pretty_breaks()) + 
  scale_colour_manual( values = season.colors[c(3,4,1,2)] ) +
  scale_shape_manual( values = c(21, 19) ) +
  ylab("Number of measurements (per 4 hours)") + xlab("Year") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

ggsave( "R/Figs/iButton_FogyCove_threshold_season.svg", height = 6, width = 6 )

