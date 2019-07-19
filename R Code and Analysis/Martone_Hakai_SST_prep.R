# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses SST data from Pine Island Lighthouse to characterize 
# surface temperature and temperature anomalies over the study period
# These data can then be used as covariates to help explain
# patterns in population and community dynamics over space and time

## load libraries
library( tidyverse )
library( lubridate )

## read data 
pine <- read_csv( "../Data/Environmental Data/Lighthouse Data/through May 2019_Peter Chandler/output from R/PineIsland_monthly_SST_anomaly.csv" )
# also read the meta data for the project so we can identify sampling dates
am <- read_csv( "../Data/R code/Output from R/Martone_Hakai_metadata.csv" )


## get relevant dates from the survey metadata
# keep it simple and average dates within each year
Dates <- am %>%
  group_by( Year ) %>%
  summarize( Date=mean(Date) ) %>%
  mutate( Date.month=format(Date,format="%Y-%m"))

Dates$month <- format(Dates$Date,format="%m")
Dates$start <- Dates$Date-335
day(Dates$start) <- 1
#
# add 2019 samples
m19 <- data.frame( Year=2019, Date=ymd("2019-06-02"), Date.month="2019-06",month="06",start=ymd("2018-06-01") )
Dates <- full_join( Dates, m19 )

# define the preceding 12 months of data prior to each sampling event, including the month of the sample
date.ranges <-  data.frame( sapply( Dates$start, function(z) seq( z,length=12,by="months" ), simplify = FALSE )  )
names(date.ranges) <- Dates$Date
date.ranges



# consider adding multiple windows for temperature anomalies (6 months, 12, 18, 24) to see which provides the best statistical fit to the data




## use the sample dates and date.ranges to extract and summarize data from pine island
# get single month temperature and anomalies
instant <- Dates$Date
day(instant) <- 1
# for 2019 only, before data fully available
instant[9] <- "2019-05-01"
pine <- pine %>% mutate( Date= ymd(paste(year,month,1)) )
inst.measure <- pine[ pine$Date %in% instant, ]
inst.measure$measurement <- "month"
inst.measure$sample.date <- as.character(Dates$Date)

# average all temperatures over the period defined in date.ranges
date.ranges
data.range  <- apply( date.ranges, 2, function(z) {
  return(pine[ pine$Date %in% as.Date(z), ]) })

# combine list elements and keep the list name as an identifying column
dr.bind <- bind_rows( data.range,.id = "name")

# summarize the data to get average Temp, and average anomaly
yearly <- dr.bind %>%
  group_by(name) %>%
  summarise( Temp.year=mean(temp,na.rm=T), Anomaly.year=mean(temp.anom) ) %>%
  mutate( sample.date=as.character(name) )

# join the monthly temperature data with the yearly average
combined <- full_join( inst.measure, yearly, but="sample.date" )


# quick plot to look at how temperature and anomolies are related
psych::pairs.panels( combined[,c(1,3,16,17)] )


## write the environmental data to disk
write_csv( combined, "output from r/PineIsland_summary.csv" )
