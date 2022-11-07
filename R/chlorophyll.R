# look at chlorophyll-a data from MODIS
# data gather by Sam Starko Nov 2022
# script started 6 November 2022 by Matt Whalen


# packages
library(tidyverse)
library(lubridate)
library(hydroTSM)


# read the data
d <- read_csv("Data/environmental_data/chlorophyll/g4.areaAvgTimeSeries.MODISA_L3m_CHL_8d_4km_2018_chlor_a.20120101-20220201.128W_51N_128W_51N.csv", skip =  8)
names(d)  <- c("time","mean_MODIS_raw")

# metadata header
read_csv("Data/environmental_data/chlorophyll/g4.areaAvgTimeSeries.MODISA_L3m_CHL_8d_4km_2018_chlor_a.20120101-20220201.128W_51N_128W_51N.csv")[1:6,]

# remove flagged values, deal with time
d <- d %>% 
  mutate( chl = ifelse(mean_MODIS_raw>0,mean_MODIS_raw,NA) ) %>% 
  mutate( tim = ymd_hms(time) ) %>%  mutate( year = year(time) )
# add sesasons
d <- d %>%
  mutate( season = time2season(d$time,out.fmt = "seasons") )
# remove 2022
d <- d %>% 
  filter(year < 2020)
# summaries
d %>% 
  filter(!is.na(chl)) %>% 
  group_by(season,year) %>% 
  summarize(nobs = length(chl) ) %>% arrange(year) %>% 
  pivot_wider(names_from = year, values_from = nobs )



# statistics
d$year.factor <- factor(d$year,ordered = T )
summary(aov( chl ~ year.factor, data = d))


# plot the time series
ggplot( d, aes(x = time, y = chl) ) + geom_point() + geom_line()

# plot summary data
# boxplot
ggplot(d, aes(x = as.factor(year), y = chl, group = year)) + geom_boxplot()
# histograms
ggplot(d, aes(x = chl)) + 
  facet_wrap( ~year ) +
  geom_histogram( stat = "bin") + 
  geom_density( bw = 0.1 ) +
  scale_x_log10()

ggplot(d, aes(x = chl)) + 
  facet_grid(season~year ) +
  geom_histogram( stat = "bin") + 
  geom_density( bw = 0.1 ) +
  scale_x_log10()


# isolate spring and summer, when we have the most data points
dss <- d %>% 
  filter( season %in% c("spring","summer") ) 
# boxplot
ggplot(dss, aes(x = as.factor(year), y = chl, group = year)) + 
  # facet_wrap(~season) + 
  geom_boxplot() +
  scale_y_log10() +
  ylab("MODIS chlorophyll-a concentration [units?]") +
  xlab("Year")

ggplot(dss, aes(x = chl)) + 
  facet_grid(season~year ) +
  geom_histogram( stat = "bin") + 
  geom_density( bw = 0.1 ) +
  scale_x_log10()


aov1 <- aov( log(chl) ~ year.factor, data = dss)
anova(aov1)
summary(aov1)
multcomp::glht(aov1,linfct = multcomp::mcp(year.factor = "Tukey"))

TukeyHSD(aov1)
