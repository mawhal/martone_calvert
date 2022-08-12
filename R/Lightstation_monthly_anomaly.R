## investigate raw data from lightstations instead of monthly averages or anomalies
library(tidyverse)
library(lubridate)
# d <- read_csv("Data/R code for Data Prep/Output from R/Lightstation_raw.csv") # see Lightstation_spectral_anomaly.R
d <- read_csv("Data/R code for Data Prep/Output from R/Lightstation_monthly_anomaly.csv") # see Lightstation_spectral_anomaly.R
# read study sites info
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )
survey.dates <- am %>% 
  filter( Year != "2011" ) %>% filter( Site != "Meay Channel") %>% 
  select( Date, Year) %>% distinct() %>% 
  summarize( date = ymd(range(Date)) ) 


# add proper date
d <- d %>% 
  mutate( date = ym(paste(year,month,sep="-")) )
# revise site names
d$Site <- factor( d$site, levels = c("mcinnes","pine","addenbroke"), labels = c("McInnes - sea surface temperature",
                                                                                "Pine - sea surface temperature",
                                                                                "Addenbroke - air temperature"))

# # take one piece of data at a time and then merge them
# dtemp <- d %>% select(date,temp,site) %>% group_by(date) %>% 
#   pivot_wider( names_from = site, values_from = temp, names_prefix = "temp_")
# dsal <- d %>% select(date,sal,site) %>% group_by(date) %>% 
#   pivot_wider( names_from = site, values_from = sal, names_prefix = "sal_")
# 
# # merge and clean
# dwide <- left_join(dtemp,dsal)# %>% select(-sal_addenbroke)



drecent <- d %>% filter( date >= "2011-07-01" & date <= "2019-10-31") %>% ungroup()


# check the data for outliers, etc.
summary(drecent)



# show monthly temperature anomalies
d <- d %>% mutate( sign = as.character(as.numeric(temp.anom>0)) )
a <- ggplot( data = d, aes( x=date, y=temp.anom,col=sign)) +
  facet_wrap(~site, ncol=1) +
  geom_vline(xintercept = survey.dates$date, lty=2, size=0.1) +
  geom_hline(yintercept = 0, size=0.33) +
  geom_segment(aes(yend=0,xend=date), size=0.33) +
  # geom_smooth( aes(group=1), method="lm") +
  scale_color_manual( values = c('blue','red') ) +
  theme_bw() + theme( legend.position = "none") +
  xlab("Date") + ylab(expression(paste("Monthly temperature anomaly (",degree,"C)"))) +
  coord_cartesian( ylim=c(-3.5,3.5))
ggsave(plot = a, "R/Figs/lighthouse_anomaly_monthly_temperature.svg",width=6, height=4)

dsal <- d %>% mutate( signsal = as.character(as.numeric(sal.anom>0)) ) %>% filter( site != "addenbroke")
b <- ggplot( data = dsal, aes( x=date, y=sal.anom,col=signsal)) +
  facet_wrap(~site, ncol=1) +
  geom_vline(xintercept = survey.dates$date, lty=2, size=0.1) +
  geom_hline(yintercept = 0, size=0.33) +
  geom_segment(aes(yend=0,xend=date), size=0.33) +
  # geom_smooth( aes(group=1), method="lm") +
  scale_color_manual( values = c('blue','red') ) +
  theme_bw() + theme( legend.position = "none") +
  xlab("Date") + ylab("Monthly salinity anomaly") +
  coord_cartesian( ylim=c(-3,3))


cowplot::plot_grid( a,b, ncol=1, rel_heights = c(3,2))
ggsave("R/Figs/lighthouse_anomaly_monthly_temperature_salinity.svg",width=6, height=8)
#
  
  
  
  
# show recentmonthly temperature anomalies
drecent <- drecent %>% mutate( sign = as.character(as.numeric(temp.anom>0)) )
# all dates
alldates <- data.frame( date = seq( ymd("2011-07-01"),ymd("2019-10-31"), by = "months"), total = 1 )
drecentall <-left_join(alldates,drecent)
drecentNA <- drecentall %>% filter( is.na(temp.anom))
monthly <- ggplot( data = drecent, aes( x=date, y=temp.anom,col=sign)) +
  facet_wrap(~Site, ncol=1) +
  # geom_vline(xintercept = survey.dates$date, lty=2, size=0.1) +
  geom_hline(yintercept = 0, size=0.33) +
  geom_segment(aes(yend=0,xend=date), size=2) +
  geom_point( data=drecentNA, aes(y=0), shape= 4,col='black') +
  geom_text(aes(x=ymd("2013-01-01"),y = -2.5,label=Site), col='black', fontface="plain",hjust = 0, size=7) +
  # geom_smooth( aes(group=1), method="lm") +
  scale_color_manual( values = c('blue','red') ) +
  xlab("") + ylab(expression(paste("Temperature anomaly (",degree,"C)"))) +
  coord_cartesian( ylim=c(-3.1,3.1), xlim=c(ymd("2011-07-01"),ymd("2019-10-31"))) +
  theme_minimal() + theme(legend.position="none",
                          panel.border = element_rect(colour = "black", fill=NA),
                          text = element_text(size = 20),strip.background = element_blank(),
                          strip.text.x = element_blank())
ggsave("R/Figs/lighthouse_anomaly_monthly_temperature_survey.svg",width=8, height=2*2)
  
  
  
  
  
  
  
  
  
  

# dates for events
# 1997-1998 El Nino
warm_times <- lubridate::ymd(c("1997-06-01","1998-06-01","201408-01","2016-12-31"))

png(file="R/Figs/BC_Lightstation_PCA1_daily_timeseries.png", res = 600, width = 6, height = 1.5, units = "in")
par( mar=c(2,4,0,1)+0.1, mfrow=c(1,1), las = 1 )
plot(x = drecent$date, y = pcscores$pca1, type = 'l', col = 'slateblue', ylab = "PCA1" ); abline(h = 0)
lines( lowess(x = drecent$date, y = pcscores$pca1, f = 1/20, iter = 10), col = "darkslateblue", lwd=2)
abline( v = warm_times, lty = 4, col = "slategrey" )
axis(1, at = lubridate::ymd(paste0(2012:2019,"-01-01")), labels = FALSE, col='magenta' )
dev.off()


dna <- bind_cols( drecent, pcscores )
dna$month = lubridate::month(dna$date)
dna$year = lubridate::year(dna$date)

# calculate winter and summer temperature anomalies
library(zoo)
yq <- as.yearqtr( as.yearmon( paste(dna$month,dna$year,sep="/"), "%m/%Y") + 1/12)
dna$season <- factor(format(yq, "%q"), levels = 1:4,
                     labels = c("winter", "spring", "summer", "fall"))
# add year groupings - for instance, seaweeds in summer 2016 would be influenced by conditions over the previous year, 
#                      so, count previous summer, fall, and current winter and spring towards a give year
dna$survey.year <- dna$year 
dna$survey.year[ dna$month %in% 6:12 ] <- dna$survey.year[ dna$month %in% 6:12 ] + 1
dna %>% select(date, season, survey.year, temp_pine ) %>% filter( survey.year > 2010)
anoms.season <- dna %>%
  group_by( survey.year, season ) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
ggplot(anoms.season, aes(x=survey.year,y=pca2,col=season)) + geom_line()
ggplot(anoms.season, aes(x=survey.year,y=pca1,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
ggplot(filter(anoms.season, survey.year>=2010), aes(x=survey.year,y=pca1,col=season)) + geom_path() + geom_point()

# include previous year as summer to summer rather than winter to winter
# anoms.season$yeargroup <-  anoms.season$year + 1
# anoms.season$yeargroup <-  ifelse(anoms.season$season %in% c("winter","spring"), anoms.season$yeargroup - 1, anoms.season$yeargroup)
# compare with means across years
anoms.annual <- anoms.season %>%  # or dna if just using annual means
  group_by( survey.year ) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)
ggplot(filter(anoms.annual, year>=2010 & survey.year <2020), aes(x=survey.year,y=pca1)) + geom_path() + geom_point()
ggplot(filter(anoms.annual, year>=2010), aes(x=survey.year,y=pca2)) + geom_path() + geom_point()
ggplot(filter(anoms.annual, year>=2010), aes(x=survey.year,y=pca3)) + geom_path() + geom_point()
ggplot(filter(anoms.annual), aes(x=survey.year,y=sal_mccinnis)) + geom_hline(yintercept = 0) +  geom_path(col = 'slateblue') + geom_point(col = 'slateblue') + ylim(c(-1.75,1.75))
ggplot(filter(anoms.annual), aes(x=survey.year,y=sal_pine)) + geom_hline(yintercept = 0) +  geom_path(col = 'slateblue') + geom_point(col = 'slateblue') + ylim(c(-1.75,1.75))
ggplot(filter(anoms.annual, year>=2010 &year<2020), aes(x=year,y=temp_pine)) + geom_path() + geom_point()


# extract data for 2011 to 2019
as.survey <- anoms.annual %>% 
  filter( survey.year>=2010 ) 

