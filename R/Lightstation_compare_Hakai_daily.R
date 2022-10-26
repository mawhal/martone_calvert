### This script compares data from BC lightstations to data from 
### the Hakai Institute
# started by Matt Whalen on 4 Nov 2021

## packages
library(tidyverse)
library(lubridate)


## data files
# lightstation - from "Lightstation_spectral_anomaly"
light <- read_csv( "Data/R code for Data Prep/output from R/Lightstation_raw.csv" )
# Hakai data - Choked inner seagrass and West Beach kelp
hakai <- read_csv( "Data/environmental_data/Hakai_nearshore/2021-11-04.1daySamples.csv",
                   skip = 3)


## prep Hakai data columns
hakai <- hakai %>% select( date, year, month, choked = WaterTemp_Avg...10, west = WaterTemp_Avg...18)


## calculate monthly averages and average the two data sources
hdaily <- hakai %>% 
  group_by( year, month, date ) %>% 
  summarize( west = mean(west,na.rm = T), choked = mean(choked,na.rm = T) ) %>% 
  pivot_longer( cols = c(choked, west), names_to = "temp") %>% 
  group_by( year, month, date) %>% 
  summarise( temp = mean(value, na.rm = T)) %>% 
  mutate( month = match(month, month.abb))
hdaily$site = "Calvert"

ldaily <- light %>%
  mutate( year=year(date), month=month(date) ) %>% 
  group_by(year, month, date ,site) %>%
  summarize( temp=mean(temp,na.rm=T) )%>% 
  mutate( site = Hmisc::capitalize(site))
  # summarize( temp=mean(temp,na.rm=T),sal=mean(sal,na.rm=T),
  #            temp_max = mean(temp_max,na.rm=T),
  #            precip = mean(precip, na.rm=T) )

# combine, pivot longer, and filter months with Hakai data
tempcombine <- bind_rows(hdaily,ldaily) %>% arrange( year, month, date) %>% 
  pivot_wider( names_from = site, values_from = temp ) %>% 
  filter( !is.na(Calvert)) 

# pairs plot
svg( "R/Figs/lighstation_hakai_pairs.svg", width = 4, height = 4 )
psych::pairs.panels( tempcombine[,-c(1,2,3)], density = F, hist.col = "whitesmoke",
                     ylim=c(-10,18), xlim=c(-10,18) )
dev.off()

# pivot longer for time series plots
temp_long <- tempcombine %>% pivot_longer( cols = c(Pine, Mcinnes, Addenbroke, Calvert), names_to = "site", values_to = "temp" )

ggplot( temp_long %>% filter(site != "Addenbroke"), aes(x = date, y = temp, col = site)) + 
  # facet_wrap( ~site, ncol = 1) +
  geom_line(aes(lwd = site), alpha = 0.5) +
  ylab( expression(paste("Temperature (", degree, "C)")) ) +
  xlab("Date") +
  scale_size_manual(values = c(1,0.5,0.5)) +
  theme_classic() + 
  theme(legend.position="top")
ggsave( "R/Figs/lightstation_hakai_timeseries.svg", width = 4, height = 4 )

range( temp_long$date )
