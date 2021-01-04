

# Show tempearture anomalies in different ways to compare years and seasons

# code by Matt Whalen


## temperature anomaly data from Pine Island
anoms <-  read_csv("Data/R code for Data Prep/output from R/PineIsland_monthly_SST_anomaly.csv")

# calculate winter and summer temperature anomalies
library(zoo)
yq <- as.yearqtr( as.yearmon( paste(anoms$month,anoms$year,sep="/"), "%m/%Y") + 1/12)
anoms$season <- factor(format(yq, "%q"), levels = 1:4, 
                       labels = c("winter", "spring", "summer", "fall"))
anoms.season <- anoms %>% 
  group_by( year, season ) %>% 
  summarize( temp.anom=mean(temp.anom) )
ggplot(anoms.season, aes(x=year,y=temp.anom,col=season)) + geom_point()
ggplot(anoms.season, aes(x=year,y=temp.anom,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
ggplot(filter(anoms.season, year>=2010), aes(x=year,y=temp.anom,col=season)) + facet_wrap(~season) + geom_path() + geom_point()
# figure out how to make an anomly plot with vertical lines from zero

# extract data for 2011 to 2019
as.survey <- anoms.season %>% 
  filter( year>=2010 ) %>% 
  spread( key = season, value=temp.anom )

test <- c("Station_1_CTD_42","Station_1_CTD_42_2")
gsub( "Station_1_CTD_42*", "Station_1", test )


# show summer versus winter temperature anomaly
as.survey.all <- anoms.season %>% 
  # filter( year>=2010 ) %>% 
  spread( key = season, value=temp.anom )

ggplot( as.survey, aes(y=summer,x=winter) ) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_smooth(se=T,size=0.5) +
  geom_point(size=2) +  geom_path(size=1) +
  geom_text_repel(label=as.survey$year, box.padding = 0.5) +
  ylab(expression(paste("Summer SST anomaly (",degree,"C)"))) +
  xlab(expression(paste("Winter SST anomaly (",degree,"C)"))) +
  theme_bw()

old.anom <- as.survey.all %>% filter( winter < -1.1 | winter >1 | summer > 1 | summer < -1) %>% filter(year<2015 )
ggplot( as.survey.all, aes(y=summer,x=winter) ) + 
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  stat_ellipse( level = 0.9 ) +
  geom_point(size=1,shape=1) +
  geom_point(data=as.survey,size=2) +  geom_path(data=as.survey,size=0.75) +
  geom_text(data=old.anom, label=old.anom$year,  size=2.5, vjust = c('top','bottom','left','right') ) +
  geom_text(data=filter(as.survey.all,year %in% c(2010,2019)), 
            label=c(2010,2019), nudge_x = 0.25, size=3 ) +
  ylab(expression(paste("Summer anomaly (",degree,"C)"))) +
  xlab(expression(paste("Winter anomaly (",degree,"C)"))) +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
ggsave( "R Code and Analysis/Figs/winter_summer_SST.svg", width=3, height=3 )
