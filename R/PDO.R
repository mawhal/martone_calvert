# Pacific Decadal Oscillation


library(tidyverse)
library(lubridate)

# read raw data
pdo <- read_csv("Data/environmetal_data/PDO/data.csv", skip = 1)
pdo$date <- ym(pdo$Date)


pdo$sign <- ifelse(pdo$Value > 0, "1", "-1")
ggplot( pdo, aes(x = date, y = Value, col = sign)) + 
  geom_point( ) +
  geom_segment( aes(x = date, xend = date, yend = Value, y = 0) ) +
  scale_color_manual( values = c('blue','red') )
  



ggplot( filter(pdo, date > ym(201001)), aes(x = date, y = Value)) + 
  geom_hline( yintercept = 0 ) +
  geom_vline( xintercept = ym( c(201410,201601)), lty=2 ) +
  geom_vline( xintercept = ym( c(201107,201907)), lty=1 ) +
  geom_point( col = 'darkblue' ) +
  geom_path( col = 'blue' ) 
ggplot( filter(pdo, date > ym(201001)), aes(x = date, y = Value, col = sign)) + 
  geom_hline( yintercept = 0 ) +
  geom_vline( xintercept = ym( c(201410,201601)), lty=2 ) +
  geom_vline( xintercept = ym( c(201107,201907)), lty=1 ) +
  geom_point( ) +
  geom_segment( aes(x = date, xend = date, yend = Value, y = 0) ) +
  scale_color_manual( values = c('blue','red') ) +
  theme( legend.position = "none" )






# calculate mean PDO for months leading up to each survey (June to June)
# calculate winter and summer temperature anomalies
# add year groupings - for instance, seaweeds in summer 2016 would be influenced by conditions over the previous year, 
#                      so, count previous summer, fall, and current winter and spring towards a give year
pdo$survey.year <- year(pdo$date) 
pdo$survey.year[ month(pdo$date) %in% 6:12 ] <- pdo$survey.year[ month(pdo$date) %in% 6:12 ] + 1

pdo.survey <- pdo %>% 
  group_by(survey.year) %>% 
  summarize( value = mean(Value) ) %>% 
  filter( survey.year %in% 2012:2019 )

# write to disk
write_csv(pdo.survey, "Data/environmetal_data/PDO/pdo_survey_years.csv")


        