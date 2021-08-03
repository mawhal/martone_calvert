# tide predictions from Adams Harbour
# example
# https://hecate.hakai.org/api/modelling/tide_prediction.csv?station=ADAMS_HARBOUR&date=2011-06-05T00:00:00.000Z&duration=30
# 

# tide predictions taken one month around the summer solstice and median sampling date in each sampled year (2011-2019)
# 

library(tidyverse)
library(lubridate)

# read data and bind together
files <- list.files(path = "Data/tides/predictions_adams_harbour/", full.names = TRUE)
dlist <- lapply( files, read_csv )
dbind <- do.call( rbind, dlist)

# add columns for year, prediction type (solstice or sampling)
d <- dbind %>% 
  mutate( year = lubridate::year(date),
          source = rep( c("sampling","solstice"), 
                        each = nrow(dbind)/length(files),
                        length.out = nrow(dbind) ) )

# scale date
d <- d %>% 
  group_by(year,source) %>% 
  mutate( date_scale = scale(date))

# plot
ggplot( data = d, aes(x = date_scale,y = sea_surface_elevation)) +
  facet_grid(year~source) +
  geom_vline(xintercept = 0, col = 'red' ) +
  geom_line(aes(col=source),lwd=0.33) +
  scale_color_manual(values = c("slategrey","slateblue")) +
  scale_y_continuous(breaks = c(1,3,5), minor_breaks = c(0:5)) +
  ylab("Sea surface elevation (m above MLLWLT)") + xlab("Time") +
  theme_classic() +
  theme(legend.position = 'none')
ggsave( "Data/tides/tides_sampling_solstice.svg", width = 5.5, height = 5)
