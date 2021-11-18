# quadrat overlap across years

# only need metadata for this

library( tidyverse )


## read data files
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

# remove Meay Channel and 2011
am <- am %>% filter( Site != "Meay Channel") %>% filter( Year != 2011 )

# define quadrats by adding Site, Zone, and meter.point
am$quadrat <- with(am, paste(Site, Zone,Meter.point))
  
# make table to compare years
with(am, table(Year,quadrat) )
# create groups of years (2012,2013,2018,2019)
am.prepost <- filter(am, Year %in% c(2012,2013,2018,2019))
am.prepost$prepost <- ifelse( am.prepost$Year %in% c(2012,2013), "pre", "post")

with( am.prepost, table(quadrat,Year) )
prepostquads <- with( am.prepost, table(quadrat,prepost) )
prepostquads <- as.data.frame( prepostquads )
prepostquads <- prepostquads %>% 
  pivot_wider( names_from = prepost, values_from = Freq) 
prepostquads <- prepostquads %>%
  mutate(both = pre * post) %>% mutate( both = both > 0) %>% 
  filter(both == TRUE)

# add back in other metadata so we can look at distribution across transects
prepostquads <- left_join( prepostquads, 
                           distinct(select(am.prepost, 
                                  Site, Zone, quadrat))  )

prepostquads %>% 
  group_by(Site, Zone) %>% 
  summarize( n = length(pre) )


# write to disk
write_csv( prepostquads, "R/output/quadrats_pre_post_heatwave.csv")

# next steps: 
# grab total seaweed cover for each and 
# our estimate of species richness
