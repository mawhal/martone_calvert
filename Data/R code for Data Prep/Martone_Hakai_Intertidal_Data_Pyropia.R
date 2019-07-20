# load libraries
library(tidyverse)


# look at some Pyropia

all.data <- read.csv( "Output from R/Martone_Hakai_data.csv")
all.meta <- read.csv("Output from R/Martone_Hakai_metadata.csv" )

# replace all commas with periods for Abundance
all.data$Abundance <- as.numeric( gsub( ",", "[.]", all.data$Abundance ) )

# combine all pyropia
pyropia.rows <- grep( "Pyropia", all.data$Taxon )
pyropia <- all.data[ pyropia.rows,]
pyropia.meta <- left_join( pyropia, all.meta, by=c("SiteHeightYear","Quadrat"="Quadrat.No.") )
    
# just look at Pyropia abbotiae
abbottiae <- all.data[ all.data$Taxon=="Pyropia abbottiae",]

# merge meta data
ab.meta <- left_join( abbottiae, all.meta, by=c("SiteHeightYear","Quadrat"="Quadrat.No.") )

# just look in the low zone
ab.low <- ab.meta[ ab.meta$Zone=="LOW", ]


# North Beach Mid
ab.north <- ab.meta[ ab.meta$Site=="North Beach", ]

# time trends in different tidal heights
ggplot( ab.north, aes(x=as.factor(Year),y=Abundance)) + facet_wrap(~Zone) + geom_boxplot() + geom_point()


# all sites
# time trends in different tidal heights
ggplot( ab.meta, aes(x=as.factor(Year),y=Abundance)) + facet_grid(Site~Zone) + geom_boxplot() + geom_point()
ggplot( pyropia.meta, aes(x=as.factor(Year),y=Abundance)) + facet_grid(Site~Zone) + geom_boxplot() + geom_point()



# save Pyropia data to file, then share with Jenn Clark
write.csv( pyropia.meta,  "Output from R/Martone_Hakai_data_Pyropia.csv", row.names=FALSE)
