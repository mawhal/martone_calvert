# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen 
# created 20 May 2022

# comparison of seaweed species diversity and cover through time


library(tidyverse)



## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
bare_raw <- ad %>% 
  filter(taxon_lumped3 == "Bare rock") %>% 
  select(UID,taxon_lumped3, Abundance) %>% 
  group_by(UID, taxon_lumped3) %>%
  summarize( cover = sum(Abundance))
  
  
# all metadata
am <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011" & am$Year <= 2019, ]
# remove Meay Channel
am <- am[ am$Site != "Meay Channel", ]

# bare rock for survey
bare <- left_join( am, bare_raw )
bare$cover[ is.na(bare$cover) ] <- 0
bare <- bare %>% 
  unite("transect", Site, Zone, sep=" ", remove=F) %>% 
  select(Year,Site,Zone,transect,cover)
transect.bare <- bare %>% 
  group_by(Year, Site, Zone, transect ) %>% 
  summarize( cover.bare = mean(cover) )

# read in data from Community_rda.R
d.simple <- read_csv("R/output/data_select_rda_HMSC.csv")

# define quadrats
d.simple <- d.simple %>% 
  mutate( quadrat = paste( Site, Zone, Meter.point) )

# using all years, unlike some of the other code looking at stability (e.g., Wang approach)
d.select.all.years <- d.simple


# just get algae, and calculate richness across sites
d.comm.algae.site <- d.select.all.years %>%
  filter( funct_2021 != "animal" ) %>% 
  select(-funct_2021) %>% 
  group_by( Year, Site,  taxon ) %>% 
  summarize( Abundance = sum(Abundance) ) %>% 
  ungroup() %>% 
  spread( taxon, Abundance, fill=0 ) 
d.comm.richness <- d.comm.algae.site %>% 
  select( - Year, -Site ) 
d.richness <- rowSums( d.comm.richness>0 )  
transect.richness <- bind_cols( select(d.comm.algae.site, Year, Site), richness=d.richness ) 
transect.richness$Site <- factor(transect.richness$Site, levels = c("Foggy Cove","Fifth Beach","North Beach") )
ggplot( transect.richness, aes(x=Year, y=richness, group=Site )) + 
  geom_line( ) +   geom_point( aes(shape=Site), bg="white", size=3 ) +
  scale_shape_manual( values = c(21,22,24) ) +
  theme_classic() + theme( legend.position = "none")

# just get algae, and calculate richness along transects
d.comm.algae.transect <- d.select.all.years %>%
  filter( funct_2021 != "animal" ) %>% 
  select(-funct_2021) %>% 
  group_by( Year, Site, Zone, taxon ) %>% 
  summarize( Abundance = sum(Abundance) ) %>% 
  ungroup() %>% 
  spread( taxon, Abundance, fill=0 ) 
d.comm.richness <- d.comm.algae.transect %>% 
  select( - Year, -Site, -Zone )
d.richness <- rowSums( d.comm.richness>0 )  
transect.richness <- bind_cols( select(d.comm.algae.transect, Year, Site, Zone), richness=d.richness ) %>%
  unite( "transect", Site, Zone, sep=" ", remove = F)
transect.richness$Zone <- factor(transect.richness$Zone, levels = c("LOW","MID","HIGH") )
transect.richness$Site <- factor(transect.richness$Site, levels = c("Foggy Cove","Fifth Beach","North Beach") )
ggplot( transect.richness, aes(x=Year, y=richness )) + 
  geom_smooth() +
  geom_point( size=3 ) +
  ylab("Total transect species richness") +
  theme_classic() + theme( legend.position = "none")
ggplot( transect.richness, aes(x=Year, y=richness, group=Zone )) + 
  facet_wrap(~Site, ncol=1 ) +
  geom_line( aes(lty=Zone, col=Zone)) +   geom_point( aes(col=Zone)) +
  scale_color_manual( values = c("black","gray","gray") ) +
  scale_linetype_manual( values = c(1,1,2) ) +
  # scale_shape_manual( values = c(21,22,24) ) +
  ylab("Total transect species richness") +
  theme_classic() + theme( legend.position = "none")
ggplot( transect.richness, aes(x=Year, y=richness,  pch=Site, group=transect )) + 
  geom_line( aes(lty=Zone, col=Zone)) +   geom_point( aes(col=Zone, fill=Zone), size=2) +
  scale_color_manual( values = c("black","gray","gray") ) +
  scale_fill_manual( values = c("black","gray","white") ) +
  scale_linetype_manual( values = c(1,1,2) ) +
  scale_shape_manual( values = c(21,22,24) ) +
  ylab("Total transect species richness") +
  theme_classic() + theme( legend.position = "none")

ggplot( transect.richness, aes(x=Year, y=richness,   group=transect )) + 
  # geom_line( aes(lty=Zone, col=Zone)) +   geom_point( aes(col=Zone, fill=Zone), size=2) +
  geom_smooth( aes(lty=Zone, col=Zone), method = 'lm', se=F) +
  scale_color_manual( values = c("black","gray","gray") ) +
  scale_fill_manual( values = c("black","gray","white") ) +
  scale_linetype_manual( values = c(1,1,2) ) +
  scale_shape_manual( values = c(21,22,24) ) +
  ylab("Total transect species richness") +
  coord_cartesian(ylim=c(0,60)) +
  theme_classic() + theme( legend.position = "none")






# make a smaller plot showing total seaweed cover over time for each transect
# first get total seaweed cover in each quadrat then average within transect
# names(comm)
d.comm.algae.quad <- d.select.all.years %>%
  filter( funct_2021 != "animal" ) %>% 
  select(-funct_2021) %>% 
  group_by( Year, Site, Zone, quadrat, taxon ) %>% 
  summarize( Abundance = sum(Abundance) ) %>% 
  ungroup() %>% 
  spread( taxon, Abundance, fill=0 ) 

d.comm.cover <- d.comm.algae.quad %>% 
  select( - Year, -Site, -Zone, -quadrat )
d.cover <- rowSums( d.comm.cover ) 
quadrat.cover <- bind_cols( select(d.comm.algae.quad, Year, Site, Zone, quadrat), cover=d.cover ) %>%
  unite( "transect", Site, Zone, sep=" ", remove = F)
transect.cover <- quadrat.cover %>% 
  group_by(Year, Site, Zone, transect ) %>% 
  summarize( cover = mean(cover) )

transect.cover$Site <- factor( transect.cover$Site, levels = c("Foggy Cove","Fifth Beach","North Beach"))
transect.cover$Zone <- factor( transect.cover$Zone, levels = c("LOW","MID","HIGH"))

# repeat for animals
d.comm.invert.quad <- d.select.all.years %>%
  filter( funct_2021 == "animal" ) %>% 
  select(-funct_2021) %>% 
  group_by( Year, Site, Zone, quadrat, taxon ) %>% 
  summarize( Abundance = sum(Abundance) ) %>% 
  ungroup() %>% 
  spread( taxon, Abundance, fill=0 ) 

d.comm.invert.quad %>% 
  select(Year, Site, Zone) %>% distinct()

d.comm.cover.invert <- d.comm.invert.quad %>% 
  select( - Year, -Site, -Zone, -quadrat )
d.cover <- rowSums( d.comm.cover.invert ) 
quadrat.cover.invert <- bind_cols( select(d.comm.invert.quad, Year, Site, Zone, quadrat), cover=d.cover ) %>%
  unite( "transect", Site, Zone, sep=" ", remove = F)
transect.cover.invert <- quadrat.cover.invert %>% 
  group_by(Year, Site, Zone, transect ) %>% 
  summarize( cover.invert = mean(cover) )

transect.cover.invert$Site <- factor( transect.cover.invert$Site, levels = c("Foggy Cove","Fifth Beach","North Beach"))
transect.cover.invert$Zone <- factor( transect.cover.invert$Zone, levels = c("LOW","MID","HIGH"))




# compare trends in richness and cover
transect.cover.richness <- left_join(left_join(left_join(transect.cover,transect.cover.invert), transect.richness), transect.bare )
transect.cover.richness$point.color <- transect.cover.richness$Year
transect.cover.richness$point.color[transect.cover.richness$Year==2012] <- "white"
transect.cover.richness$point.color[transect.cover.richness$Year==2019] <- "black"
transect.cover.richness$Site <- factor( transect.cover.richness$Site, levels = c("Foggy Cove","Fifth Beach","North Beach"))
transect.cover.richness$Zone <- factor( transect.cover.richness$Zone, levels = c("LOW","MID","HIGH"))
transect.cover.richness.points <- transect.cover.richness[ transect.cover.richness$Year %in% c(2012,2019), ]


# add heatwave duration from script 
heatwave_duration_pca <- read_csv("R/output/heatwaveR_duration_surveyyear_pca.csv")
cover.richness.heatwave <- left_join( transect.cover.richness, heatwave_duration_pca ) %>% ungroup()
cover.richness.heatwave$cover.invert[ is.na(cover.richness.heatwave$cover.invert) ] <- 0


psych::pairs.panels( select(cover.richness.heatwave, richness, cover, duration, intensity_max) )




# correlations and linear models

# library(broom)
# reglines_cover <- cover.richness.heatwave %>% 
#   nest(data = -transect) %>% 
#   mutate(
#     test = map(data, ~ lm(log(cover)~Year, data=.x)), # S3 list-col
#     tidied = map(test, tidy)
#   ) %>% 
#   unnest(tidied)
# reg_cover_sig <- reglines_cover %>% filter(term=="Year") %>% filter( p.value < 0.05 )
# 
# reglines_invert <- cover.richness.heatwave %>% 
#   nest(data = -transect) %>% 
#   mutate(
#     test = map(data, ~ lm(log(cover.invert+0.5)~Year, data=.x)), # S3 list-col
#     tidied = map(test, tidy)
#   ) %>% 
#   unnest(tidied)
# reg_invert_sig <- reglines_invert %>% filter(term=="Year") %>% filter( p.value < 0.05 )
# 
# reglines_richness <- cover.richness.heatwave %>% 
#   nest(data = -transect) %>% 
#   mutate(
#     test = map(data, ~ lm((richness)~Year, data=.x)), # S3 list-col
#     tidied = map(test, tidy)
#   ) %>% 
#   unnest(tidied)
# reg_richness_sig <- reglines_richness %>% filter(term=="Year") %>% filter( p.value < 0.05 )
# 
# reglines_bare <- cover.richness.heatwave %>% 
#   nest(data = -transect) %>% 
#   mutate(
#     test = map(data, ~ lm(log(cover.bare)~Year, data=.x)), # S3 list-col
#     tidied = map(test, tidy)
#   ) %>% 
#   unnest(tidied)
# reg_bare_sig <- reglines_bare %>% filter(term=="Year") %>% filter( p.value < 0.05 )

# regressions <- transect.cover %>%
#   nest(data = -transect) %>% 
#   mutate(
#     fit = map(data, ~ lm((cover)~Year, data = .x)),
#     tidied = map(fit, tidy),
#     glanced = map(fit, glance),
#     augmented = map(fit, augment)
#   )






psych::pairs.panels( select(cover.richness.heatwave, richness, cover, cover.invert) )
psych::pairs.panels( select(cover.richness.heatwave, Year, duration, duration5) )

summary( lm( log(cover)~duration, data=cover.richness.heatwave ) )
library(lme4)
library(lmerTest)
library(bbmle)
m11 <- ( lmer( log(cover)~duration + (1|transect), data=cover.richness.heatwave ) )
m12 <- ( lmer( log(cover)~duration5 + (1|transect), data=cover.richness.heatwave ) )
m13 <- ( lmer( log(cover)~Year + (1|transect), data=cover.richness.heatwave ) )
m14 <- ( lmer( log(cover)~scale(duration5)*scale(duration) + (1|transect), data=cover.richness.heatwave ) )
summary(m14)
AICctab( m11, m12, nobs = nrow(cover.richness.heatwave))
m21 <- ( lmer( log(cover.invert+1)~duration + (1|transect), data=cover.richness.heatwave ) )
m22 <- ( lmer( log(cover.invert+1)~duration5 + (1|transect), data=cover.richness.heatwave ) )
m23 <- ( lmer( log(cover.invert+1)~Year + (1|transect), data=cover.richness.heatwave ) )
m24 <- ( lmer( log(cover.invert+1)~scale(duration5)*scale(duration) + (1|transect), data=cover.richness.heatwave ) )
AICctab( m21, m22,  nobs = nrow(cover.richness.heatwave))
m31 <- ( lmer( log(cover.bare+1)~duration + (1|transect), data=cover.richness.heatwave ) )
m32 <- ( lmer( log(cover.bare+1)~duration5 + (1|transect), data=cover.richness.heatwave ) )
m33 <- ( lmer( log(cover.bare+1)~Year + (1|transect), data=cover.richness.heatwave ) )
m34 <- ( lmer( log(cover.bare+1)~scale(duration5)*scale(duration) + (1|transect), data=cover.richness.heatwave ) )
AICctab( m31, m32,  nobs = nrow(cover.richness.heatwave))
m41 <- ( lmer( log(richness)~duration + (1|transect), data=cover.richness.heatwave ) )
m42 <- ( lmer( log(richness)~duration5 + (1|transect), data=cover.richness.heatwave ) )
m43 <- ( lmer( log(richness)~Year + (1|transect), data=cover.richness.heatwave ) )
m44 <- ( lmer( log(richness)~scale(duration5)*scale(duration) + (1|transect), data=cover.richness.heatwave ) )
AICctab( m41, m42,  nobs = nrow(cover.richness.heatwave))

psych::pairs.panels( select( cover.richness.heatwave, 
                             duration5, duration,
                             cover, cover.bare) %>% 
                       mutate(cover=log(cover),cover.bare=log(cover.bare+1)),
                    ellipses = F, lm = T, ci=T )
cor.test( log(cover.richness.heatwave$cover), log(cover.richness.heatwave$cover.bare+1))

# for presenation use width = 4, height = 4
windows(5,5)
psych::pairs.panels( select( cover.richness.heatwave, 
                             duration5, 
                             cover, cover.invert, cover.bare) %>% 
                       mutate(seaweed.cover=(cover+1),invert.cover=(cover.invert+1),bare.rock.cover=(cover.bare+1)) %>% 
                       select( duration5, seaweed.cover, invert.cover, bare.rock.cover),
                     ellipses = F, lm = T, ci=T, log = "xy" )

windows(5,5)
psych::pairs.panels( select( cover.richness.heatwave, 
                             duration5, 
                             cover, cover.invert, cover.bare) %>% 
                       mutate(seaweed.cover=log10(cover),invert.cover=log10(cover.invert+1),bare.rock.cover=log10(cover.bare+1)) %>% 
                       select( duration5, seaweed.cover, invert.cover, bare.rock.cover),
                     ellipses = F, lm = F, ci = T, smooth = T, stars = T, cex = 1.5, cex.cor = 1 )




m11 <- ( glmer( ceiling(cover)~duration + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m12 <- ( glmer( ceiling(cover)~duration5 + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m13 <- ( glmer( ceiling(cover)~Year + (1|transect),  family="poisson", data=cover.richness.heatwave ) )
AICctab( m11, m12, m13, nobs = nrow(cover.richness.heatwave))
m21 <- ( glmer( ceiling(cover.invert)~duration + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m22 <- ( glmer( ceiling(cover.invert)~duration5 + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m23 <- ( glmer( ceiling(cover.invert)~Year + (1|transect),  family="poisson", data=cover.richness.heatwave ) )
AICctab( m21, m22, m23, nobs = nrow(cover.richness.heatwave))
m31 <- ( glmer( ceiling(cover.bare)~duration + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m32 <- ( glmer( ceiling(cover.bare)~duration5 + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m33 <- ( glmer( ceiling(cover.bare)~Year + (1|transect),  family="poisson", data=cover.richness.heatwave ) )
AICctab( m31, m32, m33, nobs = nrow(cover.richness.heatwave))
m41 <- ( glmer( ceiling(richness)~duration + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m42 <- ( glmer( ceiling(richness)~duration5 + (1|transect), family="poisson", data=cover.richness.heatwave ) )
m43 <- ( glmer( ceiling(richness)~Year + (1|transect),  family="poisson", data=cover.richness.heatwave ) )
AICctab( m41, m42, m43, nobs = nrow(cover.richness.heatwave))

mduration <- list(m11,m21,m31,m41)
mduration5 <- list(m12,m22,m32,m42)
do.call( rbind, lapply(mduration, function(z) fixef(z) ) )
do.call( rbind, lapply(mduration5, function(z) fixef(z) ) )
summary(m11)
summary(m21)
summary(m31)
summary(m41)
summary(m12)
summary(m22)
summary(m32)
summary(m42)


### model predictions
### cover
## model 1 - duration
exp( 3.4319467 + (200*-0.0002996) ) - exp( 3.4319467 + (0*-0.0002996) ) # richness
exp( 4.681079 + (200*-0.0019957079) ) - exp( 4.681079 + (0*-0.0019957079) ) # seaweed cover
exp( 2.4793016  + (200*0.0002680) ) - exp( 2.4793016  + (0*0.0002680) ) # invert cover
exp( 2.4059740 + (200*0.0023670) ) - exp( 2.4059740 + (0*0.0023670) ) # bare rock
## model 2 - five-year duration
exp( 3.425e+00 + (500*-4.906e-05) ) - exp( 3.425e+00 + (0*-4.906e-05) ) # richness
exp( 4.759e+00 + (500*-8.891099e-04) ) - exp( 4.759e+00 + (0*-8.891099e-04) ) # seaweed cover
exp( 2.3836987 + (500*0.0004574) ) - exp( 2.3836987 + (0*0.0004574) ) # invert cover
exp( 2.2190207 + (500*0.0013231) ) - exp( 2.2190207 + (0*0.0013231) ) # bare rock



## get predictions for plotting seaweed cover and bare cover as a function of duration and duration5?
#


# 
# pred2019 <- unlist(lapply( regressions$augmented, function(z) z$.fitted[ z$Year==2019] ))
# pred2019 <- data.frame( pred2019, Year = 2019, transect = regressions$transect )
# pred2019$code <- c( rep("5B",3), rep("FC",3), rep("NB",3))
# pred2019 <- pred2019 %>% 
#   separate( transect, c("Beach1", "Beach2", "Zone"), sep = " ", remove = FALSE) %>% 
#   unite( "Site", Beach1, Beach2, sep = " ")
# prod_empir_trend_transect <- ggplot( transect.cover, aes(x = Year, y = (cover), group = transect,
#                                                 col = Zone, lty = Zone)) +
#   geom_smooth( method = 'lm', se = F,   alpha=0.5, lwd = 1) +
#   geom_text( data=pred2019, aes( x = Year, y = pred2019, label = code),
#              hjust = "outward", vjust = c(0.5,0.5,0, 1,0.5,0, 1,0.5,0.5),
#              size = 2.5,show.legend = FALSE) +
#   theme_classic() +
#   theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#          legend.position = "top",
#          legend.title = element_blank(),
#          legend.text = element_text(size = 7),
#          legend.key.size = unit(0.5, "cm"),
#          legend.key = element_rect(colour = NA, fill = NA),
#          legend.box.margin=margin(-10,-10,-10,-10)) +
#   coord_cartesian( ylim = c(0,180), xlim = c(2012,2019.5) ) +
#   scale_color_manual(values = c("black","grey","grey")) +
#   scale_linetype_manual(values = c(1,1,3)) +
#   guides( lty = guide_legend(ncol=2), lwd = guide_legend(ncol=3,byrow = FALSE) ) +
#   ylab("Mean seaweed cover (%)") + xlab("Year")
# prod_empir_trend_transect
# # ggsave("R/Figs/seaweed_cover_time_transect.svg", width = 3, height = 3)

lm1 <- lm( log(cover) ~ Year*transect, data = cover.richness.heatwave )
summary(lm1)
anova(lm1)
trend.tab <- emmeans::lstrends( lm1, ~ transect, var = "Year" )
trend.cover <- trend.tab %>% as.data.frame() %>%  mutate( sig = (lower.CL*upper.CL)>0 ) %>% 
  filter( sig == T )
lm2 <- lm( log(cover.invert+1) ~ Year*transect, data = cover.richness.heatwave )
summary(lm2)
anova(lm2)
trend.tab <- emmeans::lstrends( lm2, ~ transect, var = "Year" )
trend.invert <- trend.tab %>% as.data.frame() %>%  mutate( sig = (lower.CL*upper.CL)>0 ) %>% 
  filter( sig == T )
lm3 <- lm( log(cover.bare) ~ Year*transect, data = cover.richness.heatwave )
summary(lm3)
anova(lm3)
trend.tab <- emmeans::lstrends( lm3, ~ transect, var = "Year" )
trend.bare <- trend.tab %>% as.data.frame() %>%  mutate( sig = (lower.CL*upper.CL)>0 ) %>% 
  filter( sig == T )
lm4 <- lm( log(richness) ~ Year*transect, data = cover.richness.heatwave )
summary(lm4)
anova(lm4)
trend.tab <- emmeans::lstrends( lm4, ~ transect, var = "Year" )
trend.richness <- trend.tab %>% as.data.frame() %>%  mutate( sig = (lower.CL*upper.CL)>0 ) %>% 
  filter( sig == T )



## Plot everything together
richness.keep <- cover.richness.heatwave %>% filter( transect %in% trend.richness$transect  )
a <- ggplot( cover.richness.heatwave, aes(x=Year,y=richness, col=Zone, lty=Zone, group=transect) ) + 
  facet_wrap(~Site) +
  geom_path() +
  geom_smooth(data = richness.keep, method="lm", se=F, lty=1 ) +
  geom_point(aes(size=duration,fill=Zone), shape=21,col='black') + 
  ylab("Seaweed species richness") +
  scale_color_manual( values = c('black','grey','grey85')) +
  scale_fill_manual( values = c('black','grey','white')) +
  scale_y_log10( limits = c(10,60),breaks = c(10,15,20,30,40,50)) +
  theme_classic( ) +
  theme( legend.position = "none")
cover.keep <- cover.richness.heatwave %>% filter( transect %in% trend.cover$transect  )
b <- ggplot( cover.richness.heatwave, aes(x=Year,y=cover, col=Zone, lty=Zone, group=transect) ) + 
  facet_wrap(~Site) +
  geom_hline(yintercept = 100) +
  geom_path() +
  # geom_smooth(data = cover.keep,method = "glm", formula = y ~ x, se = F,  lty = 1, size=1,
  #             method.args = list(family=gaussian(link="log")) ) +
  geom_smooth(data = cover.keep,method='lm', se=F, lty=1) +
  geom_point(aes(size=duration,fill=Zone), shape=21,col='black') + 
  ylab("Seaweed % cover") +
  scale_color_manual( values = c('black','grey','grey85')) +
  scale_fill_manual( values = c('black','grey','white')) +
  scale_y_log10( limits = c(13,210),breaks = c(13,25,50,100,150,200)) +
  theme_classic( ) +
  theme( legend.position = "right") +
  guides( color = "none", lty = "none", fill = "none", size = guide_legend(title="Heatwave\ndays in\nprevious\nyear")  ) 
invert.keep <- cover.richness.heatwave %>% filter( transect %in% trend.invert$transect  ) 
c <- ggplot( cover.richness.heatwave, aes(x=Year,y=cover.invert, col=Zone, lty=Zone, group=transect) ) + 
  facet_wrap(~Site) +
  geom_path() +
  geom_smooth(data = invert.keep,method='lm', se=F, lty=1) +
  geom_point(aes(size=duration,fill=Zone), shape=21,col='black') + 
  ylab("Invertebrate % cover") +
  scale_color_manual( values = c('black','grey','grey85')) +
  scale_fill_manual( values = c('black','grey','white')) +
  scale_y_log10( limits = c(0.5,100),breaks = c(0.5,1,5,10,25,50,100) ) +
  theme_classic( ) +
  theme( legend.position = "right") +
  guides( size = "none" )
bare.keep <- cover.richness.heatwave %>% filter( transect %in% trend.bare$transect  ) 
d <- ggplot( cover.richness.heatwave, aes(x=Year,y=cover.bare, col=Zone, lty=Zone, group=transect) ) + 
  facet_wrap(~Site) +
  geom_path() +
  geom_smooth(data = bare.keep,method='lm', se=F, lty=1) +
  geom_point(aes(size=duration,fill=Zone), shape=21,col='black') + 
  ylab("Bare rock % cover") +
  scale_color_manual( values = c('black','grey','grey85')) +
  scale_fill_manual( values = c('black','grey','white')) +
  scale_y_log10( limits = c(0.5,100),breaks = c(0.5,1,5,10,25,50,100)) +
  theme_classic( ) +
  theme( legend.position = "none")
cowplot::plot_grid( a,b,c,d, ncol=1, labels = "auto", hjust = -4, align = 'hv', axis='tblr' )
ggsave( "R/Figs/cover.richness.heatwave_trends.svg", width=6, height=8)






