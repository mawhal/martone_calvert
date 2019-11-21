# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses HMSC to model communities change over time and space

## load libraries
library( tidyverse )
library( here )
library( Hmsc )
library( corrplot )
library( viridis )

## useful references for this code
citation( "Hmsc" )
# https://github.com/hmsc-r/HMSC

## directories
ModelDir = paste0( here::here(), "/R Code and Analysis/models" )
MixingDir = paste0( here::here(), "/R Code and Analysis/mixing")



## load the model
list.files( ModelDir )
model = 7
mload <- load( paste(ModelDir,list.files( ModelDir ), sep="/")[model] )

# load the data
metacomm <- read_csv("R Code and Analysis/output from r/community.csv" )

##
## MCMC convergence

# We evaluate MCMC convergence in terms of four kinds of parameters that we are especially interested in: the species niches `Beta`, the influence of traits on species niches `Gamma`, residual species associations `Omega`, and the strength of phylogenetic signal `rho`. As there are `r ns` species, the matrix `Omega` is a `r ns` x `r ns` matrix with `r ns*ns` elements. To avoid excessive computational times, we evaluate MCMC convergence for `Omega` only for a subsample of 100 randomly selected species pairs.
mpost = convertToCodaObject(m)
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
hist(ess.beta)
hist(psrf.beta)
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma)
ns=50
sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)




## Assess model fit
preds = computePredictedValues(m)
MF <- evaluateModelFit(hM = m, predY = preds)
str(MF)
lapply( MF, summary)
par( mfrow=c(3,3), mar=c(3,3,1,1)+0.1 )
lapply( MF, boxplot )
dev.off()



## parameter estimates
postBeta = getPostEstimate(m, parName = "Beta")
# windows(5,8)
# plotBeta(m, post = postBeta, param = "Sign", supportLevel = 0.95, mar=c(7,11,0,6))
postBeta$mean[, c("Alaria.marginata","Hedophyllum.sessile")]

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "shore.elevation", "shore.elev.sq", "year",
                              "shore.year", "shore2.year" ),
                            levels = c("intercept", "shore.elevation", "shore.elev.sq",  
                                       "year", "shore.year", "shore2.year" ), 
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 6), 
                          levels = colnames(m$Y)[order(colSums(m$Y),decreasing = TRUE)],
                          ordered = TRUE)

windows(12,4)
ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")

# pull taxa with positive estimated response to temperature anomaly
taxa.pos.anom <- as.character(pos.neg[ pos.neg$parameter=='temp.anom' & pos.neg$value>0, 'species'])

# calculate mean and variance of parameter esimates
pbdf <- data.frame( t(postBeta$mean), taxon=colnames(postBeta$mean) )
names(pbdf) <- c("intercept","elev","elev2","year",
                 "elev:year","elev2:year","taxon")
## Add some basic trait information
pbdf$alga <- "alga"
pbdf$alga[c(3,15,20,56,57,62,68,82)] <- "invert"
# More specific groups


coef_plot <- pbdf %>% 
  select( -taxon ) %>%
  group_by( alga ) %>%
  gather( coefficient, value, -alga )

ggplot( coef_plot, aes(y=value, x=factor(1), col=alga)) +
  facet_wrap(~coefficient, scales="free_y") +
  geom_hline(yintercept=0) +
  stat_summary( fun.data=mean_cl_boot, geom='errorbar',
                size=1, width=0.1, position=position_dodge(width=0.9) )  +
  stat_summary( fun.y=mean, geom='point',
                size=3, position=position_dodge(width=0.9) )  +
  ylab('coefficient') + xlab('') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



## variance partitioning
VP = computeVariancePartitioning(m) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
plotVariancePartitioning(m, VP = VP)

VP.df <- as.data.frame(VP$vals) %>% 
  mutate(effect = factor(c("elevation","elev.square","year",
                           "elev:year","elev2:year",
                           "site","transect"), 
                         levels = rev(c("elevation","elev.square","year",
                                        "elev:year","elev2:year",
                                        "site","transect")), 
                         ordered = TRUE)) %>% 
  # mutate(effect = factor(c("elevation","elev.square","temp.anom.sum", "temp.anom.win",
  #                          "elev:year","elev2:year",
  #                          "site","transect","ty"), 
  #                        levels = rev(c("elevation","elev.square","temp.anom.sum", "temp.anom.win",
  #                                       "elev:year","elev2:year",
  #                                       "site","transect","ty")), 
  #                        ordered = TRUE)) %>% 
  # mutate(effect = factor(c("elevation","elev.square","temp.anom.sum", "temp.anom.win","transect","year","site"), 
  #                        levels = rev(c("elevation","elev.square","temp.anom.sum", "temp.anom.win","transect","year","site")), 
  #                        ordered = TRUE)) %>% 
  gather(key = species, value = variance, -effect) %>% 
  group_by(species) %>% 
  mutate(tempR2 = variance[effect == "year"])
  mutate(tempR2 = variance[effect == "temp.anom.sum"])

hold <- VP.df %>% filter(effect == "temp.anomaly") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(m$Y)[order(colSums(m$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(m$Y))

windows(8,5)
ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon", viridis(7)), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")





## associations
OmegaCor = computeAssociations(m)
supportLevel = 0.95
# choose the random variable to plot
rlevel = 3
toPlot = ((OmegaCor[[rlevel]]$support>supportLevel) 
          + (OmegaCor[[rlevel]]$support<(1-supportLevel))>0)*OmegaCor[[rlevel]]$mean
# reorder species matrix
plotorder <- order( postBeta$mean[5,], decreasing = TRUE )
toPlot <- toPlot[ plotorder, plotorder]
# reorder automatically
library(lessR)
mynewcor <- corReorder( toPlot, order="hclust", nclusters=4 )
# windows(12,12)
corrplot( mynewcor, method = "color", 
         col = colorRampPalette(c("blue","white","red"))(200),
           title = paste("random effect level:", m$rLNames[rlevel]), mar=c(0,0,0.5,0), tl.cex=0.6 )



# Plot predictions
# response curves####
# newXData <- data.frame(shore.height = seq(60,380,by = 1), 
#                        anom.pine.sum.1 = 1.1706,
#                        anom.pine.win = 0.8530 )
# 
# newDesign <- data.frame( site="new unit", year=c("2018","2019"), ty=rep("new unit")) 
newDF <- expand.grid( shore.height = seq(60,380,by = 1), 
             # anom.pine.sum.1 = 1.1706,
             # anom.pine.win = 0.8530,
             year=c(2012:2019),
             site="new unit", 
             transect="new unit" )
newXData   <- newDF[,1:2]
newDesign  <- newDF[,3:4]

# newDesign <- data.frame( site="new unit", year="2020", ty=rep("new unit", 321)) 
# newDesign <- data.frame( site="West Beach", 
#                          transect=rep("new unit", 321),
#                          year=2019) 
# newDesign <- purrr::map_df( newDesign,as.factor)
  
# predY <- predict(m, XData = newXData,
#                  studyDesign= newDesign,
#                  ranLevels = list(site=rL_site,transect=rL), expected = TRUE) #, transect=rL, year=rL_year
# 
# tmp = abind::abind(predY, along = 3)
# qpred = apply(tmp, c(1, 2), quantile, prob = 0.5, na.rm = TRUE)
# predictions <- bind_cols(data.frame(qpred), newXData)
# 
# qpred = apply(tmp, c(1, 2), quantile, prob = 0.025, na.rm = TRUE)
# predictions_low <- bind_cols(data.frame(qpred), newXData)
# 
# qpred = apply(tmp, c(1, 2), quantile, prob = 0.975, na.rm = TRUE)
# predictions_high <- bind_cols(data.frame(qpred), newXData)
# 
# ggplot(predictions, aes(x = shore.height, y = Fucus.distichus, col=year))+
#   geom_point() + scale_color_gradient(low = "slateblue1", high="slateblue4")
# ggplot(predictions, aes(x = shore.height, y = Barnacles, col=year))+
#   geom_point() + scale_color_gradient(low = "slateblue1", high="slateblue4")
# ggplot(predictions, aes(x = shore.height, y = Pyropia, col=year))+
#   geom_point() + scale_color_gradient(low = "slateblue1", high="slateblue4")
# ggplot(predictions, aes(x = shore.height, y = Nemalion.helminthoides, col=year))+
#   geom_point() + scale_color_gradient(low = "slateblue1", high="slateblue4")
# ggplot(predictions, aes(x = shore.height, y = Blidingia, col=year))+
#   geom_point() + scale_color_gradient(low = "slateblue1", high="slateblue4")
# 
# predictions_abund <- predictions %>% 
#   gather(key = taxon, value = N, Fucus.distichus:Colpomenia.peregrina)
# predictions_abund_low <- predictions_low %>% 
#   gather(key = taxon, value = N_low, Fucus.distichus:Colpomenia.peregrina)
# predictions_abund_high <- predictions_high %>% 
#   gather(key = taxon, value = N_high, Fucus.distichus:Colpomenia.peregrina)
# 
# predictions_abund <- left_join(left_join(predictions_abund, predictions_abund_low), predictions_abund_high)
# 
# predictions_abund$taxon <- factor(predictions_abund$taxon, levels = colnames(Y)[order(colSums(Y),decreasing = T)], ordered = FALSE)
# # predictions_abund$dispersal <- factor(predictions_abund$dispersal, levels = c("none", "low", "high"), ordered = TRUE)
# # predictions_abund$year <- rep( levels(newDesign$year), each=nrow(predictions_abund)/length(levels(newDesign$year))/97 )
# # save this
# save( predictions_abund, file = paste("R Code and Analysis/output from r/hmsc_pred",list.files( ModelDir )[model], sep="_") )

load( file=paste("R Code and Analysis/output from r/hmsc_pred",list.files( ModelDir )[model], sep="_") )

# use actual data
comm_select <- metacomm %>% 
  gather(key = taxon, value = N, Acrochaetium.sp.:Unknown.crust) %>%
  select( UID, shore.height=Shore_height_cm, taxon, N ) %>%
  filter( N>0 )
comm_select <- comm_select %>% 
  separate( UID,c("beach","b","zone","year","quad") ) %>% 
  unite( site, beach, b )
# remove rare taxa
comm_rare <- comm_select %>% 
  group_by( taxon ) %>% 
  summarize( N=sum(N) ) %>% 
  filter( N<=10 )
# merge
comm_final <- comm_select %>% 
  filter( !(taxon %in% comm_rare$taxon) )

# reorder levels of taxon
comm_final$taxon <- factor(comm_final$taxon, levels = colnames(Y)[order(colSums(Y),decreasing = T)], ordered = TRUE)

inverts <- c("Barnacles","Mytilus.sp.","Anemone","Bryozoan","Tunicata/Porifera","Pollicipes.polymerus","Tube.worms","Hydroid" )
comm_final$alga <- "alga"
comm_final$alga[ comm_final$taxon %in% inverts ] <- "invert"
predictions_abund$alga <- "alga"
predictions_abund$alga[ predictions_abund$taxon %in% inverts ] <- "invert"
# merge with other traits

# occurrences across years
commtab <- with(comm_select, table(taxon, year))
commtabdf <- as.data.frame( commtab )
commtab_wide <-commtabdf %>% 
  spread( taxon, Freq )
# commtab <- commtab[ order(rowSums(commtab),decreasing = T), ]
write.table( commtab, "R Code and Analysis/output from r/occurrence_table.txt")


## add trait information
taxa  <- read_csv( "Data/taxa/TaxonList_corrected_lumped_unique.csv" )
trait <- read_csv( "Data/taxa/Algae_functional_groups.csv" )

tt <- left_join( taxa, trait, by = c("taxon_revised"="taxon") )
anti_join( trait, taxa, by = c("taxon"="taxon_revised") )
tt$kelp_fucoid_turf[ tt$`non-alga flag` == "Animal" ]  <- "Animal"
tt$littler_groups[ tt$`non-alga flag` == "Animal" ]  <- "Animal"
tt$new_cat[ tt$`non-alga flag` == "Animal" ]  <- "Animal"
tt$new_cat_simple3[ tt$`non-alga flag` == "Animal" ]  <- "Animal"
tt$kelp_fucoid_turf[tt$taxon_lumped2=="Unknown CCA"] <- "crust"
tt$littler_groups[tt$taxon_lumped2=="Unknown CCA"] <- "crustose"
tt$kelp_fucoid_turf[tt$taxon_revised=="Ectocarpus sp."] <- "filament"
tt$littler_groups[tt$taxon_revised=="Ectocarpus sp."] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Melobesia sp."] <- "crustose coralline"
tt$littler_groups[tt$taxon_revised=="Melobesia sp."] <- "coralline"
tt$kelp_fucoid_turf[tt$taxon_revised=="Savoiea robusta"] <- "filament_turf"
tt$littler_groups[tt$taxon_revised=="Savoiea robusta"] <- "filament"
tt$new_cat_simple3[tt$taxon_revised=="Savoiea robusta"] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Symphyocladia plumosa"] <- "filament_turf"
tt$littler_groups[tt$taxon_revised=="Symphyocladia plumosa"] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Hedophyllum nigripes"] <- "Kelp"
tt$new_cat_simple3[tt$taxon_revised=="Hedophyllum nigripes"] <- "Kelp"
tt$littler_groups[tt$taxon_revised=="Hedophyllum nigripes"] <- "thick leathery"
tt$kelp_fucoid_turf[tt$taxon_revised=="Hedophyllum sessile"] <- "Kelp"
tt$new_cat_simple3[tt$taxon_revised=="Hedophyllum sessile"] <- "Kelp"
tt$littler_groups[tt$taxon_revised=="Hedophyllum sessile"] <- "thick leathery"
tt$kelp_fucoid_turf[tt$taxon_revised=="Chamberlainium tumidum"] <- "crustose coralline"
tt$new_cat_simple3[tt$taxon_revised=="Chamberlainium tumidum"] <- "crustose coralline"
tt$littler_groups[tt$taxon_revised=="Chamberlainium tumidum"] <- "coralline"
tt$kelp_fucoid_turf[tt$taxon_revised=="Ulothrix-Urospora sp."] <- "green turf"
tt$new_cat_simple3[tt$taxon_revised=="Ulothrix-Urospora sp."] <- "filament"
tt$littler_groups[tt$taxon_revised=="Ulothrix-Urospora sp."] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Prionitis lanceoloata"] <- "red turf"
tt$littler_groups[tt$taxon_revised=="Prionitis lanceoloata"] <- "finely branching"
tt$kelp_fucoid_turf[tt$taxon_revised=="Derbesia marina"] <- "green turf"
tt$littler_groups[tt$taxon_revised=="Derbesia marina"] <- "filament"
tt$new_cat_simple3[tt$taxon_revised=="Derbesia marina"] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Acrosiphonia spp."] <- "filament_turf"
tt$littler_groups[tt$taxon_revised=="Acrosiphonia spp."] <- "filament"
tt$new_cat_simple3[tt$taxon_revised=="Ectocarpus sp."] <- "filament"
tt$new_cat_simple3[tt$taxon_revised=="Melobesia sp."] <- "crustose coralline"
tt$new_cat_simple3[tt$taxon_revised=="Symphyocladia plumosa"] <- "filament"
tt$new_cat_simple3[tt$taxon_revised=="Chamerlaini"] <- "filament"



tt <- tt[ !is.na(tt$kelp_fucoid_turf),]
tt$taxon <- gsub( " ",".",tt$taxon_lumped2 )
tt2 <- tt %>% 
  select( taxon, funct=kelp_fucoid_turf, littler_groups, new_cat, 
          new_cat_simple1, new_cat_simple2, new_cat_simple3 ) %>% 
  distinct()
# tt2$funct[ tt2$funct=="kelp"] <- "Kelp"
# tt2$funct[ tt2$funct=="green turf"] <- "filament_turf"

tt2 %>% 
  select(taxon, )
questions <- c(4, 20, 76, 213, 215, 219)
ttuse <- tt2[ -c(4, 213, 215, 219), ]
ttuse$new_cat_simple3[ttuse$taxon=="Unknown.CCA"] <- "crustose coralline"
ttuse$new_cat_simple3[ttuse$new_cat_simple3=="Animal"] <- "animal"
ttuse$new_cat_simple3[ttuse$new_cat_simple3=="Kelp"] <- "kelp"

ttuse %>% 
  select(taxon) %>% 
  duplicated()

# colors for functional groups
unique(ttuse$new_cat_simple3)
ttuse$fill <- "whitesmoke"
ttuse$fill[ttuse$new_cat_simple3=="animal"] <- "deepskyblue"


predictions_abund_trait <- left_join( predictions_abund, ttuse)

predictions_abund_trait$funct[predictions_abund_trait$taxon=="Tunicata.Porifera"] <- "Animal"

predictions_use <- predictions_abund_trait %>% 
  filter( !is.na(new_cat_simple3) )




# only pull the 10 most common taxa
pick4 <- colnames(Y)[c(27,35,48,76)] 
top6 <- colnames(Y)[c(1:4,6,14,15,29, 36)]
# 10 rarest taxa
bot10 <- colnames(Y)[(length(colnames(Y))-9):length(colnames(Y))]

windows(6,4)
ggplot( filter(predictions_abund,taxon %in% pick4, year %in% c(2012,2019) ), aes(x = shore.height, y = N,
                                                        fill=factor(year) ))+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, col="gray50")+
  geom_line(size = 1)+
  facet_wrap(~taxon, scales = "free_y")+
  theme_classic() +#+
  # scale_fill_viridis_d() +
  # scale_color_manual(values=c("black","black"))+
  scale_fill_manual(values=c("whitesmoke", "mediumslateblue"))+
  geom_point( data = filter( comm_final, taxon %in% pick4, year %in% c(2012,2019)), pch=21 ) +
  scale_y_sqrt() + ylab("Percent cover") + xlab("Shore height (cm)") +
  labs( fill="Year" )
ggsave("R Code and Analysis/Figs/hmsc_response_curves_WSN.pdf", width = 5, height = 4)




# Find the peak for each instance
peaks <- predictions_use %>% 
  group_by( year, taxon, new_cat_simple3, fill ) %>% 
  summarize( peak = shore.height[which(N==max(N))] )
ggplot( filter(peaks, taxon %in% pick4), aes(x=year,y=peak) ) + facet_wrap(~taxon) +
  geom_point()

# get difference between peaks for 2012 and 2019
peak_shift <- peaks %>% 
  group_by( taxon, new_cat_simple3, fill ) %>% 
  filter( year %in% c(2012,2019) ) %>% 
  summarize( shift = diff(peak))
peak_shift %>% arrange(-shift)

# plot by functional group
windows(6,4)
ggplot( peak_shift, aes(x=reorder(new_cat_simple3, shift, FUN = median), y=shift/100,
                        fill=fill) ) + 
  geom_hline (yintercept=0, lty=2 ) +
  geom_boxplot()  + geom_point() +
  xlab("Functional group") + ylab("Peak shift (meters)") +
  scale_fill_manual(values=rev(unique(peak_shift$fill))) +
  theme_classic() +
  theme( axis.text.x = element_text(angle=45,hjust=1,vjust=1) ) +
  theme(legend.position = "none")
ggsave("R Code and Analysis/Figs/hmsc_peak_shift.pdf", width = 6, height = 4)



 
summary( lm(shift~1,peak_shift) )
summary( lm(shift~new_cat_simple3,peak_shift) )
summary( lme4::lmer(shift~1+(1|new_cat_simple3),peak_shift) )

# find the integral of the function for each year
integ <- predictions_use %>% 
  group_by( year, taxon, new_cat_simple3, fill ) %>% 
  summarize( integral = sum(N) )
ggplot( filter(integ, taxon %in% top6), aes(x=year,y=integral) ) + facet_wrap(~taxon) +
  geom_point()
# get difference between integrand for 2012 and 2019
int_shift <- integ %>% 
  group_by( taxon, new_cat_simple3, fill ) %>% 
  filter( year %in% c(2012,2019) ) %>% 
  summarize( shift = integral[2]/integral[1] )

int_shift %>% arrange(-shift)
int_shift %>% arrange(shift)
filter(int_shift,new_cat_simple3=="animal") %>% arrange(shift)

summary( lm(shift~1,int_shift) )
summary( lm(shift~1,filter(int_shift,shift<16)) )

mod <- lm(log(shift,base=2)~1, int_shift )
mod <- lm(log(shift,base=2)~1, filter(int_shift,shift<60) )
summary(mod) 

# plot by functional group
ggplot( filter(int_shift,shift<60), aes(x=reorder(new_cat_simple3, shift, FUN = median), 
                       y=log(shift,base=2), fill=fill ) ) + 
  geom_hline (yintercept=0, lty=2 ) +
  geom_boxplot()  + geom_point() +
  xlab("Functional group") + ylab("Abundance shift") +
  theme_classic() +
  scale_y_continuous( breaks=c(6,4,2,0,-2,-4), 
                      labels=c('64x','16x','4x','0','1/4x','1/16x')) +
  scale_fill_manual(values=rev(unique(peak_shift$fill))) +
  theme( axis.text.x = element_text(angle=45,hjust=1,vjust=1) ) +
  theme(legend.position = "none")
ggsave("R Code and Analysis/Figs/hmsc_int_shift.pdf", width = 6, height = 4)

int_shift %>% arrange(shift) 



## How to solve the peak and integral problem a little more elegantly?
p <- postBeta$mean[, "Alaria.marginata"]
elev <- 200
y1 <- 2012
p[1] + y1*p[5]*p[2]*elev + y1*p[6]+p[3]*I(elev^2) + y1*p[4] 


























###
####
#####
# Patrick Thompson's Code Below



ggplot( data = filter( comm_final, taxon %in% top10, year %in% levels(newDesign$year)), aes(x = shore.height, y = N, group=taxon, col=year ) ) +
  facet_wrap(~taxon, scales = "free_y")+
  theme_classic() +#+
  # scale_color_viridis_d(end = 0.8)+
  # scale_fill_viridis_d(end = 0.8)+
  geom_point(  ) +
  # geom_smooth( method="glm", aes(group=year), se=F, method.args=list(family="quasipoisson"), formula=y~I(x^2) )+
  scale_y_sqrt()

filter(predictions_abund, week == 32, dispersal == "none", species != "simocephalus_serrulatus") %>% 
  group_by(species) %>% 
  mutate(max_temp = temp[N == max(N)]) %>%
  filter(temp<=30, temp >=18) %>% 
  mutate(N_norm = decostand(N, method = "range")) %>%
  ggplot(aes(x = temp, y = N_norm, group = species, color = fct_reorder(species, max_temp, .desc = TRUE)))+
  #geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  ylab("Normalized predicted abundance")+
  scale_color_brewer(palette = "RdYlBu", name = "")+
  theme_classic()+
  xlab("temperature Â°C")
ggsave("./figures/HMSC/response_curves_normalized.pdf", width = 6, height = 5)

filter(predictions_abund, week == 32, dispersal == "none", species != "simocephalus_serrulatus") %>% 
  group_by(species) %>% 
  mutate(max_temp = temp[N == max(N)]) %>%
  filter(temp<=30, temp >=18) %>% 
  ggplot(aes(x = temp, y = N, group = species, color = fct_reorder(species, max_temp, .desc = TRUE)))+
  #geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  ylab("Predicted abundance")+
  scale_color_brewer(palette = "RdYlBu", name = "")+
  theme_classic()+
  xlab("temperature Â°C")+
  scale_y_sqrt()
ggsave("./figures/HMSC/response_curves_base.pdf", width = 6, height = 5)
