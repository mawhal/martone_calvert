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
library( cowplot )
library( ggrepel )

## useful references for this code
citation( "Hmsc" )
# https://github.com/hmsc-r/HMSC

## directories
ModelDir = paste0( here::here(), "/R Code and Analysis/models" )
MixingDir = paste0( here::here(), "/R Code and Analysis/mixing")



## load the model
list.files( ModelDir )
model = 13
mload <- load( paste(ModelDir,list.files( ModelDir ), sep="/")[model] )

# load the data
metacomm <- read_csv("R Code and Analysis/output from r/community.csv" )

##
## MCMC convergence

# We evaluate MCMC convergence in terms of four kinds of parameters that we are especially interested in: the species niches `Beta`, the influence of traits on species niches `Gamma`, residual species associations `Omega`, and the strength of phylogenetic signal `rho`. As there are `r ns` species, the matrix `Omega` is a `r ns` x `r ns` matrix with `r ns*ns` elements. To avoid excessive computational times, we evaluate MCMC convergence for `Omega` only for a subsample of 100 randomly selected species pairs.
mpost = convertToCodaObject(m)
# plot traces for particular species 
taxa <- rep(colnames(m$Y), each=12)
toplot <- which(taxa %in% "Hedophyllum.sessile")
toplot <- which(taxa %in% "Polysiphonia")
toplot <- which(taxa %in% "Fucus.distichus")
betaplot <- lapply(mpost$Beta,function(z) z[,toplot])
lapply(betaplot, colnames)
# betaplot <- lapply( betaplot, function(z) {
#   colnames(z) <- c('int','elev','elev2','year','elev:year','elev2:year')
# })
params <- c('int','elev','elev2','year','elev:year','elev2:year')
betaplot <- as.mcmc.list(betaplot)
par( mfrow=c(6,4),mar=c(2,2,0,0)+0.1)
plot(betaplot, auto.layout = FALSE, cex.main=0.8, cex.sub=0.8)
# plot(mpost$Beta)
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
# windows(5,5)
lapply( MF, summary)
par( mfrow=c(3,3), mar=c(3,3,1,1)+0.1 )
lapply( MF, boxplot )



## parameter estimates
postBeta = getPostEstimate(m, parName = "Beta")
# windows(5,8)
# plotBeta(m, post = postBeta, param = "Sign", supportLevel = 0.95, mar=c(7,11,0,6))
postBeta$mean[, c("Alaria.marginata","Hedophyllum.sessile","Polysiphonia")]

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

windows(6,4)
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
# plotVariancePartitioning(m, VP = VP)

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





#####
# Model  predictions
# response curves
newXData <- data.frame(shore.height = seq(60,380,by = 1),
                       anom.pine.sum.1 = 1.1706,
                       anom.pine.win = 0.8530 )
newXData <- data.frame(shore.height = seq(60,380,by = 1),
                       anom.pine.sum.1 = 1.1706,
                       anom.pine.win = 0.8530 )
years <- as.data.frame( m$XData %>% select(year1,year2) %>% mutate(year1=round(year1,7),year2=round(year2,7)) %>% distinct() )
elevs <- as.data.frame( m$XData %>% select(elev1,elev2) %>% mutate(elev1=round(elev1,5),elev2=round(elev2,5)) %>% distinct() %>% arrange(elev1) )
plot(elevs)
plot(years)

newDF <- expand.grid( shore.height = seq(60,380,by = 1),
             year=c(2012:2019),
             site="new unit",
             transect="new unit" )
newDF <- data.frame( merge(years, elevs),
             site="new unit",
             transect="new unit" )
newXData   <- newDF[,1:(ncol(newDF)-2)]
newDesign  <- newDF[,(ncol(newDF)-1):ncol(newDF)]

predY <- predict(m, XData = newXData,
                 studyDesign= newDesign,
                 ranLevels = list(site=rL_site,transect=rL), expected = TRUE) #, transect=rL, year=rL_year

tmp = abind::abind(predY, along = 3)
qpred = apply(tmp, c(1, 2), quantile, prob = 0.5, na.rm = TRUE)
predictions <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.025, na.rm = TRUE)
predictions_low <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.975, na.rm = TRUE)
predictions_high <- bind_cols(data.frame(qpred), newXData)

ggplot(predictions, aes(x = shore.height, y = Fucus.distichus, col=year))+
  geom_point() + scale_color_gradient(low = "slateblue1", high="slateblue4")
ggplot(predictions, aes(x = shore.height, y = Fucus.distichus, col=year))+
  geom_point() + scale_color_gradient(low = "slateblue1", high="slateblue4")
#

predictions_abund <- predictions %>%
  gather(key = taxon, value = N, Fucus.distichus:Mazzaella.parvula)
predictions_abund_low <- predictions_low %>%
  gather(key = taxon, value = N_low, Fucus.distichus:Mazzaella.parvula)
predictions_abund_high <- predictions_high %>%
  gather(key = taxon, value = N_high, Fucus.distichus:Mazzaella.parvula)

predictions_abund <- left_join(left_join(predictions_abund, predictions_abund_low), predictions_abund_high)

predictions_abund$taxon <- factor(predictions_abund$taxon, levels = colnames(Y)[order(colSums(Y),decreasing = T)], ordered = FALSE)
# predictions_abund$dispersal <- factor(predictions_abund$dispersal, levels = c("none", "low", "high"), ordered = TRUE)
# predictions_abund$year <- rep( levels(newDesign$year), each=nrow(predictions_abund)/length(levels(newDesign$year))/97 )
# save this
# save( predictions_abund, file = paste("R Code and Analysis/output from r/hmsc_pred",list.files( ModelDir )[model], sep="_") )

# load( file = paste("R Code and Analysis/output from r/hmsc_pred",list.files( ModelDir )[model], sep="_") )
# load( file = paste("R Code and Analysis/output from r/hmsc_pred_model_5_chains_4_thin_100_samples_1000.Rdata") )

# use actual data
comm_select <- metacomm %>% 
  gather(key = taxon, value = N, Acrochaetium.sp.:Unknown.crust) %>%
  select( UID, shore.height=Shore_height_cm, taxon, N )# %>%
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
comm_final$taxon <- factor(comm_final$taxon, levels = colnames(m$Y)[order(colSums(m$Y),decreasing = T)], ordered = TRUE)

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
# # commtab <- commtab[ order(rowSums(commtab),decreasing = T), ]
# write.table( commtab, "R Code and Analysis/output from r/occurrence_table.txt")


## add trait information
taxa  <- read_csv( "Data/taxa/TaxonList_corrected_lumped_unique.csv" )
trait <- read_csv( "Data/taxa/Algae_functional_groups.csv" )

tt <- left_join( taxa, trait, by = c("taxon_revised"="taxon") )
anti_join( trait, taxa, by = c("taxon"="taxon_revised") )
tt$kelp_fucoid_turf[ tt$non.alga.flag == "Animal" ]  <- "Animal"
tt$littler_groups[ tt$non.alga.flag == "Animal" ]  <- "Animal"
tt$kelp_fucoid_turf[tt$taxon_lumped2=="Unknown CCA"] <- "crust"
tt$littler_groups[tt$taxon_lumped2=="Unknown CCA"] <- "crustose"
tt$kelp_fucoid_turf[tt$taxon_revised=="Ectocarpus sp."] <- "filament"
tt$littler_groups[tt$taxon_revised=="Ectocarpus sp."] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Melobesia sp."] <- "crustose coralline"
tt$littler_groups[tt$taxon_revised=="Melobesia sp."] <- "coralline"
tt$kelp_fucoid_turf[tt$taxon_revised=="Savoiea robusta"] <- "filament_turf"
tt$littler_groups[tt$taxon_revised=="Savoiea robusta"] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Symphyocladia plumosa"] <- "filament_turf"
tt$littler_groups[tt$taxon_revised=="Symphyocladia plumosa"] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Hedophyllum nigripes"] <- "Kelp"
tt$littler_groups[tt$taxon_revised=="Hedophyllum nigripes"] <- "thick leathery"
tt$kelp_fucoid_turf[tt$taxon_revised=="Hedophyllum sessile"] <- "Kelp"
tt$littler_groups[tt$taxon_revised=="Hedophyllum sessile"] <- "thick leathery"
tt$kelp_fucoid_turf[tt$taxon_revised=="Chamberlainium tumidum"] <- "crustose coralline"
tt$littler_groups[tt$taxon_revised=="Chamberlainium tumidum"] <- "coralline"
tt$kelp_fucoid_turf[tt$taxon_revised=="Ulothrix-Urospora sp."] <- "green turf"
tt$littler_groups[tt$taxon_revised=="Ulothrix-Urospora sp."] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Prionitis lanceoloata"] <- "red turf"
tt$littler_groups[tt$taxon_revised=="Prionitis lanceoloata"] <- "finely branching"
tt$kelp_fucoid_turf[tt$taxon_revised=="Derbesia marina"] <- "green turf"
tt$littler_groups[tt$taxon_revised=="Derbesia marina"] <- "filament"
tt$kelp_fucoid_turf[tt$taxon_revised=="Acrosiphonia spp."] <- "filament turf"
tt$littler_groups[tt$taxon_revised=="Acrosiphonia spp."] <- "filament"



tt <- tt[ !is.na(tt$kelp_fucoid_turf),]
tt$taxon <- gsub( " ",".",tt$taxon_lumped3 )
tt <- tt %>% 
  select( taxon, funct=kelp_fucoid_turf, littler_groups ) %>% 
  distinct()
tt$funct[ tt$funct=="kelp"] <- "Kelp"
tt$funct[ tt$funct=="green turf"] <- "filament_turf"

tt[ tt$funct == "crustose coralline",]

ttuse <- tt
# tt <- tt[-c(4,21,213),]
ttuse$fill <- "whitesmoke"
ttuse$fill[ttuse$funct=="Animal"] <- "deepskyblue"

predictions_abund_trait <- left_join( predictions_abund, ttuse)
sort(unique(predictions_abund_trait$taxon))
sort(unique(predictions_abund_trait$taxon[ predictions_abund_trait$funct== "crustose coralline"]))
predictions_abund_trait$funct[predictions_abund_trait$taxon=="Tunicata.Porifera"] <- "Animal"
predictions_abund_trait$funct[predictions_abund_trait$taxon=="coralline.crust"] <- "crustose coralline"
predictions_abund_trait$fill[predictions_abund_trait$taxon=="coralline.crust"] <- "whitesmoke"







# only pull the 10 most common taxa
top6 <- colnames(m$Y)[c(1:6,9,14,15)]
upX <- colnames(m$Y)[c(8,15,18,23,27,29,35,37,48)]
downX <- colnames(m$Y)[c(43,36,5,24,28,33,30)]
upY <- colnames(m$Y)[c(15,42,39,47,44,46,20,22,9)]
customXY <- colnames(m$Y)[c(1,3,4,5,6,8,9,14,15,27,30,37,42)]
# 10 rarest taxa
bot10 <- colnames(m$Y)[(length(colnames(m$Y))-8):length(colnames(m$Y))]

taxa2plot <- top6

# windows(6,4)
ggplot( filter(predictions_abund,taxon %in% taxa2plot & year %in% c(2012,2019)), aes(x = shore.height, y = N,
                                                        fill=factor(year) ))+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, col="gray75")+
  facet_wrap(~taxon, scales = "free_y")+
  theme_classic() +#+
  # scale_fill_viridis_d() +
  # scale_color_manual(values=c("black","black"))+
  scale_fill_manual(values=c("whitesmoke", "mediumslateblue"))+
  geom_point( data = filter( comm_final, taxon %in% taxa2plot, year %in% c(2012,2019)), pch=21 ) +
  geom_line(size = 0.5)+
  scale_y_sqrt() + ylab("Percent cover") + xlab("Shore height (cm)")
# ggsave("R Code and Analysis/Figs/hmsc_response_curves.pdf", width = 6, height = 4)

# just show Fucus
fuc <- "Fucus.distichus"
fuc <- "Polysiphonia"
fuc <- "Mazzaella.parvula"
fuc <- "Microcladia.borealis"
fuc <- "Palmaria.hecatensis"
fuc <- "Farlowia.mollis"
fuc <- "Cladophora.columbiana"
fuc <- "coralline.crust"
fuc <- "Egregia.menziesii"
fuc <- "Hedophyllum.sessile"
fuc <- "Crusticorallina.muricata"
fuc <- "Lithothamnion.phymatodeum"
fuc <- "Mytilus.sp."
fuc <- "Barnacles"
fuc <- "Alaria.marginata"
fuc <- "Elachista.fucicola"
# windows(5,4)
ggplot( filter(predictions_abund,taxon %in% fuc ), 
        aes(x = shore.height, y = N ))+
  facet_wrap(~year)+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, fill="gray75")+
  geom_line(size = 0.5)+
  theme_classic() +
  geom_point( data = filter( comm_final, taxon %in% fuc), pch=21 ) +
  scale_y_sqrt(breaks=c(1,10,50,100,200)) + ylab("Percent cover") + xlab("Shore height (cm)") +
  coord_cartesian(ylim = c(-0, 100)) 
#



# Find the predicted peak for each instance
peaks <- predictions_abund_trait %>% 
  group_by( year, taxon, funct, fill ) %>% 
  summarize( peak = mean(shore.height[which(N==max(N))]) )
peaks$peak
ggplot( filter(peaks, taxon %in% top6), aes(x=year,y=peak) ) + facet_wrap(~taxon) +
  geom_point()

# get difference between peaks for 2012 and 2019
peak_shift <- peaks %>% 
  group_by( taxon, funct, fill ) %>% 
  filter( year %in% c(2012,2019) ) %>% 
  summarize( shift = diff(peak))

# merge 2012 peaks with peak shift to compare shift relative to starting point
peak_initial <- peaks %>% filter( year==2012 )
peak_compare <- left_join( peak_shift, peak_initial )

cutoff <- 20
peak_compare %>% filter(shift <= -cutoff)
peak_compare %>% filter(shift >= cutoff)
peak_compare %>% filter(shift < cutoff & shift > -cutoff)
peak_compare %>% filter(shift == 0 )


# plot by functional group
ggplot( peak_shift, aes(x=funct, y=shift) ) + 
  geom_hline (yintercept=0, lty=2 ) +
  geom_boxplot(fill="whitesmoke")  + geom_point() +
  xlab("Functional group") + ylab("Peak shift (cm)") +
  theme_classic() +
  theme( axis.text.x = element_text(angle=45,hjust=1) )
trait.p <- ggplot( peak_shift, aes(x=reorder(funct, shift, FUN = median), y=shift/100,
                        fill=fill) ) + 
  geom_hline (yintercept=0, lty=2 ) +
  geom_boxplot()  + geom_point() +
  xlab("Functional group") + ylab("Peak shift (meters)") +
  scale_fill_manual(values=rev(unique(peak_shift$fill))) +
  theme_classic() +
  theme( axis.text.x = element_text(angle=45,hjust=1,vjust=1) ) +
  theme(legend.position = "none")

peak_shift %>% arrange(-shift)
peak_shift %>% arrange(shift)
peak_shift %>% filter(shift >-10 & shift < 10)
shift_increase <- peak_shift %>% arrange(-shift)
choose <- shift_increase$taxon[1:6]
ggplot( filter(peaks, taxon %in% choose), aes(x=year,y=peak) ) + facet_wrap(~taxon) +
  geom_point()

ggplot( filter(predictions_abund,taxon %in% "Alaria.marginata" ), aes(x = shore.height, y = N,
                                                                                 fill=factor(year), col=factor(year) ))+
  # geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, col="gray50")+
  geom_line(size = 1)+
  facet_wrap(~taxon, scales = "free_y")+
  theme_classic() +#+
  # scale_fill_viridis_d() +
  # scale_color_manual(values=c("black","black"))+
  # scale_fill_manual(values=c("whitesmoke", "mediumslateblue"))+
  geom_point( data = filter( comm_final, taxon %in% "Alaria.marginata"), pch=21 ) +
  scale_y_sqrt() + ylab("Percent cover") + xlab("Shore height (cm)") +
  labs( fill="Year" )

lm_shift <- lm( shift~1, peak_shift)
summary(lm_shift)



# need to add lines showing the realm of possible shifts
xs <- range( m$XData$shore.height )
y1 <- c(0,diff(xs))
y2 <- c(-diff(xs),0)
df.bound <- data.frame( x1=xs[1],x2=xs[2],y1,y2 )
df.poly <- data.frame( x=rep(xs,each=2), y=c(0,diff(xs),0,-diff(xs)) )



(a <- ggplot( peak_compare, aes(x=peak,y=shift)) +
    geom_polygon( data=df.poly, aes(x=x,y=y), fill='whitesmoke', col='slategray', lty=2) +
    geom_hline( yintercept = 0, lty=2 ) +
     geom_smooth(method='lm', se=F, col='black') +
    geom_point( size=3, pch=1, col='slateblue' ) +
  ylab("peak elevationshift (cm)") + xlab("initial peak elevation (cm)") +
  theme_classic() )

  
# find the integral of the function for each year
integ <- predictions_abund_trait %>% 
  group_by( year, taxon, funct, fill ) %>% 
  summarize( integral = sum(N) )
ggplot( filter(integ, taxon %in% top6), aes(x=year,y=integral) ) + facet_wrap(~taxon) +
  geom_point()
# get difference between integrand for 2012 and 2019
int_shift <- integ %>% 
  group_by( taxon, funct, fill ) %>% 
  filter( year %in% c(2012,2019) ) %>% 
  summarize( shift = integral[2]/integral[1] )

# plot by functional group
int_shift %>% arrange(-shift) 

trait.a <- ggplot( int_shift, aes(x=reorder(funct, shift, FUN = median), 
                                        y=log(shift,base=2), fill=fill ) ) + 
  geom_hline (yintercept=0, lty=2 ) +
  geom_boxplot()  + geom_point() +
  xlab("Functional group") + ylab("Abundance shift") +
  theme_classic() +
  scale_y_continuous( breaks=c(4,2,0,-2,-4), 
                      labels=c('16x','4x','0','1/4x','1/16x')) +
  scale_fill_manual(values=rev(unique(peak_shift$fill))) +
  theme( axis.text.x = element_text(angle=45,hjust=1,vjust=1) ) +
  theme(legend.position = "none")

cowplot::plot_grid( trait.p, trait.a, ncol=1, align='hv', labels="AUTO" )

summary( lm(shift~1,int_shift) )
mod <- lm(log(shift,base=2)~1, int_shift )
summary(mod) 

# by short height
int_initial <- integ %>% filter( year==2012 )
peak_initial <- peaks %>% filter( year==2012 )
peak_compare <- left_join( peak_shift, peak_initial )
abund_compare <- left_join( int_shift, int_initial )

# combined figure of mean elevation shift and abundance shift
compare_all <- left_join( peak_compare, abund_compare, by=c("taxon","funct","fill","year") )

# shift in abundance ~ initial peak elevation
(b <- ggplot( compare_all, aes(x=peak,y=log(shift.y,base=2))) + 
  geom_hline( yintercept=0, lty=2 ) +
  geom_point(size=3, pch=1, col='slateblue') + 
  ylab("abundance shift") + xlab("initial peak elevation (cm)") +
  scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
                      labels=c('10x','2x','0','1/2x','1/10x')) +
  theme_classic() )

# peak shift ~ initial abundance
(c <- ggplot( compare_all, aes(x=integral/30,y=shift.x)) + 
    geom_hline( yintercept=0, lty=2 ) +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("peak elevationshift (cm)") + xlab("appox. initial cover (%)") +
    # scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
    #                     labels=c('10x','2x','0','1/2x','1/10x')) +
    scale_x_continuous(trans = "log2") +
    theme_classic() )

# abund shift ~ initial abundance
summary(lm( log(shift.y,base=2)~log(integral,base=2), compare_all ))
(d <- ggplot( compare_all, aes(x=integral/30,y=log(shift.y,base=2))) + 
    geom_hline( yintercept=0, lty=2 ) +
    geom_smooth(method = 'lm', se = T, col='black' ) +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("Abundance shift") + xlab("appox. initial cover (%)") +
    scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
                        labels=c('10x','2x','0','1/2x','1/10x')) +
    scale_x_continuous(trans = "log2") +
    theme_classic() )

cowplot::plot_grid( a, b, c, d, ncol=2, align = 'hv', labels = "AUTO" )
ggsave(file="R Code and Analysis/Figs/abundance+peak_shift_intial.svg",width = 6, height = 6)

compare_all %>% arrange(shift.y)
compare_all %>% arrange(-shift.x)

# basically no relationship between peak and abundance shift
ggplot( compare_all, aes(x=shift.y,y=shift.x) ) + 
  geom_point() + geom_smooth( method='lm' )
ggplot( compare_all, aes(x=log(shift.y,base=2),y=shift.x) ) + 
  geom_point() + geom_smooth( method='lm' )

with( compare_all, cor.test( shift.x,shift.y) )
with( compare_all, cor.test( shift.x,log(shift.y,base=2)) )


# boxplots
ylimits1 <- c(-max(abs(range(compare_all$shift.x))),max(abs(range(compare_all$shift.x))))
ylimits2 <- c(-6,6) #c(2^-max(abs(log(range(compare_all$shift.y),base=2))), 2^max(abs(log(range(compare_all$shift.y),base=2))))
(bpa <- ggplot( compare_all, aes(y=shift.x,x=1) ) +
    # geom_boxplot(notch = T,outlier.color = NA) + 
    geom_violin(outlier.color = NA, draw_quantiles = 0.5, trim=FALSE ) +
    geom_hline(yintercept=0,lty=2)+
    geom_point(alpha=0.25) +
    ylab( "Peak elevation shift (cm)" )   +
    scale_y_continuous(limits=ylimits1,breaks=c(-300,-200,-100,-50,0,50,100,200,300),
                       position="left") +
    theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) )
(bpb <- ggplot( compare_all, aes(y=log(shift.y,base=2),x=1) ) +
    # geom_boxplot(notch=TRUE) + 
    geom_violin(outlier.color = NA, draw_quantiles = 0.5, trim=FALSE ) +
    geom_hline(yintercept=0,lty=2)+
    geom_point(alpha=0.25) +
    ylab( "Abundance shift" ) +
    scale_y_continuous(limits=ylimits2,breaks=c(log(50,base=2),log(10,base=2),log(2,base=2),0,log(0.5,base=2),log(1/10,base=2),log(1/50,base=2)),
                      labels=c('50x','10x','2x','0','1/2x','1/10x','1/50x'),
                      position="right") +
    theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) )
# windows(3,2.5)
cowplot::plot_grid(bpa,bpb, align = "h", axis='tblr', labels = "AUTO")
#

# add nice scatteplot
# taxa to plot
taxalabel <- top6
taxalabel <- customXY
xy <- ggplot( compare_all, aes(x=log(shift.y,base=2),y=shift.x)) + 
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(pch=21,size=2,fill="whitesmoke") +
  geom_text_repel(data=filter(compare_all,taxon%in%taxalabel),aes(label=taxon),
                  box.padding = 0.3, point.padding = 0.1, size=3) +
  geom_point( data=filter(compare_all,taxon%in%taxalabel)) +
  theme_classic() +
  scale_x_continuous(breaks=c(log(50,base=2),log(10,base=2),log(5,base=2),log(2,base=2),0,log(0.5,base=2),log(1/5,base=2),log(1/10,base=2),log(1/50,base=2)),
                     labels=c('50x','10x','5x','2x','0','1/2x','1/5x','1/10x','1/50x'),
                     position="bottom")+
  xlab("Abundance shift") + ylab("Elevation shift (cm)")

# a densities not violin plots
ydens <- axis_canvas(xy, axis = "y", coord_flip = TRUE)+
  geom_vline(xintercept=mean(compare_all$shift.x), col='red' ) +
  geom_vline(xintercept=0) +
  geom_density(data = compare_all, aes(x = shift.x),
               alpha = 0.7, size = 0.5) +
  coord_flip()
  
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
xdens <- axis_canvas(xy, axis = "x")+
  geom_vline(xintercept=mean(log(compare_all$shift.y,base=2)), col='red' ) +
  geom_vline(xintercept=0) +
  geom_density(data = compare_all, aes(x = log(shift.y,base=2) ),
               alpha = 0.7, size = 0.5)

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

p1 <- insert_xaxis_grob(xy, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
ggsave(file="R Code and Analysis/Figs/abundance~peak.svg",width = 4, height = 4)

write_csv( compare_all, "R Code and Analysis/output from r/shifts_predicted.csv")
#
  

#


































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
