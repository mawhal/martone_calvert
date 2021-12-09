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
library( colorspace )
# library( vioplot )
library( rstan )
library( scales )
library( lessR )
source( "R/mcmc.list2array.R")

## useful references for this code
# citation( "Hmsc" )
# https://github.com/hmsc-r/HMSC

## directories
ModelDir = paste0( here::here(), "/R/models" )
MixingDir = paste0( here::here(), "/R/mixing")



## load the model
list.files( ModelDir )
model = "model_elevxyear_hurdle_chains_4_thin_100_samples_250_quad.Rdata"
# model = "model_elevxyear_hurdle_chains_4_thin_100_samples_250.Rdata"
# model = "model_elevxyear_hurdle_test_chains_1_thin_1_samples_5.Rdata"
mload <- load( paste(ModelDir,model, sep="/") )

# load the data
metacomm <- read_csv("R/output/community_all.csv" )
comm <- read_csv("R/output/community.csv")
##
## MCMC convergence

# We evaluate MCMC convergence in terms of four kinds of parameters that we are especially interested in: the species niches `Beta`, the influence of traits on species niches `Gamma`, residual species associations `Omega`, and the strength of phylogenetic signal `rho`. As there are `r ns` species, the matrix `Omega` is a `r ns` x `r ns` matrix with `r ns*ns` elements. To avoid excessive computational times, we evaluate MCMC convergence for `Omega` only for a subsample of 100 randomly selected species pairs.
mpost_pa   = convertToCodaObject(models[[1]])
mpost_abun = convertToCodaObject(models[[2]])

# trace plots #####
# plot traces for particular species
taxa <- rep(colnames(models[[1]]$Y), each=7)
toplot <- which(taxa %in% "Hedophyllum.sessile")
# toplot <- which(taxa %in% "Mazzaella.parvula")
betaplot <- lapply(mpost_pa$Beta,function(z) z[,toplot])
betaplot2 <- lapply(mpost_abun$Beta,function(z) z[,toplot])
lapply(betaplot, colnames)
# betaplot <- lapply( betaplot, function(z) {
#   colnames(z) <- c('int','elev','elev2','year','elev:year','elev2:year')
# })
params <- c('int','elev','elev2','year','elev:year','elev2:year')
betaplot <- as.mcmc.list(betaplot)
betaplot2 <- as.mcmc.list(betaplot2)
par( mfrow=c(7,2),mar=c(2,3,0,0)+0.1)
plot(betaplot, auto.layout = FALSE, cex.main=0.8, cex.sub=0.8)
par( mfrow=c(7,2),mar=c(2,3,0,0)+0.1)
plot(betaplot2, auto.layout = FALSE, cex.main=0.8, cex.sub=0.8)
dev.off()

# Gelman statistic for Beta values (slopes)
psrf.beta.pa = gelman.diag(mpost_pa$Beta, multivariate=FALSE)$psrf
hist(psrf.beta.pa)
psrf.beta.abun = gelman.diag(mpost_abun$Beta, multivariate=FALSE)$psrf
hist(psrf.beta.abun)

nm = length(models)
ma = NULL
for(j in 1:nm){
  mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
  beta <- mcmc.list2array(mpost$Beta)
  ge.beta <- sapply(1:dim(beta)[3], FUN = function(x) Rhat(sims = beta[,,x]))
  # psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  tmp = summary(ge.beta)
  if(is.null(ma)){
    ma = ge.beta#[,1]
  } else {
    ma = cbind(ma,ge.beta)#[,1])
  }
}
png(file="R/Figs/MCMC_convergence_hurdle.png", res = 600, width = 4.5, height = 6, units = "in")
# windows(4.5,6)
par(mfrow=c(2,1),las=1, mar=c(3,4,2,1)+0.1)
vioplot::vioplot(ma,col=rainbow_hcl(nm),names=c("presence-absence","abundance-cop"),
        ylim=c(0.99,max(ma)),main="psrf(beta)", ylab="Gelman statistic" )
abline(a=1.05,b=0)
vioplot::vioplot(ma,col=rainbow_hcl(nm),names=c("presence-absence","abundance-cop"),
        ylim=c(0.99,1.1),main="psrf(beta)", ylab="Gelman statistic" )
abline(a=1.05,b=0)
dev.off()
# what proportion is above 1.05?
1-apply(ma, 2, function(z) sum((z > 1.05)) ) / nrow(ma)



# ## Assess model fit #####
# preds = lapply( models, computePredictedValues )
# MF <- mapply( evaluateModelFit, models, preds )
# # windows(5,5)
# for(i in 1:length(models)){
#   print(lapply( MF[[i]], summary))
# }




# ## parameter estimates ####
# postBeta = lapply( models, getPostEstimate, parName = "Beta")
# # windows(5,8)
# # plotBeta(m, post = postBeta, param = "Sign", supportLevel = 0.95, mar=c(7,11,0,6))
# postBeta[[1]]$mean[, c("Alaria.marginata","Hedophyllum.sessile","Polysiphonia")]
# postBeta[[2]]$mean[, c("Alaria.marginata","Hedophyllum.sessile","Polysiphonia")]
# 
# cor( as.vector(postBeta[[1]]$mean), as.vector(postBeta[[2]]$mean) )
# plot( as.vector(postBeta[[1]]$mean), as.vector(postBeta[[2]]$mean) )
# 
# pos.negs <- NULL
# for(i in 1:length(models)){
#   pos.neg <- data.frame(pos = c(postBeta[[i]]$support), neg = c(postBeta[[i]]$supportNeg))
#   pos.neg[pos.neg< 0.95] <- 0
#   pos.neg$neg <- -pos.neg$neg
#   pos.neg$value <- pos.neg$pos + pos.neg$neg
#   names(pos.neg) <- paste0( c('pos','neg','value'),i)
#   if(is.null(pos.negs)){
#     pos.negs=pos.neg#[,1]
#   } else {
#     pos.negs = cbind(pos.negs,pos.neg)#[,1])
#   }
#   pos.negs$parameter <- factor(c("intercept", "year1", "year2", "elev1",
#                                 "elev2", "elev1:year1", "elev1:year2" ),
#                               levels = c("intercept", "year1", "year2", "elev1",
#                                          "elev2", "elev1:year1", "elev1:year2" ),
#                               ordered = TRUE)
#   pos.negs$species <- factor(rep(colnames(postBeta[[i]]$mean), each = 7),
#                             levels = colnames(models[[i]]$Y)[order(colSums(models[[i]]$Y),decreasing = TRUE)],
#                             ordered = TRUE)
# }
# 
# 
# windows(12,8)
# support1 <- ggplot(pos.negs, aes(y = parameter, x = species, fill = value1))+
#   geom_tile()+
#   scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
#   theme(axis.text.x = element_text(angle = 90,size=7.5))+
#   xlab(label = "")+
#   ylab(label = "")
# support2 <- ggplot(pos.negs, aes(y = parameter, x = species, fill = value2))+
#   geom_tile()+
#   scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
#   theme(axis.text.x = element_text(angle = 90,size=7.5))+
#   xlab(label = "")+
#   ylab(label = "")
# plot_grid( support1, support2, ncol=1 )
# 
# 
# # # pull taxa with positive estimated response to temperature anomaly
# # taxa.pos.anom <- as.character(pos.neg[ pos.neg$parameter=='temp.anom' & pos.neg$value>0, 'species'])
# 
# # calculate mean and variance of parameter esimates
# pbdf <- data.frame( t(postBeta$mean), taxon=colnames(postBeta$mean) )
# names(pbdf) <- c("intercept","elev","elev2","year",
#                  "elev:year","elev2:year","taxon")
# ## Add some basic trait information
# pbdf$alga <- "alga"
# pbdf$alga[c(3,15,20,56,57,62,68,82)] <- "invert"
# # More specific groups
# 
# 
# coef_plot <- pbdf %>%
#   select( -taxon ) %>%
#   group_by( alga ) %>%
#   gather( coefficient, value, -alga )
# 
# windows(6,4)
# ggplot( coef_plot, aes(y=value, x=factor(1), col=alga)) +
#   facet_wrap(~coefficient, scales="free_y") +
#   geom_hline(yintercept=0) +
#   stat_summary( fun.data=mean_cl_boot, geom='errorbar',
#                 size=1, width=0.1, position=position_dodge(width=0.9) )  +
#   stat_summary( fun.y=mean, geom='point',
#                 size=3, position=position_dodge(width=0.9) )  +
#   ylab('coefficient') + xlab('') +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())





## variance partitioning ####
VP = lapply( models, computeVariancePartitioning ) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP)
VP.dfs <- NULL
for( i in 1:length(models) ){
  VP.df <- as.data.frame(VP[[i]]$vals) %>%
    mutate(effect = factor(c("year1","year2","elev1","elev2",
                             "elev1:year1","elev2:year1",
                             "site","transect","quadrat"),
                           levels = rev(c("year1","year2","elev1","elev2",
                                          "elev1:year1","elev2:year1",
                                          "site","transect","quadrat")),
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
    mutate(tempR2 = variance[effect == "year1"])
  if(is.null(VP.dfs)){
    VP.dfs=VP.df#[,1]
  } else {
    VP.dfs = list(VP.dfs,VP.df)#[,1])
  }
  # mutate(tempR2 = variance[effect == "temp.anom.sum"])
}
VP.df <- do.call(rbind, VP.dfs)
VP.df$model <- gl( n = 2,k = nrow(VP.df)/2,labels = c("presence-absence","abunance_cop") )

hold <- VP.df %>% filter(effect == "year1") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species,
                        levels = colnames(models[[i]]$Y)[order(colSums(models[[i]]$Y),decreasing = TRUE)],
                        ordered = TRUE)


R2.df2 <- data.frame(R2 = round(MF[[2]]$R2,1), species = colnames(models[[1]]$Y))
R2.df1 <- data.frame(R2 = round(MF[[1]]$TjurR2,1), species = colnames(models[[1]]$Y))

# windows(8,5)
ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  facet_wrap(~model,ncol=1) +
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", viridis(6)), name = "")+
  # geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  # geom_point(data = R2.df, aes(y = -0.07, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")
ggsave("R/Figs/hmsc_varpart.svg", width = 8, height = 6)

# how much variation explained by fixed and random effects on average?
VP.df$ranfix <- ifelse( VP.df$effect %in% c("quadrat","transect","site"),"random","fixed")
VP.df %>% 
  group_by( model, ranfix, species ) %>% 
  summarize( variance = sum(variance) ) %>% 
  group_by( model, ranfix ) %>% 
  summarize( variance = mean(variance))

# how much variation explained by site
VP.df %>% 
  filter( effect == "site" ) %>% 
  group_by( model, ranfix ) %>% 
  summarize( variance = mean(variance))




## associations ####
OmegaCor = lapply( models, computeAssociations )
supportLevel = 0.95
# choose the random variable to plot
rlevel = 3
pick <- 1
toPlot = ((OmegaCor[[pick]][[rlevel]]$support>supportLevel)
          + (OmegaCor[[pick]][[rlevel]]$support<(1-supportLevel))>0)*OmegaCor[[pick]][[rlevel]]$mean
# reorder species matrix
plotorder <- order( postBeta[[pick]]$mean[2,], decreasing = TRUE )
toPlot <- toPlot[ plotorder, plotorder]
# rename rows for easier plotting
# newnames <- vegan::make.cepnames( rownames( toPlot) )
# newnames[[11]] <- "Boss_art"
# rownames( toPlot ) <- newnames
# colnames( toPlot ) <- newnames
# reorder automatically
mynewcor <- corReorder( toPlot, order="hclust", nclusters = 3, plot = F )
# windows(5.75,5.75)
corrplot( mynewcor, method = "color", type = "upper", tl.col="black",  
         col = colorRampPalette(c("blue","white","red"))(200),
           title = paste("random effect level:", models[[1]]$rLNames[rlevel]), mar=c(0,0,0.5,0), tl.cex=0.6 )







#







#








##### merge community data and traits #####
comm_select <- metacomm %>% 
  gather(key = taxon, value = N, Acrosiphonia:Unknown.red.blade) %>%
  select( UID, shore.height=Shore_height_cm, taxon, N ) %>% 
  separate( UID,c("beach","b","zone","year","quad") ) %>% 
  unite( site, beach, b )
comm_select_cond <- metacomm %>% 
  gather(key = taxon, value = N, Acrosiphonia:Unknown.red.blade) %>%
  select( UID, shore.height=Shore_height_cm, taxon, N ) %>%
  filter( N>0 ) %>% 
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
comm_final$taxon <- factor(comm_final$taxon, 
                           levels = colnames(models[[1]]$Y)[order(colSums(models[[2]]$Y),decreasing = T)], ordered = TRUE)

inverts <- c("Barnacles","Mytilus.sp.","Anemone","Bryozoan","Tunicata/Porifera","Pollicipes.polymerus","Tube.worms","Hydroid" )
comm_final$alga <- "alga"
comm_final$alga[ comm_final$taxon %in% inverts ] <- "invert"
comm_final$elev <- comm_final$shore.height

# filter out zeros to get conditional abundances
comm_final_cond <- comm_final %>%
  filter( N>0 )
# convert to presence absence
comm_final_pa <- comm_final %>%
  mutate( N = ifelse(N>0,1,0) )
# predictions_abund$alga <- "alga"
# predictions_abund$alga[ predictions_abund$taxon %in% inverts ] <- "invert"
# merge with other traits

# occurrences across years
commtab <- with(comm_select, table(taxon, year))
commtabdf <- as.data.frame( commtab )
commtab_wide <-commtabdf %>% 
  spread( taxon, Freq )
# # commtab <- commtab[ order(rowSums(commtab),decreasing = T), ]
# write.table( commtab, "R/output/occurrence_table.txt")

## add trait information
taxa  <- read_csv( "Data/taxa/TaxonList_corrected_lumped_unique.csv" )
FunGroups <- read_csv( "Data/taxa/Algae_functional_groups.csv" )
FunGroups$taxon<-gsub(" ",".",FunGroups$taxon)
# group seagrass with large browns
# FunGroups$funct_Sep2020 <- gsub("large_brown","canopy",FunGroups$funct_Sep2020)
# FunGroups$funct_Sep2020 <- gsub("seagrass","canopy",FunGroups$funct_Sep2020)
# FunGroups$taxon[ FunGroups$taxon == "Bossiella.articulate" ] <- "Bossiella_articulate"
# fix naming discrepancies
# names(comm)[colnames(comm)=="Bossiella_articulate"] <- "Bossiella.articulate"

# Creat2 matrix of functional groups summed
taxon<-data.frame(taxon = colnames(comm))
taxon.key<-left_join(taxon, FunGroups, by = "taxon") # animals should be missing
taxon.key$funct <- taxon.key$funct_2021
taxon.key$funct [ is.na(taxon.key$funct )] <- "animal"
taxon.key <- left_join( data.frame(taxon = colnames(models[[1]]$Y)), taxon.key )



#####
# Model  predictions across elevations
# predict any elevation
XData_choose <- models[[1]]$XData %>%
  select(year,year1,year2) %>% mutate_all(round,7) %>%  distinct()
XData_choose_elev <- models[[1]]$XData %>%
  select(elev,elev1,elev2) %>% mutate_all(round,7) %>%  distinct() %>% arrange(elev)
# get a model for elev2 as functions of elev and elev1 so we can predict any elevation
plot(XData_choose_elev)
lm_raw <- lm( elev2 ~ poly(elev, 2, raw = T), data = XData_choose_elev )
lm_scale <- lm( elev2 ~ poly(elev1, 2, raw = T), data = XData_choose_elev )
lm_linear <- lm( elev1 ~ elev, data = XData_choose_elev )
predict( lm_raw, newdata = data.frame(elev = 61))
predict( lm_scale, newdata = data.frame(elev1 = predict(lm_linear, newdata = data.frame(elev = 61))))
# response curves across all combinations
years <- as.data.frame( models[[1]]$XData %>% select(year,year1,year2) %>% mutate(year1=round(year1,7),year2=round(year2,7)) %>% distinct() )
elevs <- as.data.frame( models[[1]]$XData %>% select(elev,elev1,elev2) %>% mutate(elev1=round(elev1,5),elev2=round(elev2,5)) %>% distinct() %>% arrange(elev1) )
# use lm_linear, etc. to predict any elevation
quantile(models[[1]]$XData$elev, probs = c(0.01,0.25, 0.5, 0.75, 0.99))
elev_length = 100
# elev_grad <- round((seq((69), (379), length = elev_length)))
elev_grad <- round((seq(min(models[[1]]$XData$elev), max(models[[1]]$XData$elev), length = elev_length)))
hist(models[[1]]$XData$elev)
abline(v = elev_grad, col='blue', lty = 4)
# elev_grad_bins <- elev_grad[seq(1,elev_length, by = 2)]
# elev_bin_prop <- table(cut(models[[1]]$XData$elev, elev_grad_bins, right = FALSE))
# elev_bin_prop <- elev_bin_prop/sum(elev_bin_prop)
# elev_grad <-  elev_grad[seq(2,elev_length, by = 2)]
# abline(v = elev_grad, col='slateblue', lty=2)
# elev_grad <- -50:500
# elev_grad <- elev_grad[seq(1,length(elev_grad), by = 6)]
elev_grad2 <- data.frame( elev = elev_grad,
                          elev1 = predict( lm_linear, newdata = data.frame(elev = elev_grad)),
                          elev2 = predict( lm_scale, newdata = data.frame(elev1 = predict(lm_linear, newdata = data.frame(elev = elev_grad)))) )

# Gradient <- constructGradient(models[[1]], focalVariable = 'year1', non.focalVariables = list('elev1' = list(1)),
#                               ngrid = 8 )
# Gradient$XDataNew$year2 <- XData_choose$year2
# Gradient$XDataNew$elev1 <- elev_grad2$elev1[i]
# Gradient$XDataNew$elev2 <- elev_grad2$elev2[i]

newDF <- data.frame( merge(XData_choose, elev_grad2),
                     site = "new unit",
                     transect = "new unit",
                     quadrat = "new unit")

# newDF <- expand.grid( shore.height = seq(60,380,by = 1),
#              year=c(2012:2019),
#              site="new unit",
#              transect="new unit" )
# newDF <- data.frame( merge(years, elevs),
#              site="new unit",
#              transect="new unit" )
# trim dataset to every nth observation, only particular years
# newDF <- data.frame( merge(years[years$year %in% c(2012:2019),], elevs[seq(1,170, by = 3),]),
#                      site="new unit",
#                      transect="new unit" )

# newDFsel <- newDF %>% select(year1,year2,elev1,elev2,site,transect,quadrat)
newXData   <- newDF[,1:(ncol(newDF)-3)]
newDesign  <- data.frame( site = as.factor(newDF$site), transect = as.factor(newDF$transect),
                          quadrat = as.factor(newDF$quadrat))
names(models[[1]]$XData)
names(newXData)
newXData <- newXData %>% 
  select(year,elev,year1,year2,elev1,elev2)


## predictions of individual models
predY_pa <- predict(models[[1]], XData = newXData,
                 studyDesign = newDesign,
                 ranLevels = list(site = models[[1]]$rL$site, transect = models[[1]]$rL$transect, quadrat = models[[1]]$rL$quadrat), expected = TRUE)
predY_cop <- predict(models[[2]], XData = newXData,
                     studyDesign= newDesign,
                     ranLevels = list(site = models[[1]]$rL$site, transect = models[[1]]$rL$transect, quadrat = models[[1]]$rL$quadrat), expected = TRUE) #, transect=rL, year=rL_year
## predictions of both models multiplied together
predY_abun <- Map('*', predY_pa, lapply(predY_cop,exp) )
# predY_pa[[1]][1,1] * predY_cop[[1]][1,1]
# predY_abun[[1]][1,1]

# ###
### other ways to summarize the model results --------------------------------------------------------------------
tmp = abind::abind(predY_pa, along = 3)
qpred = apply(tmp, c(1, 2), quantile, prob = 0.5, na.rm = TRUE)
predictions_pa <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.025, na.rm = TRUE)
predictions_pa_low <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.975, na.rm = TRUE)
predictions_pa_high <- bind_cols(data.frame(qpred), newXData)

tmp = abind::abind(lapply(predY_cop,exp), along = 3)
qpred = apply(tmp, c(1, 2), quantile, prob = 0.5, na.rm = TRUE)
predictions_cop <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.025, na.rm = TRUE)
predictions_cop_low <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.975, na.rm = TRUE)
predictions_cop_high <- bind_cols(data.frame(qpred), newXData)

tmp = abind::abind(predY_abun, along = 3)
qpred = apply(tmp, c(1, 2), quantile, prob = 0.5, na.rm = TRUE)
predictions_abun <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.025, na.rm = TRUE)
predictions_abun_low <- bind_cols(data.frame(qpred), newXData)

qpred = apply(tmp, c(1, 2), quantile, prob = 0.975, na.rm = TRUE)
predictions_abun_high <- bind_cols(data.frame(qpred), newXData)




# color scheme
as.survey <- read_csv(  "R/output/sst_anoms_survey.csv" )
as.survey$year <-  as.survey$survey.year
library( RColorBrewer )
anom.range <- c(-2,2)
n=9
cols <- brewer.pal(n,"RdBu")
pal <- colorRampPalette(rev(cols))

predictions_pa <- left_join(predictions_pa, as.survey)
predictions_cop <- left_join(predictions_cop, as.survey)
predictions_abun <- left_join(predictions_abun, as.survey)

### colors
cols.two <- c( rgb( 211,230,240, maxColorValue=255), rgb( 232,139,110, maxColorValue=255))
cols.two <- c( "#BDCFD8",  "#D07D63" )
# use info from the pca1 in script "Community_rda.R"
# 2012 pca1 is -1.21
# 2019 pca1 is 1.52
data.frame( pca1 = seq(-2.84,2.84, length=100), hex = pal2(100) )
cols.two <- c( "#A2CDE2",  "#EF9B7A" )
cols.two <- c( "#A2CDE2",  "#BF7C61" ) # make the redder one darker

# add empirical means
models[[1]]$studyDesign %>% 
  group_by(year) %>% 
  summarize(ntrans=length(transect))
comm <- comm[ !is.na(comm$Shore_height_cm), ]
ogd <- bind_cols(  models[[1]]$studyDesign, comm )  #data.frame(m$Y)
ogd$year <- ogd$year

ogd.long <- ogd %>% 
  pivot_longer( cols="Acrosiphonia":"Unknown.crust", 
                names_to = "taxon", values_to = "N" ) %>% 
  select( UID, year, elev = Shore_height_cm, taxon, N )
ogd.long$N[ ogd.long$N %in% c(0.1,0.25) ] <- 0.5
ogd.cop <- ogd.long
ogd.cop$N[ ogd.cop$N == 0 ] <-  NA
ogd.pa <- ogd.long
ogd.pa$N <- as.numeric( ogd.pa$N > 0 )

# make wide again t4o matching plotting style below
ogd.wide <- ogd.long %>% 
  pivot_wider( names_from = taxon, values_from = N )
ogd.cop.wide <- ogd.cop %>% 
  pivot_wider( names_from = taxon, values_from = N )
ogd.pa.wide <- ogd.pa %>% 
  pivot_wider( names_from = taxon, values_from = N )


# plot responses
predictions_pa$yearplot <- as.character(factor(predictions_pa$year1, levels = unique(predictions_pa$year1), labels = unique(newDF$year)))
predictions_cop$yearplot <- as.character(factor(predictions_pa$year1, levels = unique(predictions_pa$year1), labels = unique(newDF$year)))
predictions_abun$yearplot <- as.character(factor(predictions_pa$year1, levels = unique(predictions_pa$year1), labels = unique(newDF$year)))
# predictions_pa$elev <- as.numeric(as.character(factor(predictions_pa$elev1, levels = unique(predictions_pa$elev1), labels = unique(newDF$elev))))
# predictions_cop$elev <- as.numeric(as.character(factor(predictions_pa$elev1, levels = unique(predictions_pa$elev1), labels = unique(newDF$elev))))
# predictions_abun$elev <- as.numeric(as.character(factor(predictions_pa$elev1, levels = unique(predictions_pa$elev1), labels = unique(newDF$elev))))


hist( comm$Fucus.distichus )
# taxon <- "Hedophyllum.sessile"
taxon <- "Alaria.marginata"
# taxon <- "Fucus.distichus"
# taxon <- "Mytilus.sp."

lwds = 0.75
sizes = 1

a <- ggplot(filter(predictions_pa, year %in% c(2012,2019), elev >= range(models[[1]]$XData$elev)[1], elev <= range(models[[1]]$XData$elev)[2] ), 
            aes_string(x = 'elev', y = taxon, col='yearplot'))+
  geom_point( data = filter(ogd.pa.wide, year %in% c(2012,2019), elev >= range(models[[1]]$XData$elev)[1], elev <= range(models[[1]]$XData$elev)[2] ), 
              aes_string(x = 'elev', y = taxon, col='year'), alpha = 0.5, size = sizes ) +
  geom_line(lwd = lwds) +
  scale_color_manual(values=cols.two) +
  coord_cartesian( ylim = c(0,1) ) +
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "none") +
  # theme(legend.position = c(1,0.25), legend.justification = c(1,0), legend.title=element_blank() ) +
  # guides(col=guide_legend("year"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
b <- ggplot(filter(predictions_cop, year %in% c(2012,2019), elev >= range(models[[1]]$XData$elev)[1], elev <= range(models[[1]]$XData$elev)[2] ), 
            aes_string(x = 'elev', y = taxon, col='yearplot'))+
  geom_point( data = filter(ogd.cop.wide, year %in% c(2012,2019), elev >= range(models[[1]]$XData$elev)[1], elev <= range(models[[1]]$XData$elev)[2] ),
              aes_string(x = 'elev', y = taxon, col='year'), alpha = 0.5, size = sizes ) +
  geom_line(lwd = lwds) +
  # scale_y_log10() +
  scale_color_manual(values=cols.two) +
  coord_cartesian( ylim = c(0,100) ) +
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
c <- ggplot(filter(predictions_abun, year %in% c(2012,2019), elev >= range(models[[1]]$XData$elev)[1], elev <= range(models[[1]]$XData$elev)[2] ), 
            aes_string(x = 'elev', y = taxon, col='yearplot'))+
  geom_point( data = filter(ogd.wide, year %in% c(2012,2019), elev >= range(models[[1]]$XData$elev)[1], elev <= range(models[[1]]$XData$elev)[2] ),
              aes_string(x = 'elev', y = taxon, col='year'), alpha = 0.5, size = sizes ) +
  geom_line(lwd = lwds) +
  # scale_y_log10() +
  scale_color_manual(values=cols.two) +
  coord_cartesian( ylim = c(0,100) ) +
  ylab("") +
  xlab("") +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
plot_grid(c,a,b,ncol=1, align='hv')
ggsave( paste0("R/Figs/hmsc_elev_metrics_",taxon,".svg"), height = 4.5, width = 1.5 ) # width 1.6 for Hedophyllum, 1.5 for others

# a <- ggplot(predictions_pa, aes_string(x = 'factor(year1,labels=2012:2019)',
#                                 y = taxon))+
#   geom_boxplot() +  ylab("probability of occurrence") + xlab('')
# b <- ggplot(predictions_cop, aes_string(x = 'factor(year1,labels=2012:2019)',
#                                         y = taxon))+
#   geom_boxplot() +  ylab("percent cover | occurrence") + xlab('')
# c <- ggplot(predictions_abun, aes_string(x = 'factor(year1,labels=2012:2019)',
#                                          y = taxon))+
#   geom_boxplot() +  ylab("percent cover") + xlab('year')
# plot_grid(a,b,c,ncol=1)
#

## tidy the data and combine
predictions_abund <- predictions_abun %>%
  gather(key = taxon, value = N, Fucus.distichus:Mazzaella.parvula)
predictions_abund_low <- predictions_abun_low %>%
  gather(key = taxon, value = N_low, Fucus.distichus:Mazzaella.parvula)
predictions_abund_high <- predictions_abun_high %>%
  gather(key = taxon, value = N_high, Fucus.distichus:Mazzaella.parvula)
predictions_abund <- left_join(left_join(predictions_abund, predictions_abund_low), predictions_abund_high)
predictions_abund$taxon <- factor(predictions_abund$taxon, levels = colnames(Y)[order(colSums(Y),decreasing = T)], ordered = FALSE)
# abundance conditional on presence
predictions_copp <- predictions_cop %>%
  gather(key = taxon, value = N, Fucus.distichus:Mazzaella.parvula)
predictions_copp_low <- predictions_cop_low %>%
  gather(key = taxon, value = N_low, Fucus.distichus:Mazzaella.parvula)
predictions_copp_high <- predictions_cop_high %>%
  gather(key = taxon, value = N_high, Fucus.distichus:Mazzaella.parvula)
predictions_copp <- left_join(left_join(predictions_copp, predictions_copp_low), predictions_copp_high)
predictions_copp$taxon <- factor(predictions_copp$taxon, levels = colnames(Y)[order(colSums(Y),decreasing = T)], ordered = FALSE)
# presence-absence
predictions_pab <- predictions_pa %>%
  gather(key = taxon, value = N, Fucus.distichus:Mazzaella.parvula)
predictions_pab_low <- predictions_pa_low %>%
  gather(key = taxon, value = N_low, Fucus.distichus:Mazzaella.parvula)
predictions_pab_high <- predictions_pa_high %>%
  gather(key = taxon, value = N_high, Fucus.distichus:Mazzaella.parvula)
predictions_pab <- left_join(left_join(predictions_pab, predictions_pab_low), predictions_pab_high)
predictions_pab$taxon <- factor(predictions_pab$taxon, levels = colnames(models[[1]]$Y)[order(colSums(models[[1]]$Y),decreasing = T)], ordered = FALSE)

###


# load( file = paste("R/output/hmsc_pred",model, sep="_") )






##### Counter-factual predictions ####
## predictions over on factor only
# get predictions over time
Gradient <- constructGradient(models[[1]], focalVariable = "year1", non.focalVariables = list("elev1" = list(1)),
                              ngrid = 8) #length(unique(models[[1]]$XData$year1)))

Gradient$XDataNew$year1

XData_choose_elev %>% filter( elev2==min(elev2))
mean(XData_choose_elev$elev)
Gradient$XDataNew$year2 <- XData_choose$year2#[ XData_choose$year1 == round(Gradient$XDataNew$year1,7)]

# middle of low zone
mean(quantile( elevs$elev, c(0,0.33) )) # middle of lower zone (although for Fifth Beach this is very close to the lower limit)
Gradient$XDataNew$elev1 <- -0.04676
Gradient$XDataNew$elev2 <- 0.02619

# middle of shore
Gradient$XDataNew$elev1 <- 0.00011
Gradient$XDataNew$elev2 <- -0.03570

years <- data.frame( year1 = round(Gradient$XDataNew$year1,7), year = 2012:2019 )

# Gradient$studyDesignNew$Time <- unique(mod_list[[1]]$studyDesign$Time)
# Gradient$rLNew$Time <- mod_list[[1]]$rL$Time

Time_predY1 <- predict(models[[1]], XData = Gradient$XDataNew, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = TRUE)
Time_predY2 <- predict(models[[2]], XData = Gradient$XDataNew, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = TRUE)

Time_predY1 <- abind::abind(Time_predY1, along = 3)
Time_predY2 <- exp(abind::abind(Time_predY2, along = 3))
Time_predY3 <- Time_predY1*Time_predY2

scaled_occur <- Time_predY1
log_con_bmass <- log(Time_predY2)
log_biomass <- log(Time_predY3)



# gammas
# Gammas
postGamma = lapply( models, getPostEstimate, parName = "Gamma")
plotGamma(models[[1]], post=postGamma[[1]], param="Support", supportLevel = 0.95, covNamesNumbers = c(T,F), trNamesNumbers = c(T,F), colorLevels = 3 )
plotGamma(models[[2]], post=postGamma[[2]], param="Support", supportLevel = 0.95, covNamesNumbers = c(T,F), trNamesNumbers = c(T,F), colorLevels = 3 )
plotGamma(models[[2]], post=postGamma[[2]], param="Mean", supportLevel = 0.5, covNamesNumbers = c(T,F), trNamesNumbers = c(T,F), colorLevels = 3)
#
pick = 1
hM <- models[[pick]]
# Gradient = constructGradient(models[[pick]],focalVariable = "year1")
predY = predict(hM,Gradient = Gradient, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = TRUE)
prob = c(0.25,0.5,0.75)
# plotGradient
png(file="R/Figs/HMSC_species_richness.png", res = 600, width = 4, height = 4, units = "in")
plotGradient(hM, Gradient, pred=predY, measure="S", showData = TRUE, q = prob, las = 1, axes = F, xlab="Year")
axis(2,las = 1)
axis(1, at = Gradient$XDataNew$year1, labels = 2012:2019 )
dev.off()
plotGradient(hM, Gradient, pred=predY, measure="S", index=1, showData = TRUE, q = prob, axes=F) # prob should be q
axis(1, at = years$year1, labels = years$year )
axis(2, las = 1)
plotGradient(hM, Gradient, pred=predY, measure="Y", index=1, showData = TRUE, q = prob) # prob should be q
plotGradient(hM, Gradient, pred=predY, measure="T", index=1, showData = TRUE,  q = prob) # prob should be q
plotGradient(hM, Gradient, pred=predY, measure="T", index=2, showData = TRUE,  q = prob) # prob should be q
plotGradient(hM, Gradient, pred=predY, measure="T", index=3, showData = TRUE,  q = prob) # prob should be q
plotGradient(hM, Gradient, pred=predY, measure="T", index=4, showData = TRUE,  q = prob) # prob should be q
plotGradient(hM, Gradient, pred=predY, measure="T", index=5, showData = TRUE,  q = prob) # prob should be q
plotGradient(hM, Gradient, pred=predY, measure="T", index=6, showData = TRUE,  q = prob) # prob should be q








#####
# code from Patrick Thompson
# 
quantile(models[[1]]$XData$elev, probs = c(0.01,0.25, 0.5, 0.75, 0.99))
hist(models[[1]]$XData$elev)
elev_length = 8
elev_grad <- round((seq((69), (379), length = elev_length)))
# elev_grad_bins <- elev_grad[seq(1,elev_length, by = 2)]
elev_grad_bins <- elev_grad
elev_bin_prop <- table(cut(models[[1]]$XData$elev, elev_grad_bins, right = FALSE))
elev_bin_prop <- elev_bin_prop/sum(elev_bin_prop)
round(elev_bin_prop,2)
# elev_grad <-  elev_grad[seq(2,elev_length, by = 2)]
abline(v = elev_grad, col='slateblue', lty=2)
elev_grad2 <- data.frame( elev = elev_grad,
                          elev1 = predict( lm_linear, newdata = data.frame(elev = elev_grad)),
                          elev2 = predict( lm_scale, newdata = data.frame(elev1 = predict(lm_linear, newdata = data.frame(elev = elev_grad)))) )
occur_all_list <- list()
cond_bmass_list <- list()
bmass_list <- list()
for(i in 1:length(elev_grad)){
  Gradient <- constructGradient(models[[1]], focalVariable = 'year1', non.focalVariables = list('elev1' = list(1)),
                                ngrid = 8 )
  Gradient$XDataNew$year2 <- XData_choose$year2
  Gradient$XDataNew$elev1 <- elev_grad2$elev1[i]
  Gradient$XDataNew$elev2 <- elev_grad2$elev2[i]
  Time_predY1 <- predict(models[[1]], XData = Gradient$XDataNew, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = TRUE)
  Time_predY2 <- predict(models[[2]], XData = Gradient$XDataNew, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = TRUE)
  Time_predY1 <- abind::abind(Time_predY1, along = 3)
  Time_predY2 <- exp(abind::abind(Time_predY2, along = 3))
  Time_predY3 <- Time_predY1*Time_predY2
  occur_all_list[[i]] <- Time_predY1
  cond_bmass_list[[i]] <- Time_predY2
  bmass_list[[i]] <- Time_predY3
}
occur_mean <- occur_all_list[[1]]
cond_bmass_mean <- cond_bmass_list[[1]]
bmass_mean <- bmass_list[[1]]
for(i in 2:length(elev_grad)){
  occur_mean <- occur_mean + occur_all_list[[i]] * elev_bin_prop[i-1]
  cond_bmass_mean <- cond_bmass_mean + cond_bmass_list[[i]] * elev_bin_prop[i-1]
  bmass_mean <- bmass_mean + bmass_list[[i]] * elev_bin_prop[i-1]
}
scaled_occur <- occur_mean
log_con_bmass <- log(cond_bmass_mean)
log_biomass <- log(bmass_mean)

species_occur <- left_join(apply(scaled_occur, c(1,2), median) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "occurrence") %>% gather(key = Species, value = median, -year, -metric),
                           apply(scaled_occur, c(1,2), quantile, prob = 0.25) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "occurrence") %>% gather(key = Species, value = quant_0.25, -year, -metric)) %>%
  left_join(apply(scaled_occur, c(1,2), quantile, prob = 0.75) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "occurrence") %>% gather(key = Species, value = quant_0.75, -year, -metric))

species_occur_scale <- species_occur %>%
  gather(key = quant, value = value, median:quant_0.75) %>%
  mutate(value = scale(value)) %>%
  spread(key = quant, value = value)

species_con_bmass <- left_join(apply(log_con_bmass, c(1,2), median) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "conditional cover") %>% gather(key = Species, value = median, -year, -metric),
                               apply(log_con_bmass, c(1,2), quantile, prob = 0.25) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "conditional cover") %>% gather(key = Species, value = quant_0.25, -year, -metric)) %>%
  left_join(apply(log_con_bmass, c(1,2), quantile, prob = 0.75) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "conditional cover") %>% gather(key = Species, value = quant_0.75, -year, -metric))

species_con_bmass_scale <- species_con_bmass %>%
  gather(key = quant, value = value, median:quant_0.75) %>%
  mutate(value = scale(value)) %>%
  spread(key = quant, value = value)

species_bmass <- left_join(apply(log_biomass, c(1,2), median) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "cover") %>% gather(key = Species, value = median, -year, -metric),
                           apply(log_biomass, c(1,2), quantile, prob = 0.25) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "cover") %>% gather(key = Species, value = quant_0.25, -year, -metric)) %>%
  left_join(apply(log_biomass, c(1,2), quantile, prob = 0.75) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "cover") %>% gather(key = Species, value = quant_0.75, -year, -metric))

species_bmass_scale <- species_bmass %>%
  gather(key = quant, value = value, median:quant_0.75) %>%
  mutate(value = scale(value)) %>%
  spread(key = quant, value = value)

species_temporal.df <- rbind(species_occur_scale, species_con_bmass_scale, species_bmass)
species_temporal.df <- left_join( mutate(species_temporal.df,year=round(year,7)), years )
species_temporal.df$Species <- gsub(pattern = "_", replacement = ".", x = species_temporal.df$Species )


Years <- 2012:2019
slopes <- sapply(X = 1:dim(scaled_occur)[2], FUN = function(i){
  sapply(X = 1:dim(scaled_occur)[3], FUN = function(x){
    c(coef(lm(scaled_occur[,i,x]~Years))[2])
  })
})
slopes_con_biomass_ln <- sapply(X = 1:dim(log_con_bmass)[2], FUN = function(i){
  sapply(X = 1:dim(log_con_bmass)[3], FUN = function(x){
    c(coef(lm(log_con_bmass[,i,x]~Years))[2])
  })
})
slopes_biomass_ln <- sapply(X = 1:dim(log_biomass)[2], FUN = function(i){
  sapply(X = 1:dim(log_biomass)[3], FUN = function(x){
    c(coef(lm(log_biomass[,i,x]~Years))[2])
  })
})
occur_trends <- data.frame(species = colnames(models[[1]]$Y),
                           median = apply(slopes, 2, median),
                           lower_0.25 = apply(slopes, 2, quantile, prob = 0.25),
                           upper_0.75 = apply(slopes, 2, quantile, prob = 0.75),
                           lower = apply(slopes, 2, quantile, prob = 0.025),
                           upper = apply(slopes, 2, quantile, prob = 0.975),
                           measure = 'Occurrence prob.')
occur_trends_taxa <- occur_trends
biomass_cond_ln_trends <- data.frame(species = colnames(models[[1]]$Y),
                                     median = apply(slopes_con_biomass_ln, 2, median),
                                     lower_0.25 = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.25),
                                     upper_0.75 = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.75),
                                     lower = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.025),
                                     upper = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.975),
                                     measure = 'Conditional cover (log)')
cop_trends_taxa <- biomass_cond_ln_trends
biomass_ln_trends <- data.frame(species = colnames(models[[1]]$Y),
                                median = apply(slopes_biomass_ln, 2, median),
                                lower_0.25 = apply(slopes_biomass_ln, 2, quantile, prob = 0.25),
                                upper_0.75 = apply(slopes_biomass_ln, 2, quantile, prob = 0.75),
                                lower = apply(slopes_biomass_ln, 2, quantile, prob = 0.025),
                                upper = apply(slopes_biomass_ln, 2, quantile, prob = 0.975),
                                measure = 'Cover (log)')
cover_trends_taxa <- biomass_ln_trends
# save biomass_ln_trends for later
biomass_ln_slopes <- slopes_biomass_ln
all_trends <- bind_rows(occur_trends, biomass_cond_ln_trends, biomass_ln_trends) %>%
  mutate(sig = sign(upper) == sign(lower)) %>%
  mutate(sign = ifelse(median > 0 & sig == TRUE, 'increasing', ifelse(median < 0 & sig == TRUE, 'decreasing', 'no trend'))) %>%
  mutate(sign = factor(sign, levels = c('decreasing', 'no trend', 'increasing'), ordered = TRUE)) %>%
  mutate(species = gsub("_",".",species)) %>%
  mutate(species = factor(species, levels = gsub("_",".",occur_trends$species[order(occur_trends$median)]), ordered = TRUE)) %>%
  mutate(measure = factor(measure, levels = c('Occurrence prob.', 'Conditional cover (log)', 'Cover (log)')))
all_trends %>%
  ggplot(aes(x = species, y = median, color = sign))+
  geom_vline(xintercept = c(levels(all_trends$species)[seq(3, models[[1]]$ns, by = 3)]), col = 'grey', size = 0.1, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin = lower_0.25, ymax = upper_0.75), width = 0, size = 1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
  geom_point()+
  coord_flip()+
  facet_wrap(~measure, scales = 'free_x')+
  scale_color_manual(values = c('blue', 'grey', 'red'), guide = FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  ylab('Change per year')+
  xlab('')
ggsave("R/Figs/hmsc_species_trends_summary.svg",width = 6.5, height = 7)



## plot trends for all species
# species_temporal.df$Species <- factor(species_temporal.df$Species, levels = order_occur$Species, ordered = TRUE)
# species_temporal.df$Species <- factor(species_temporal.df$Species, levels = order_occur$Species, labels = gsub("[.]","\n",order_occur$Species), ordered = TRUE)
occur_trends$species <- gsub("_",".", occur_trends$species)
species_temporal.df$metric  <- factor(species_temporal.df$metric, levels = c("occurrence","conditional cover","cover"))
# arrange by change in cover or occurrence
species_temporal.df$Species2 <- factor( species_temporal.df$Species, 
                                        # levels = gsub("_",".",occur_trends$species[order(occur_trends$median)]),
                                        levels = occur_trends$species[order(occur_trends$median)],
                                        labels = gsub( pattern = "[.]", replacement = "\n", x = occur_trends$species[order(occur_trends$median)] ),
                                        ordered = TRUE )
# adjust color scheme so that cover is darker than conditional cover
sp.trends <- ggplot(species_temporal.df, aes(x = year, y = median, color = metric,  fill = metric, group = metric))+ #fill = metric,
  geom_ribbon(aes(ymin = quant_0.25, ymax= quant_0.75), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  facet_wrap(~Species2, scales = "free", ncol = 9)+
  # scale_color_brewer(type = "qual", palette = "Dark2", name = "")+
  # scale_fill_brewer(type = "qual", palette = "Dark2", name = "", guide = FALSE)+
  scale_color_manual( values = c('gray25','seagreen','springgreen') ) +
  scale_fill_manual( values = c('gray50','seagreen','springgreen3') ) +
  xlab("Year")+
  ylab("Scaled estimate")+
  theme_classic() + 
  theme(legend.justification=c(1,0), legend.position=c(1,0),
        legend.direction="horizontal",
        legend.title=element_blank(), 
        legend.text=element_text(size=22) )
ggsave("R/Figs/hmsc_sp_trends.svg", height = 5.75*1.5, width = 9*1.5)
#


# predict total cover of producers
# remove invert columns, then total them up
colnames(models[[1]]$Y)[-c(3,15,22)]
slopes_biomass_ln_prod <- colSums(exp(log_biomass[,-c(3,15,22)]))
log_biomass_producer <- log_biomass[,-c(3,15,22),]
producer_bmass <- left_join(apply(log_biomass, c(1,2), median) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "cover") %>% gather(key = Species, value = median, -year, -metric),
                           apply(log_biomass, c(1,2), quantile, prob = 0.25) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "cover") %>% gather(key = Species, value = quant_0.25, -year, -metric)) %>%
  left_join(apply(log_biomass, c(1,2), quantile, prob = 0.75) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "cover") %>% gather(key = Species, value = quant_0.75, -year, -metric))
producer_bmass_total <- producer_bmass %>% 
  group_by( year ) %>% 
  summarize( median = sum(exp(median)),
             quant_0.25 = sum(exp(quant_0.25)),
             quant_0.75 = sum(exp(quant_0.75)) )




# # now predict what each transect did in each year
# Gradient_transect <- Gradient
# studyDesignNew_transect <- models[[1]]$studyDesign %>% 
#   select( site, transect ) %>% 
#   distinct() %>% 
#   mutate( quadrat = "new unit" )
# 
# bmass_list <- list()
# bmass_list_tmp <- list()
# for(i in 1:length(elev_grad)){
#   for(j in 1:length( unique(models[[1]]$studyDesign$transect) ) ){
#     Gradient_transect <- constructGradient(models[[1]], focalVariable = 'year1', non.focalVariables = list('elev1' = list(1)),
#                                   ngrid = 8 )
#     Gradient_transect$XDataNew$year2 <- XData_choose$year2
#     Gradient_transect$XDataNew$elev1 <- elev_grad2$elev1[i]
#     Gradient_transect$XDataNew$elev2 <- elev_grad2$elev2[i]
#     Gradient_transect$studyDesignNew$site <- studyDesignNew_transect$site[j]
#     Gradient_transect$studyDesignNew$transect <- studyDesignNew_transect$transect[j]
#     Time_predY1 <- predict(models[[1]], XData = Gradient_transect$XDataNew, studyDesign = Gradient_transect$studyDesignNew, ranLevels = Gradient_transect$rLNew, expected = TRUE)
#     Time_predY2 <- predict(models[[2]], XData = Gradient_transect$XDataNew, studyDesign = Gradient_transect$studyDesignNew, ranLevels = Gradient_transect$rLNew, expected = TRUE)
#     Time_predY1 <- abind::abind(Time_predY1, along = 3)
#     Time_predY2 <- exp(abind::abind(Time_predY2, along = 3))
#     Time_predY3 <- Time_predY1*Time_predY2
#     log_biomass_producer <- Time_predY3[,-c(3,15,22),]
#     producer_bmass <- apply(log_biomass_producer, c(1,2), median) %>%  as.data.frame() %>% mutate(year = 2012:2019, metric = "cover") %>% gather(key = Species, value = median, -year, -metric)
#     producer_bmass_total_tmp <- producer_bmass %>% 
#       group_by( year ) %>% 
#       summarize( median = sum((median)) )
#     producer_bmass_total_tmp$transect <- unique(models[[1]]$studyDesign$transect)[j]
#     bmass_list_tmp[[j]] <- producer_bmass_total_tmp
#   }
#   bmass_list[[i]] <- do.call(rbind,bmass_list_tmp)
# }
# # bmass_list <- lapply(bmass_list, function(z) do.call(rbind,z) )
# bmass_mean <- bmass_list[[1]]
# for(i in 2:length(elev_grad)){
#   bmass_mean$median <- bmass_mean$median + bmass_list[[i]]$median * elev_bin_prop[i-1]
# }
# #
# 
# ggplot( bmass_mean, aes(x = year, y = median, group = transect)) +
#   # geom_path() +
#   geom_smooth( method = 'lm', se = F) +
#   geom_ribbon( data = producer_bmass_total, aes(ymin = quant_0.25, ymax = quant_0.75,group=1), 
#                col = "slategrey", alpha = 0.5) +
#   geom_path( data = producer_bmass_total, aes(ymin = quant_0.25, ymax = quant_0.75,group=1),
#              lwd=2) +
#   # scale_y_continuous(trans = 'log') +
#   theme_classic()




## now predict traits ####
##### code modified from function plotGradient() 
hM = models[[1]]
all(hM$distr[, 1] == 1)
# make a new version of the intercept here to show canopy taxa
which( hM$TrData == "canopy")
new_traits <- hM$Tr
new_traits[,1] <- 0
new_traits[which( hM$TrData == "canopy"),1] <- 1

split.along.dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}

predYO <- split.along.dim(scaled_occur, n=3)
predYC <- split.along.dim(log_con_bmass, n=3)
predYT <- split.along.dim(log_biomass, n=3)

# presence-absence
predTO = lapply(predYO, function(a) (a %*% new_traits)/matrix(rep(rowSums(a),
                                                                  hM$nt), ncol = hM$nt))
predTO_scale <- lapply( predTO, function(z) data.frame(matrix(scale(c(z)), ncol=6) )) 
# conditional cover/biomass/abundance 
predTC = lapply(predYC, function(a) (exp(a) %*% new_traits)/matrix(rep(rowSums(exp(a)), 
                                                                     hM$nt), ncol = hM$nt))
predTC = lapply(predYC, function(a) (exp(a) %*% new_traits))
predTC_log_scale <- lapply( predTC, function(z) data.frame(matrix(scale(log(c(z))), ncol=6) )) 
# total cover
predTT = lapply(predYT, function(a) (exp(a) %*% new_traits))
predTT_log_scale <- lapply( predTT, function(z) data.frame(matrix(scale(log(c(z))), ncol=6) )) 

# index = 2
predTO_scale = abind::abind(predTO_scale, along = 3)
predTC_log_scale = abind::abind(predTC_log_scale, along = 3)
predTT_log_scale = abind::abind(predTT_log_scale, along = 3)
predTO = abind::abind(predTO, along = 3)
predTC = abind::abind(predTC, along = 3)
predTT = abind::abind(predTT, along = 3)
# xx = Gradient$XDataNew[, 1]
# ngrid = length(xx)
# Pr = mean(predT[ngrid, 1, ] > predT[1, 1, ])
q = prob = c(0.25,0.5,0.75)
# qpred = apply(predT, c(1, 2), quantile, probs = q, na.rm = TRUE)
# qpred = qpred[, , index]
# lo = qpred[1, ]
# hi = qpred[3, ]
# me = qpred[2, ]
# lo1 = min(lo, na.rm = TRUE)
# hi1 = max(hi,  na.rm = TRUE)
# ylabel = hM$trNames[[index]]
# xlabel = colnames(Gradient$XDataNew)[[1]]
# plot(2012:2019, qpred[2, ], ylim = c(lo1, hi1), type = "l", 
#      xlab = "Year", ylab = ylabel)
# polygon(c(2012:2019, rev(2012:2019)), c(qpred[1, ], rev(qpred[3, ])), 
#         col = 'slateblue1', border = FALSE)
# lines(2012:2019, qpred[2, ], lwd = 2)

# combine qpreds to make some different plots
qpred = apply(predTT, c(1, 2), quantile, probs = q, na.rm = TRUE)
qpred_list <- split.along.dim( qpred, n = 3 )
qpred_bind <- do.call(rbind, qpred_list)
qpred_long <- qpred_bind %>% c() %>% 
  data.frame() 
names(qpred_long) = "value"
qpred_long$year <- gl(8,3*6, labels = 2012:2019)
qpred_long$FG <- gl(6,3, labels = levels(models[[1]]$TrData$FG) )
qpred_long$quantile <- gl(3,1, labels = c("25%","50%","75%"))
qpred_median <- qpred_long %>% 
  filter( quantile == "50%" )
FGcolor <- c("midnightblue","darkred","red","pink","darkgrey","#996633")
ggplot( data = qpred_median, aes(x = year, y = value, group = FG, fill = FG)) +
  geom_area(col="black") +
  scale_fill_manual( values = rev(FGcolor) )
  
# make slopes for traits like with species predictions
Years <- 2012:2019

slopes <- sapply(X = 1:dim(predTO_scale)[2], FUN = function(i){
  sapply(X = 1:dim(predTO_scale)[3], FUN = function(x){
    c(coef(lm(predTO_scale[,i,x]~Years))[2])
  })
})
slopes_con_biomass_ln <- sapply(X = 1:dim(predTC_log_scale)[2], FUN = function(i){
  sapply(X = 1:dim(predTC_log_scale)[3], FUN = function(x){
    c(coef(lm(predTC_log_scale[,i,x]~Years))[2])
  })
})
slopes_biomass_ln <- sapply(X = 1:dim(predTT_log_scale)[2], FUN = function(i){
  sapply(X = 1:dim(predTT_log_scale)[3], FUN = function(x){
    c(coef(lm(predTT_log_scale[,i,x]~Years))[2])
  })
})
occur_trends <- data.frame(FG = levels(models[[1]]$TrData$FG) ,
                           median = apply(slopes, 2, median),
                           lower_0.25 = apply(slopes, 2, quantile, prob = 0.25),
                           upper_0.75 = apply(slopes, 2, quantile, prob = 0.75),
                           lower = apply(slopes, 2, quantile, prob = 0.025),
                           upper = apply(slopes, 2, quantile, prob = 0.975),
                           measure = 'Occurrence prob.')
biomass_cond_ln_trends <- data.frame(FG = levels(models[[1]]$TrData$FG),
                                     median = apply(slopes_con_biomass_ln, 2, median),
                                     lower_0.25 = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.25),
                                     upper_0.75 = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.75),
                                     lower = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.025),
                                     upper = apply(slopes_con_biomass_ln, 2, quantile, prob = 0.975),
                                     measure = 'Conditional cover (log)')
biomass_ln_trends <- data.frame(FG = levels(models[[1]]$TrData$FG),
                                median = apply(slopes_biomass_ln, 2, median),
                                lower_0.25 = apply(slopes_biomass_ln, 2, quantile, prob = 0.25),
                                upper_0.75 = apply(slopes_biomass_ln, 2, quantile, prob = 0.75),
                                lower = apply(slopes_biomass_ln, 2, quantile, prob = 0.025),
                                upper = apply(slopes_biomass_ln, 2, quantile, prob = 0.975),
                                measure = 'Cover (log)')
all_trends <- bind_rows(occur_trends, biomass_cond_ln_trends, biomass_ln_trends) %>%
  mutate(sig = sign(upper) == sign(lower)) %>%
  mutate(sign = ifelse(median > 0 & sig == TRUE, 'increasing', ifelse(median < 0 & sig == TRUE, 'decreasing', 'no trend'))) %>%
  mutate(sign = factor(sign, levels = c('decreasing', 'no trend', 'increasing'), ordered = TRUE)) %>%
  mutate(FG = gsub("_"," ",FG)) %>%
  mutate(FG = factor(FG, levels = gsub("_"," ",biomass_ln_trends$FG[order(biomass_ln_trends$median)]), ordered = TRUE)) %>%
  mutate(measure = factor(measure, levels = c('Occurrence prob.', 'Conditional cover (log)', 'Cover (log)')))
all_trends %>%
  ggplot(aes(x = FG, y = median, color = sign))+
  geom_vline(xintercept = c(levels(all_trends$species)[seq(3, models[[1]]$ns, by = 3)]), col = 'grey', size = 0.1, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin = lower_0.25, ymax = upper_0.75), width = 0, size = 1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
  geom_point()+
  coord_flip()+
  facet_wrap(~measure, scales = 'free_x')+
  scale_color_manual(values = c('blue', 'grey', 'red'), guide = FALSE)+
  ylab('Change per year')+
  xlab('')
# horizontal; group by FG
measure_cols <- c("black","mediumspringgreen","mediumseagreen")
measure_cols <- c("black","springgreen2","seagreen")
reps <- models[[1]]$TrData %>% 
  mutate(FG = gsub("_"," ",FG)) %>% 
  group_by(FG) %>% 
  summarize(n=length(FG)) %>% 
  mutate(median = -0.35, measure = 'Conditional cover (log)')

all_trends %>% 
  select(FG, measure, median, sign) %>% 
  arrange(FG)

# all_trends <- left_join(all_trends, select(reps,FG,n) )
# all_trends <- all_trends %>% 
#   mutate( FG = paste0(FG," (",n,")") ) %>% 
#   mutate( FG = factor(FG, 
#                       levels = c("thin turf (11)","blade (6)",
#                       "turf (14)","canopy (5)","crust (7)","animal (3)"),
#                       labels = c("thin turf (11)","blade (6)",
#                                  "turf (14)","canopy (5)","crust (7)","animal (3)")
#                       ))
# reps <- all_trends %>% 
#   select( FG, n ) %>% 
#   mutate( median = -0.35, measure = "Conditional cover (log)") %>% 
#   distinct()

hmsc_FG <- all_trends %>%
  ggplot(aes(x = FG, y = median, group = measure, col = measure))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_errorbar(aes(ymin = lower_0.25, ymax = upper_0.75), width = 0, size = 1, position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(width = 0.5))+
  geom_point(position = position_dodge(width = 0.5))+
  # geom_point(data = mutate(reps,median = median), aes(fill = FG), size = 8, col="black", shape = 21 ) +
  geom_text(data = reps, aes(label = paste0("(",n,")")), size = 2.5, col = 'black') +
  scale_color_manual(values = measure_cols) +
  scale_fill_manual(values = c("midnightblue","darkgrey","#996633","pink","red","darkred" ), guide = F) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  guides( color =  guide_legend(nrow = 2, byrow = F, reverse = F) ) +
  theme_bw() +
  theme( legend.position="top",
         panel.grid = element_blank(),
         legend.title = element_blank(),
         legend.text = element_text(size = 7),
         legend.key.size = unit(0.5, "cm"),
         legend.key = element_rect(colour = NA, fill = NA),
         legend.box.margin=margin(-10,-10,-10,-10) ) +
  ylab('Change per year (scaled)')+
  xlab('Functional Group')
hmsc_FG
ggsave( "R/Figs/hmsc_scale_change_FG.svg", plot = hmsc_FG, height = 3.5, width = 4.02)


#####
#






# add empirical means
models[[1]]$studyDesign %>% 
  group_by(year) %>% 
  summarize(ntrans=length(transect))
ogd <- bind_cols(  models[[1]]$studyDesign, comm.all )  #data.frame(m$Y)
ogd$year <- as.numeric(as.character(ogd$year))
ogd.long <- ogd.mean <- ogd %>% 
  pivot_longer( cols=names(comm.all)[1]:names(comm.all)[length(names(comm.all))], 
                names_to = "taxon", values_to = "N" )
ogd.mean <- ogd %>% 
  pivot_longer( cols=names(comm.all)[1]:names(comm.all)[length(names(comm.all))], 
                names_to = "taxon", values_to = "N" ) %>% 
  group_by( year, transect , taxon) %>% 
  summarize( N = mean(N) )
ogd.pa <- ogd %>% 
  pivot_longer( cols=names(comm.all)[1]:names(comm.all)[length(names(comm.all))], 
                names_to = "taxon", values_to = "N" ) %>% 
  group_by( transect , taxon ) %>% 
  summarize( pa = sum(N) ) %>% 
  filter( pa > 0 )

# zero presence in entire survey
ogd.zeros <- ogd %>% 
  pivot_longer( cols=names(comm.all)[1]:names(comm.all)[length(names(comm.all))], 
                names_to = "taxon", values_to = "N" ) %>% 
  group_by( taxon ) %>% 
  mutate( Nmean=mean(N) ) %>% 
  group_by( year, taxon, Nmean ) %>% 
  summarize( pa = sum(N) ) %>% 
  filter( pa == 0 ) %>% 
  arrange( -Nmean )

ogd.mean.pa <- full_join(ogd.mean, ogd.pa)
ogd.mean.pa <- ogd.mean.pa %>% 
  ungroup() %>% 
  # filter( !is.na(pa) ) %>%
  mutate( year =  as.numeric(as.character(year)) ) 
ogd.mean.pa.mean <- ogd.mean.pa %>% 
  group_by(year,taxon) %>% 
  summarize(N=mean(N))


species <- c("Fucus.distichus","Elachista.fucicola")
species <- c("Palmaria.hecatensis","Palmaria.mollis") 
species <- c("Polysiphonia","Lithothamnion.phymatodeum") 
species <- "Pyropia"
species <- c("Alaria.marginata","Hedophyllum.sessile","Egregia.menziesii")
species <- c("Alaria.marginata","Hedophyllum.sessile")
species <- c( "Mazzaella.parvula", "Mazzaella.oregona", "Mazzaella.splendens" )
species <- "Cladophora.columbiana"
species <- "Corallina"
species <- c("Colpomenia.bullosa","Colpomenia.peregrina")
species <- c("Mastocarpus","Petrocelis")
species <- c("Phyllospadix.sp.")
species <- c("Gloiopeltis.furcata")
species <- c("Barnacles","Mytilus.sp.")
species <- c("Costaria.costata","Osmundea.spectabilis","Nemalion.helminthoides")
species <- c("Scytosiphon.lomentaria","Lomentaria.hakodatensis","Salishia.firma")
species <- c("Erythrotrichia.carnea","Ectocarpus.commensalis","Elachista.fucicola")
species <- list( c("Fucus.distichus","Elachista.fucicola"),c("Palmaria.hecatensis","Palmaria.mollis"),
                 c("Polysiphonia","Lithothamnion.phymatodeum"),"Pyropia",c("Alaria.marginata","Hedophyllum.sessile"),
                 c( "Mazzaella.parvula", "Mazzaella.oregona", "Mazzaella.splendens" ),"Cladophora.columbiana",
                 "Corallina",c("Dactylosiphon.bullosus","Colpomenia.peregrina"),c("Mastocarpus","Petrocelis"),
                 c("Phyllospadix.sp."),c("Gloiopeltis.furcata"),c("Barnacles","Mytilus.sp."),
                 c("Costaria.costata","Osmundea.spectabilis","Nemalion.helminthoides"),
                 c("Scytosiphon.lomentaria","Lomentaria.hakodatensis","Salishia.firma"),
                 c("Erythrotrichia.carnea","Ectocarpus.commensalis","Elachista.fucicola"))

# species_temporal.df$taxon <- gsub("\n", ".", species_temporal.df$Species)
# predictions_abund_mean <- predictions_abund %>% 
#   group_by( year,  taxon) %>%  #site, transect,
#   summarize( N = mean(N) )

alpha_choose = 0.25
alpha_choose2 = 0.75
size_choose = 2.5
size_choose2 = 3
point_colors <- c("slateblue","firebrick", "goldenrod" )
for(i in 1:length(species)){
  abun.time.plot <- species_bmass %>% 
  filter(Species %in% species[[i]] ) %>% ungroup() #%>% 
  
  ggplot( abun.time.plot, aes( x=year, y=exp(median), group=Species, col=Species ) ) + 
  # geom_ribbon( aes(x=year, ymin=N_high/90, ymax=N_low/90), fill="grey70", alpha=0.5 ) +
  # facet_wrap(~taxon, scales="free_y",ncol=1) +
  geom_path(size=size_choose, alpha=alpha_choose2) +
  # geom_point( data=filter(ogd.mean.pa, taxon %in% species ), aes(y=N+0.01111111), size=5, alpha = 0.1 ) +
  geom_point( data=filter(ogd.mean, taxon %in% species[[i]] ), aes(y=N+0.01111111, group=taxon, col=taxon), size=size_choose2, 
              alpha = alpha_choose, position = position_dodge(width = 0.25) ) +
  geom_path( data=filter(ogd.mean.pa.mean, taxon %in% species[[i]] ), aes(y=N+0.01111111, group=taxon, col=taxon), lwd=0.75, alpha = 1 ) +
  # stat_summary( data=filter(ogd.mean.pa, taxon %in% species ), aes(y=N), fun = "median", geom="line", size = 0.5 ) +
  scale_color_manual(values=point_colors[1:length(species[[i]])]) +
  ylab("Percent Cover") +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank()) +
  # ylim( c(0,18) ) +
  # scale_y_sqrt(labels = comma, breaks = c(0, 1.01111111, 10, 25, 50, 75, 100)) +
  scale_y_log10(labels = comma) +
  # scale_y_continuous(labels = comma) +
  guides(color=guide_legend(nrow=3,ncol=1,byrow=F))  
  
  ggsave( paste0("R/Figs/temporal_trends_select_taxa/",paste0(species[[i]],collapse="_"),"_model+data.svg"), 
        width=3,height=3.5)
}
#

















predictions_pab_trait <- left_join( predictions_pab, taxon.key )
predictions_copp_trait <- left_join( predictions_copp, taxon.key )
predictions_abund_trait <- left_join( predictions_abund, taxon.key )


# only pull the 10 most common taxa
taxa <- colnames(models[[1]]$Y)
top6 <- taxa[c(1:6,9,14,15)]

custom_occur <- taxa[c(4,33,29,18,22,30,34,3,34,40,15,11,5,38)]
custom_cop   <- taxa[c(33,25,14,36,7,1,13,8,39,27,42,2,37,15)]
custom_abun  <- taxa[c(33,4,29,18,25,14,36,30,7,1,13,27,45,22,15,38)]
# 10 rarest taxa
bot10 <- taxa[(length(taxa)-8):length(taxa)]

taxa2plot <- customXY

# windows(6,4)
ggplot( filter(predictions_pab,taxon %in% custom_occur & year %in% c(2012,2019)), 
        aes(x = elev, y = N,
            fill=factor(year), col=factor(year) ))+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, col="gray75")+
  facet_wrap(~taxon, scales = "free_y", ncol=4)+
  theme_classic() +
  scale_fill_manual(values=c("gold1", "darkslategray4"), name = "Year")+
  geom_point( data = filter( comm_final_pa, taxon %in% custom_occur, year %in% c(2012,2019)), 
              pch=21, alpha = 0.75 ) +
  geom_line(size = 0.5 ) +
  scale_color_manual(values = c("darkslategray","black"), name = "Year") +
  ylab("Probability of occurrence") + xlab("Shore height (cm)") +
  theme( strip.text.x = element_text(size = 7)) +
  coord_cartesian(ylim=c(0,1))
ggsave("R/Figs/hmsc_response_curves_hurdle_occur.svg", width = 7, height = 5)
ggplot( filter(predictions_copp,taxon %in% custom_cop & year %in% c(2012,2019)), 
        aes(x = elev, y = N,
            fill=factor(year), col=factor(year) ))+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, col="gray75")+
  facet_wrap(~taxon, scales = "free_y", ncol=4)+
  theme_classic() +
  scale_fill_manual(values=c("gold1", "darkslategray4"), name = "Year")+
  geom_point( data = filter( comm_final_cond, taxon %in% custom_cop, year %in% c(2012,2019)), 
              pch=21, alpha = 0.75 ) +
  geom_line(size = 0.5 ) +
  scale_color_manual(values = c("darkslategray","black"), name = "Year") +
  scale_y_sqrt(breaks=c(0,1,10,25,50,100)) +
  ylab("Cover conditional on presence (%)") + xlab("Shore height (cm)") +
  theme( strip.text.x = element_text(size = 7)) +
  coord_cartesian(ylim=c(0,100))
ggsave("R/Figs/hmsc_response_curves_hurdle_cop.svg", width = 7, height = 5)
ggplot( filter(predictions_abund,taxon %in% custom_abun & year %in% c(2012,2019)), 
        aes(x = elev, y = N,
            fill=factor(year), col=factor(year) ))+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, col="gray75")+
  facet_wrap(~taxon, scales = "free_y", ncol=4)+
  theme_classic() +
  scale_fill_manual(values=c("gold1", "darkslategray4"), name = "Year")+
  geom_point( data = filter( comm_final, taxon %in% custom_abun, year %in% c(2012,2019)), 
              pch=21, alpha = 0.75 ) +
  geom_line(size = 0.5 ) +
  scale_color_manual(values = c("darkslategray","black"), name = "Year") +
  scale_y_sqrt(breaks=c(0,1,10,25,50,100)) +
  ylab("Cover (%)") + xlab("Shore height (cm)") +
  theme( strip.text.x = element_text(size = 7)) +
  coord_cartesian(ylim=c(0,100))
ggsave("R/Figs/hmsc_response_curves_hurdle_total.svg", width = 7, height = 5)

# just show Fucus
fuc <- "Fucus.distichus"
fuc <- "Mytilus.sp."
fuc <- "Hedophyllum.sessile"
# windows(5,4)
ggplot( filter(predictions_copp,taxon %in% fuc ), 
        aes(x = elev, y = N ))+
  facet_wrap(~year)+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, fill="gray75")+
  geom_line(aes(group=year),size = 0.5)+
  theme_classic() +
  ggtitle(fuc) +
  geom_point( data = filter( comm_final, taxon %in% fuc), pch=21 ) +
  scale_y_sqrt(breaks=c(1,10,50,100,200)) +
  ylab("Percent cover") + xlab("Shore height (cm)") +
  coord_cartesian(ylim = c(-0, 100))
ggsave(paste0("R/Figs/model+data_",fuc,".svg"),width=5,height=5 )
#




## get peak elevation and "abundance" for each species in YEAR in each run #####
# elevation peaks for each species in each YEAR in each run
peaks  <- lapply( predY_pa, 
                  function(i) lapply( split( i, newDF$year ), 
                                      function(l) apply(matrix(l,byrow = F,ncol = ncol(models[[1]]$Y)), 2, 
                                                        function(z) mean(unique(newDF$elev)[which(z==max(z))])) ) )
peaks_bind <- lapply( peaks, function(z) do.call(rbind,z) )
peaks_array <- abind::abind(peaks_bind, along=3)

diff_half_prob <- lapply( predY_pa, 
                     function(i) lapply( split( i, newDF$year ), 
                                         function(l) apply(matrix(l,byrow = F,ncol = ncol(models[[1]]$Y)), 2, 
                                                           function(z) diff(range(unique(newDF$elev)[which(z>=max(z)/2)]))) ) )
diffp_bind <- lapply( diff_half_prob, function(z) do.call(rbind,z) )
diffp_array <- abind::abind(diffp_bind, along=3)

peak_half_prob <- lapply( predY_pa, 
                     function(i) lapply( split( i, newDF$year ), 
                                         function(l) apply(matrix(l,byrow = F,ncol = ncol(models[[1]]$Y)), 2, 
                                                           function(z) mean(range(unique(newDF$elev)[which(z>=max(z)/2)]))) ) )
halfp_bind <- lapply( peak_half_prob, function(z) do.call(rbind,z) )
halfp_array <- abind::abind(halfp_bind, along=3)
# ALSO THINK ABOUT CHANGES IN UPPER VERSUS LOWER DIST. (max or min at 50% of max occurrence prob)
min_half_prob <- lapply( predY_pa, 
                          function(i) lapply( split( i, newDF$year ), 
                                              function(l) apply(matrix(l,byrow = F,ncol = ncol(models[[1]]$Y)), 2, 
                                                                function(z) min(unique(newDF$elev)[which(z>=max(z)/2)])) ) )
minp_bind <- lapply( min_half_prob, function(z) do.call(rbind,z) )
minp_array <- abind::abind(minp_bind, along=3)
max_half_prob <- lapply( predY_pa, 
                          function(i) lapply( split( i, newDF$year ), 
                                              function(l) apply(matrix(l,byrow = F,ncol = ncol(models[[1]]$Y)), 2, 
                                                                function(z) max(unique(newDF$elev)[which(z>=max(z)/2)])) ) )
maxp_bind <- lapply( max_half_prob, function(z) do.call(rbind,z) )
maxp_array <- abind::abind(maxp_bind, along=3)

array_use <- halfp_array

## calculate shifts as slope of peaks over time
Years <- 2012:2019
dim(array_use)
slopes_peak <- sapply(X = 1:dim(array_use)[2], FUN = function(i){
  sapply(X = 1:dim(array_use)[3], FUN = function(x){
    c(coef(lm(array_use[,i,x]~Years))[2])
  })
})
# change in elevational breadth of the distribution
# slopes_range <- sapply(X = 1:dim(diffp_array)[2], FUN = function(i){
#   sapply(X = 1:dim(diffp_array)[3], FUN = function(x){
#     c(coef(lm(diffp_array[,i,x]~Years))[2])
#   })
# })
# slopes.range <- slopes_range
# range.init           <- apply( diffp_array, c(2,3), function(z) z[1] )
# range.init.med       <- apply(range.init, 1, quantile, prob = 0.5, na.rm = TRUE)
# range.shifts.med     <- apply(slopes.range*8, 2, quantile, prob = 0.5, na.rm = TRUE)
# range.shifts.high    <- apply(slopes.range*8, 2, quantile, prob = 0.025, na.rm = TRUE)
# range.shifts.low     <- apply(slopes.range*8, 2, quantile, prob = 0.975, na.rm = TRUE)
# range.shifts.25     <- apply(slopes.range*8, 2, quantile, prob = 0.25, na.rm = TRUE)
# range.shifts.75     <- apply(slopes.range*8, 2, quantile, prob = 0.75, na.rm = TRUE)
# range.shifts.summary <- data.frame(range.init.med, range.shifts.med, range.shifts.low, range.shifts.high,range.shifts.25,range.shifts.75)
# range.shifts.summary.algae <- range.shifts.summary[taxon.key$funct != "animal",]
# range.shift <- data.frame( shift = c(slopes.range*8) )
# hist(c(slopes.range*8))
# boxplot(c(slopes.range*8), notch = T); abline(h = 0, col = 'red' )
# mean(c(slopes.range*8))
# quantile(c(slopes.range*8),probs = 0.5)
# summary( lm(c(slopes.range*8)~1) )
# slopes_min <- sapply(X = 1:dim(minp_array)[2], FUN = function(i){
#   sapply(X = 1:dim(minp_array)[3], FUN = function(x){
#     c(coef(lm(minp_array[,i,x]~Years))[2])
#   })
# })
# slopes.min <- slopes_min
# min.init           <- apply( minp_array, c(2,3), function(z) z[1] )
# min.init.med       <- apply(min.init, 1, quantile, prob = 0.5, na.rm = TRUE)
# min.shifts.med     <- apply(slopes.min*8, 2, quantile, prob = 0.5, na.rm = TRUE)
# min.shifts.high    <- apply(slopes.min*8, 2, quantile, prob = 0.025, na.rm = TRUE)
# min.shifts.low     <- apply(slopes.min*8, 2, quantile, prob = 0.975, na.rm = TRUE)
# min.shifts.25     <- apply(slopes.min*8, 2, quantile, prob = 0.25, na.rm = TRUE)
# min.shifts.75     <- apply(slopes.min*8, 2, quantile, prob = 0.75, na.rm = TRUE)
# min.shifts.summary <- data.frame(min.init.med, min.shifts.med, min.shifts.low, min.shifts.high,min.shifts.25,min.shifts.75)
# min.shifts.summary.algae <- min.shifts.summary[taxon.key$funct != "animal",]
# min.shift <- data.frame( shift = c(slopes.min*8) )
# hist(c(slopes.min*8))
# boxplot(c(slopes.min*8), notch = T); abline(h = 0, col = 'red' )
# mean(c(slopes.min*8))
# quantile(c(slopes.min*8),probs = 0.5)
# summary( lm(c(slopes.min*8)~1) )
# slopes_max <- sapply(X = 1:dim(maxp_array)[2], FUN = function(i){
#   sapply(X = 1:dim(maxp_array)[3], FUN = function(x){
#     c(coef(lm(maxp_array[,i,x]~Years))[2])
#   })
# })
# slopes.max <- slopes_max
# max.init           <- apply( maxp_array, c(2,3), function(z) z[1] )
# max.init.med       <- apply(max.init, 1, quantile, prob = 0.5, na.rm = TRUE)
# max.shifts.med     <- apply(slopes.max*8, 2, quantile, prob = 0.5, na.rm = TRUE)
# max.shifts.high    <- apply(slopes.max*8, 2, quantile, prob = 0.025, na.rm = TRUE)
# max.shifts.low     <- apply(slopes.max*8, 2, quantile, prob = 0.975, na.rm = TRUE)
# max.shifts.25     <- apply(slopes.max*8, 2, quantile, prob = 0.25, na.rm = TRUE)
# max.shifts.75     <- apply(slopes.max*8, 2, quantile, prob = 0.75, na.rm = TRUE)
# max.shifts.summary <- data.frame(max.init.med, max.shifts.med, max.shifts.low, max.shifts.high,max.shifts.25,max.shifts.75)
# max.shifts.summary.algae <- max.shifts.summary[taxon.key$funct != "animal",]
# max.shift <- data.frame( shift = c(slopes.max*8) )
# hist(c(slopes.max*8))
# boxplot(c(slopes.max*8), notch = T); abline(h = 0, col = 'red' )
# mean(c(slopes.max*8))
# quantile(c(slopes.max*8),probs = 0.5)
# summary( lm(c(slopes.max*8)~1) )

# ## calculate differences between 2012 and 2019 to get shift in end member states
# # elevation peak
# elev.shifts.run     <- apply( peaks_array, c(2,3), function(z) z[length(z)]-z[1] )
slopes.peak.algae     <- slopes_peak[,taxon.key$funct != "animal"]
slopes.peak     <- slopes_peak
elev.init           <- apply( array_use, c(2,3), function(z) z[1] )
# elev.init.algae           <- elev.init[taxon.key$funct != "animal",]
# elev.init.algae           <- elev.init
elev.init.med       <- apply(elev.init, 1, quantile, prob = 0.5, na.rm = TRUE)
elev.shifts.med     <- apply(slopes.peak*8, 2, quantile, prob = 0.5, na.rm = TRUE)
elev.shifts.high    <- apply(slopes.peak*8, 2, quantile, prob = 0.025, na.rm = TRUE)
elev.shifts.low     <- apply(slopes.peak*8, 2, quantile, prob = 0.975, na.rm = TRUE)
elev.shifts.25     <- apply(slopes.peak*8, 2, quantile, prob = 0.25, na.rm = TRUE)
elev.shifts.75     <- apply(slopes.peak*8, 2, quantile, prob = 0.75, na.rm = TRUE)
elev.shifts.summary <- data.frame(elev.init.med, elev.shifts.med, elev.shifts.low, elev.shifts.high,elev.shifts.25,elev.shifts.75)
elev.shifts.summary.algae <- elev.shifts.summary[taxon.key$funct != "animal",]
elev.shift <- data.frame( shift = c(slopes.peak*8) )
hist(c(slopes.peak*8))
boxplot(c(slopes.peak*8), notch = T); abline(h = 0, col = 'red' )
mean(c(slopes.peak*8))
quantile(c(slopes.peak*8),probs = 0.5)
summary( lm(c(slopes.peak*8)~1) )

quantile( c(slopes.peak.algae*8), prob = rev(c(0.025,0.25,0.5,0.75,0.975)) )
# library(easystats)
# library(rstanarm)
# stan1 <- stan_glm(shift ~ 1, data = elev.shift)
# posteriors <- describe_posterior(stan1)
# # for a nicer table
# print_md(posteriors, digits = 2)
# plot(stan1)

# abundance
# cover shifts
biomass_ln_slopes
# abun.shifts.run     <- apply( abunds_array, c(2,3), function(z) z[length(z)]/z[1] )
# abun.shifts.run.algae     <- biomass_ln_slopes[,taxon.key$funct != "animal"]
abun.shifts.run     <- biomass_ln_slopes
## use cover slopes for this part
abun.init <- apply( log_biomass, c(2,3), function(z) exp(z[1]) )
# abun.init           <- apply( abunds_array, c(2,3), function(z) z[1] )
# abun.init.algae           <- abun.init[taxon.key$funct != "animal",]
abun.init.algae           <- abun.init
abun.init.med       <- apply(abun.init.algae, 1, quantile, prob = 0.5, na.rm = TRUE)
abun.shifts.run <- exp(t(log(abun.init.algae))+(abun.shifts.run*8)) / t(abun.init.algae)
abun.shifts.run.algae     <- abun.shifts.run[,taxon.key$funct != "animal"]
abun.shifts.med     <- apply(abun.shifts.run, 2, quantile, prob = 0.5, na.rm = TRUE)
abun.shifts.high    <- apply(abun.shifts.run, 2, quantile, prob = 0.025, na.rm = TRUE)
abun.shifts.low     <- apply(abun.shifts.run, 2, quantile, prob = 0.975, na.rm = TRUE)
abun.shifts.25     <- apply(abun.shifts.run, 2, quantile, prob = 0.25, na.rm = TRUE)
abun.shifts.75     <- apply(abun.shifts.run, 2, quantile, prob = 0.75, na.rm = TRUE)
abun.shifts.summary <- data.frame(abun.init.med,abun.shifts.med, abun.shifts.low, abun.shifts.high,abun.shifts.25,abun.shifts.75)
abun.shifts.summary.algae <- abun.shifts.summary[taxon.key$funct != "animal",]
hist( log(c(abun.shifts.run), base=2) )
hist( log(c(abun.shifts.run.algae), base=2) )
summary( lm(log(c(abun.shifts.run.algae),base=2)~1) )

2^quantile( log(abun.shifts.run.algae,base=2), prob = (c(0.025,0.25,0.5,0.75,0.975)) )-1


## combine these results
shift.summary <- data.frame( elev.shifts.summary, abun.shifts.summary) #, range.shifts.summary, min.shifts.summary, max.shifts.summary )
shift.summary$taxon <- colnames(models[[1]]$Y)#[taxon.key$funct != "animal"]
shift.summary %>% arrange(-abun.init.med)
shift.summary$rank <- 1:46
ggplot( shift.summary, aes(x = abun.shifts.med, y = elev.shifts.med) ) +
  # geom_hline( yintercept=0, lty=2 ) + geom_vline( xintercept = 0, lty=2 ) +
  geom_linerange( aes( ymin = elev.shifts.low, ymax= elev.shifts.high), alpha=0.25) +
  geom_linerange( aes(xmin = abun.shifts.low, 
                      xmax = abun.shifts.high), alpha=0.25 ) +
  geom_point() + 
  scale_x_continuous( trans = "log2" ) +
  # scale_x_continuous( breaks=c(4,2,0,-2,-4), 
                      # labels=c('16x','4x','0','1/4x','1/16x')) +
  ylab( "Elevation shift (cm)" ) + xlab( "Abundance shift" ) +
  theme_bw() + theme( panel.grid.minor = element_blank() )
# ggsave( "R/Figs/shifts_error.pdf", width=4, height=4 )
# ggplot( shift.summary, aes(x = range.shifts.med, y = elev.shifts.med) ) + 
#   geom_hline(yintercept = 0) + geom_vline( xintercept = 0 ) +
#   geom_point() + geom_smooth()
# ggplot( shift.summary, aes(x = min.shifts.med, y = elev.shifts.med) ) + 
#   geom_hline(yintercept = 0) + geom_vline( xintercept = 0 ) +
#   geom_point() + geom_smooth()
# ggplot( shift.summary, aes(x = max.shifts.med, y = elev.shifts.med) ) + 
#   geom_hline(yintercept = 0) + geom_vline( xintercept = 0 ) +
#   geom_point() + geom_smooth()
# ggplot( shift.summary, aes(x = max.shifts.med, y = min.shifts.med) ) + 
#   geom_hline(yintercept = 0) + geom_vline( xintercept = 0 ) +
#   geom_point() + geom_smooth()
# psych::pairs.panels( shift.summary %>% select( elev.shifts.med, range.shifts.med, min.shifts.med, max.shifts.med) )

# correlation of median responses
with( shift.summary, cor.test( elev.shifts.med, abun.shifts.med, method = 'spearman' ) )
with( shift.summary, cor.test( elev.shifts.med, abun.shifts.med, method = 'pearson' ) )
filter( shift.summary, taxon=="Fucus.distichus" )

# "significant" shifts as those that did not include zero
shift.summary %>% 
  mutate( peak.sig = elev.shifts.high*elev.shifts.low > 0 ) %>% 
  mutate( abun.high.one = abun.shifts.high > 1,
          abun.low.one = abun.shifts.low > 1 ) %>% 
  mutate( abun.sig = abun.high.one == abun.low.one ) %>% 
  select( rank, taxon, peak.sig, abun.sig, elev.shifts.med, abun.shifts.med ) %>% 
  filter( peak.sig == TRUE | abun.sig == TRUE )

# individual shifts
shift.summary <- left_join(shift.summary, taxon.key)
shift.summary <- shift.summary %>% 
  mutate(order_FG = factor(funct, levels = biomass_ln_trends$FG[order(biomass_ln_trends$median)], ordered = TRUE)) %>% 
  mutate(order_taxa = factor(taxon, levels = cover_trends_taxa$species[order(cover_trends_taxa$median)], ordered = TRUE)) %>% 
  mutate(rank = 1:nrow(shift.summary)) %>% 
  arrange(order_FG, order_taxa) %>% 
  mutate(newrank = 1:nrow(shift.summary) )
shift.summary$funct <- factor( shift.summary$funct, levels =  )
data.frame( x = seq(3,46,3))

shift.summary.algae <- shift.summary %>% filter(funct != "animal")
  
elev.shift.plot <- ggplot( shift.summary, aes(x=newrank,
                                                    y=elev.shifts.med, col = order_FG)) + 
  # geom_vline( xintercept = seq(4,46,5), lty = 2, col = "gray" ) +
  geom_hline( yintercept=0 ) +
  geom_errorbar( aes(ymin=elev.shifts.low,ymax=elev.shifts.high), width = 0 ) +
  geom_errorbar( aes(ymin=elev.shifts.25,ymax=elev.shifts.75), width = 0, lwd = 1 ) +
  geom_point() +
  ylab( "Peak elevation shift (cm)" ) + xlab("") +
  scale_x_continuous(breaks = seq(1,46,1), labels = shift.summary$taxon ) + #c(1,10,20,30,40,46)) +
  scale_y_continuous(breaks = seq(-800,400,100) ) + #c(1,10,20,30,40,46)) +
  scale_color_manual(values = c("red","darkgrey","darkred","#996633","pink","midnightblue") ) +
  theme( panel.grid.minor.x = element_blank(),
         legend.position = "none",
         panel.background = element_rect(fill = "whitesmoke",
                                         colour = "whitesmoke",
                                         size = 0.5, linetype = "solid")  ) +
  coord_cartesian( xlim = c(1,47)) +
  coord_flip()
abun.shift.plot <- ggplot( shift.summary, aes(x=newrank,
                                                    y=log(abun.shifts.med,base=2), col = order_FG) ) + 
  # geom_vline( xintercept = seq(4,46,5), lty = 2, col = "gray" ) +
  geom_hline( yintercept=0 ) +
  geom_errorbar( aes(ymin=log(abun.shifts.low,base=2),ymax=log(abun.shifts.high,base=2)), width = 0 ) +
  geom_errorbar( aes(ymin=log(abun.shifts.25,base=2),ymax=log(abun.shifts.75,base=2)), width = 0, lwd = 1 ) +
  geom_point() +
  ylab( "Percent cover shift" ) + xlab("") +
  scale_y_continuous( breaks=c(log(64,base=2),log(16,base=2),log(4,base=2),0,
                               log(1/4,base=2),log(1/16,base=2),log(1/64,base=2)),
                      labels=c('64x','16x','4x','0','1/4x','1/16x','1/64x') ) +
  scale_x_continuous( breaks = seq(1,46,1), labels = NULL ) +
  scale_color_manual(values = c("red","darkgrey","darkred","#996633","pink","midnightblue") ) +
  theme( panel.grid.minor.x = element_blank(),
         legend.position = "none",
         panel.background = element_rect(fill = "whitesmoke",
                                         colour = "whitesmoke",
                                         size = 0.5, linetype = "solid")  ) +
  theme( axis.text.y = element_blank(),
    axis.ticks.y  = element_blank()) +
  coord_flip()
plot_grid( elev.shift.plot, abun.shift.plot, ncol=2, rel_widths = c(1.75,1) )
ggsave( "R/Figs/shifts_2panel.svg", width=7, height=6 )

# # Find the predicted peak for each instance
# # filter out animals
# peaks <- predictions_pab_trait %>%
#   filter( funct != "animals" ) %>% 
#   group_by( year, taxon, funct ) %>%
#   summarize( peak = mean(elev[which(N==max(N))]), peaksd = mean(elev[which(N==max(N))],na.rm=T) )
# peaks$peak
# ggplot( filter(peaks, taxon %in% custom_occur), aes(x=year,y=peak) ) + facet_wrap(~taxon) +
#   geom_point()
# ggplot( peaks, aes(x=year,y=peak) ) + #facet_wrap(~taxon) +
#   # geom_path(aes(group=taxon), alpha=0.5) + 
#   geom_smooth(method="lm") + geom_smooth(aes(group=taxon),col='black',se=F,lwd=0.5,alpha=0.5)
# summary(lm(peak~1+year,data=peaks))
# 
# # get difference between peaks for 2012 and 2019
# peak_shift <- peaks %>% 
#   group_by( taxon, funct ) %>% 
#   filter( year %in% c(2012,2019) ) %>% 
#   summarize( shift = diff(peak))
# 
# # merge 2012 peaks with peak shift to compare shift relative to starting point
# peak_initial <- peaks %>% filter( year==2012 )
# peak_compare <- left_join( peak_shift, peak_initial )
# 
# cutoff <- 20
# peak_compare %>% filter(shift <= -cutoff)
# peak_compare %>% filter(shift >= cutoff)
# peak_compare %>% filter(shift < cutoff & shift > -cutoff)
# peak_compare %>% filter(shift == 0 )
# peak_compare %>% arrange(shift)
# peak_compare %>% arrange(-shift)
# 
# # plot by functional group
# trait.p <- ggplot( peak_shift, aes(x=reorder(funct, shift, FUN = median), y=shift/100) ) + 
#   geom_hline (yintercept=0, lty=2 ) +
#   geom_boxplot()  + geom_point() +
#   xlab("Functional group") + ylab("Peak shift (meters)") +
#   # scale_fill_manual(values=c("whitesmoke","dodgerblue")) +
#   theme_classic() +
#   theme( axis.text.x = element_text(angle=45,hjust=1,vjust=1) ) +
#   theme(legend.position = "none")
# 
# peak_shift %>% arrange(-shift)
# peak_shift %>% arrange(shift)
# peak_shift %>% filter(shift==0)
# peak_shift %>% filter(shift >-10 & shift < 10)
# shift_increase <- peak_shift %>% arrange(-shift)
# choose <- shift_increase$taxon[1:6]
# ggplot( filter(peaks, taxon %in% choose), aes(x=year,y=peak) ) + facet_wrap(~taxon) +
#   geom_point()






# need to add lines showing the realm of possible shifts
xs <- range( models[[1]]$XData$elev )
y1 <- c(0,diff(xs))
y2 <- c(-diff(xs),0)
df.bound <- data.frame( x1=xs[1],x2=xs[2],y1,y2 )
df.poly <- data.frame( x=rep(xs,each=2), y=c(0,diff(xs),0,-diff(xs)) )



# FOUR panels to compare peak and abundance shifts
summary(lm( elev.shifts.med ~ elev.init.med, data=shift.summary.algae ))
cor.test( x=shift.summary.algae$elev.shifts.med, y = shift.summary.algae$elev.init.med )
(a <- ggplot( shift.summary.algae, aes(x=elev.init.med,y=elev.shifts.med)) +
    geom_polygon( data=df.poly, aes(x=x,y=y), fill='whitesmoke', col='slategray', lty=2) +
    geom_hline( yintercept = 0, lty=2 ) +
    geom_smooth(method='lm', se=T, col='black') +
    geom_point( size=3, pch=1, col='slateblue' ) +
    ylab("Peak elevationshift (cm)") + xlab("Initial peak elevation (cm)") +
    theme_classic() +
    coord_cartesian(ylim = c(-125,20),xlim = c(60,380)))


# shift in abundance ~ initial peak elevation
(b <- ggplot( shift.summary.algae, aes(x=elev.init.med,y=log(abun.shifts.med,base = 2))) + 
    geom_hline( yintercept=1, lty=2 ) +
    # geom_smooth(method='lm', se=T, col='black') +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("Cover shift") + xlab("Initial peak elevation (cm)") +
    scale_y_continuous( breaks=c(log(4,base = 2),log(2,base=2),0,log(1/2,base=2),log(1/4, base = 2),log(1/8, base = 2),log(1/16, base = 2)),
                        labels=c('4x','2x','0','1/2x','1/4x','1/8x','1/16x')) +
    theme_classic()  +
  coord_cartesian(xlim = c(60,380)) )

# peak shift ~ initial abundance
(c <- ggplot( shift.summary.algae, aes(x=abun.init.med,y=elev.shifts.med)) + 
    geom_hline( yintercept=0, lty=2 ) +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("Peak elevation shift (cm)") + xlab("Initial cover (%)") +
    scale_x_continuous(trans = "log2") +
    theme_classic() + 
    coord_cartesian(ylim = c(-125,20)) )

# abund shift ~ initial abundance
summary(lm( log(abun.shifts.med, base = 2)~log(abun.init.med,base=2), shift.summary.algae ))
cor.test( x=log(shift.summary.algae$abun.shifts.med, base = 2), y = log(shift.summary.algae$abun.init.med,base=2) )
(d <- ggplot( shift.summary.algae, aes(x=abun.init.med,y=log(abun.shifts.med, base = 2))) + 
    geom_hline( yintercept=0, lty=2 ) +
    geom_smooth(method = 'lm', se = T, col='black' ) +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("Cover shift") + xlab("Initial cover (%)") +
    scale_y_continuous( breaks=c(log(4,base = 2),log(2,base=2),0,log(1/2,base=2),log(1/4, base = 2),log(1/8, base = 2),log(1/16, base = 2)),
                        labels=c('4x','2x','0','1/2x','1/4x','1/8x','1/16x')) +
    scale_x_continuous(trans = "log2") +
    theme_classic() )

cowplot::plot_grid( a, c, b, d, ncol=2, align = 'hv', labels = "AUTO" )
ggsave(file="R/Figs/abundance+peak_shift_intial.svg",width = 6, height = 6)

head(shift.summary %>% arrange(-elev.init.med))
head(shift.summary %>% arrange(elev.shifts.med),20)
head(shift.summary %>% arrange(-abun.init.med))
head(shift.summary %>% arrange(-abun.shifts.med))
head(shift.summary %>% filter(elev.init.med < 200) %>% arrange, 20)

noxshift <- compare_all$taxon[ compare_all$shift.x == 0 ]

# basically no relationship between peak and abundance shift
ggplot( compare_all, aes(x=log(shift.y,base=2),y=shift.x) ) + 
  geom_point() + geom_smooth( method='lm' )


with( compare_all, cor.test( shift.x,shift.y) )
with( compare_all, cor.test( shift.x,log(shift.y,base=2)) )


# boxplots
ylimits1 <- c(-max(abs(range(shift.summary$elev.shifts.med))),max(abs(range(shift.summary$elev.shifts.med))))
ylimits2 <- c(-6,6) #c(2^-max(abs(log(range(compare_all$shift.y),base=2))), 2^max(abs(log(range(compare_all$shift.y),base=2))))
# (bpa <- ggplot( shift.summary, aes(y=elev.shifts.med,x=1) ) +
#     # geom_boxplot(notch = T,outlier.color = NA) + 
#     geom_violin(outlier.color = NA, draw_quantiles = 0.5, trim=FALSE ) +
#     geom_hline(yintercept=0,lty=2)+
#     geom_point(alpha=0.25) +
#     ylab( "Peak elevation shift (cm)" )   +
#     scale_y_continuous(limits=ylimits1,breaks=c(-300,-200,-100,-50,0,50,100,200,300),
#                        position="left") +
#     theme_classic() +
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) )
# (bpb <- ggplot( shift.summary, aes(y=log(abun.shifts.med,base=2),x=1) ) +
#     # geom_boxplot(notch=TRUE) + 
#     geom_violin(outlier.color = NA, draw_quantiles = 0.5, trim=FALSE ) +
#     geom_hline(yintercept=0,lty=2)+
#     geom_point(alpha=0.25) +
#     ylab( "Abundance shift" ) +
#     scale_y_continuous(limits=ylimits2,breaks=c(log(50,base=2),log(10,base=2),log(2,base=2),0,log(0.5,base=2),log(1/10,base=2),log(1/50,base=2)),
#                        labels=c('50x','10x','2x','0','1/2x','1/10x','1/50x'),
#                        position="right") +
#     theme_classic() +
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) )
# # windows(3,2.5)
# cowplot::plot_grid(bpa,bpb, align = "h", axis='tblr', labels = "AUTO")
# #

# add nice scatteplot
# taxa to plot
compare_all_plot <- shift.summary
compare_all_plot$labels <- factor( compare_all_plot$taxon, labels=1:46 )
# taxalabel <- top6
# taxalabel <- customXY
# # color points by kingdom (red, green, brown)
# compare_all_plot$first <- unlist(lapply( strsplit(compare_all_plot$taxon, split="[.]|_"), function(z) z[1] ))
# tax <- read_csv( "Data/taxa/algal_taxonomy.csv" )
# compare_all_plot <- left_join( compare_all_plot, select(tax,first=query,phylum) )
# # fill in gaps
# compare_all_plot$phylum[ compare_all_plot$first == "Ralfsioid"] <-  "Ochrophyta"
# compare_all_plot$phylum[ compare_all_plot$first == "Petrocelis"] <-  "Rhodophyta"
# compare_all_plot$phylum[ compare_all_plot$first == "coralline"] <-  "Rhodophyta"
# compare_all_plot$phylum[ is.na(compare_all_plot$phylum) ] <-  "animal"
# compare_all_plot$group <- compare_all_plot$phylum
# deleted code to label points with taxon names
# geom_text_repel(data=filter(compare_all,taxon%in%taxalabel),aes(label=taxon),
#                 box.padding = 0.3, point.padding = 0.1, size=3) +
# color points by function group
compare_all_plot_fun <- left_join(compare_all_plot, taxon.key)
compare_all_plot_fun$group <- compare_all_plot_fun$funct_Sep2020
compare_all_plot_fun$group[ is.na(compare_all_plot_fun$group) ] <- 'animal'
compare_all_plot_fun$group <- factor( compare_all_plot_fun$group, 
                                      levels = rev(c('canopy','blade','crust','thin_turf','turf','animal')))
compare_all_plot_fun <- compare_all_plot_fun %>% 
  arrange(group)
  
# functional group colors
c("darkred", "red","pink", "darkgrey", "#996633","whitesmoke")

cor( elev.shift.all.df$elev.shift, log(abun.shift.all.df$abun.shift,base = 2) )
plot( elev.shift.all.df$elev.shift, log(abun.shift.all.df$abun.shift,base = 2) )
cor( elev.shift.all.df$elev.shift, log(abun.shift.all.df$abun.shift,base = 2), method='spearman' )
cor( elev.shifts.med, log(abun.shifts.med, base = 2) )
plot( elev.shifts.med, log(abun.shifts.med, base = 2) )
cor( elev.shifts.med, log(abun.shifts.med, base = 2), method = 'spearman' )
compare_all_plot_fun %>%  arrange(elev.shifts.med)
compare_all_plot_algae <- filter( compare_all_plot_fun, funct != "animal" )
compare_all_plot_algae$group <- factor(as.character(compare_all_plot_algae$funct), levels = c('turf','thin_turf','crust','blade','canopy'))

windows(3,3)
# use full set here
psych::pairs.panels( data.frame(elevation = elev.shift.all.df$elev.shift, 
                                cover = log(abun.shift.all.df$abun.shift,base = 2)),
                     method = "pearson", smoother = TRUE,
                     hist.col = "whitesmoke", density = FALSE, rug = FALSE )
dev.off()
xy <- ggplot( compare_all_plot_algae, aes(x=log(abun.shifts.med,base = 2),y=elev.shifts.med)) + 
    geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
    # geom_hline(yintercept = quantile( c(slopes.peak.algae*8), prob = 0.5 ), lty = 1, col = "lightgray", lwd = 0.33)+
    # geom_vline(xintercept = quantile( log(abun.shifts.run.algae,base=2), prob = 0.5 ), lty = 1, col = "seagreen", lwd = 0.33)+
    # data = abun.shift.all.df.summary, aes(y = ybar, x = median )
    geom_point( pch = 21, size = 1.5, fill = "whitesmoke" ) + #aes(fill=group)
    theme_classic() +
    scale_x_continuous(breaks=c(log(16,base=2),log(4,base=2),log(2,base=2),0,log(0.5,base=2),
                                log(1/8,base=2),log(1/32,base=2)),
                       labels=c('16x','4x','2x','0','1/2x','1/8x','1/32x'),
                       # labels=c('50x','10x','4x','2x','0','0.5x','0.25x','0.0125x','0.0625x','0.05x'),
                       position="bottom") +
    # scale_x_continuous(breaks=c(log(50,base=2),log(10,base=2),log(5,base=2),log(2,base=2),0,log(0.5,base=2),
    #                             log(1/5,base=2),log(1/10,base=2),log(1/50,base=2)),
    #                    labels=c('50x','10x','5x','2x','0','0.5x','0.2x','0.1x','0.02x'),
    #                    position="bottom") +
    # ylim( c(-150,30) ) +
    scale_fill_manual( values=(c("darkred", "red","pink", "darkgrey", "#996633")), guide='none' ) +
    annotate("text", label = expression(paste(italic("r "),"= -0.05")), x = log(3,base=2), y = 100, hjust = 0) +
    xlab("Cover shift") + ylab("Elevation shift (cm)") +
    coord_cartesian( xlim = c(log(1/32,base=2), log(32,base=2)),ylim = c(-125,125))
xy

# densities of all shifts
elev.shift.all.df <- data.frame( elev.shift = c(slopes_peak[,taxon.key$funct != "animal"]*8) )
hist(elev.shift.all.df$elev.shift, freq = F, las = 1)
summary(elev.shift.all.df$elev.shift)
elev.shift.all.df.summary <- elev.shift.all.df %>% 
  summarize( median = quantile(elev.shift,prob = 0.5),
             upper = quantile(elev.shift,prob = 0.975),
             lower = quantile(elev.shift,prob = 0.025),
             lower25 = quantile(elev.shift,prob = 0.25),
             upper75 = quantile(elev.shift,prob = 0.75))
# calculate density of median shifts to median bar goes to height of density plot
# remove the extra lines
# make median lines same color and interquartile range
elev.dens.max = density(elev.shift.all.df$elev.shift)$y[ floor(density(elev.shift.all.df$elev.shift)$x) == ceiling(quantile( c(slopes.peak.algae*8), prob = 0.5 )) ] / max(density(elev.shift.all.df$elev.shift)$y)
elev.dens <- data.frame( x = quantile( c(slopes.peak.algae*8), prob = 0.5 ),
                         xend = quantile( c(slopes.peak.algae*8), prob = 0.5 ),
                         y = 0,
                         yend = elev.dens.max )

abun.shift.all.df <- data.frame( abun.shift = c(abun.shifts.run.algae) )
minabundiff <- abs(density(log(abun.shift.all.df$abun.shift,base = 2))$x*100 - quantile( log(abun.shifts.run.algae,base=2), prob = 0.5 )*100)
density(log(abun.shift.all.df$abun.shift,base = 2))$y[ which( minabundiff == min(minabundiff) ) ]
abun.dens.max = density(log(abun.shift.all.df$abun.shift,base = 2))$y[ which( minabundiff == min(minabundiff) ) ]  / max(density(log(abun.shift.all.df$abun.shift,base = 2))$y)
abun.dens <- data.frame( x = quantile( log(abun.shifts.run.algae,base=2), prob = 0.5 ), 
                         xend = quantile( log(abun.shifts.run.algae,base=2), prob = 0.5 ),
                         y = 0,
                         yend = abun.dens.max )

#
ybar = 0
error_color = 'black'
dens_height = 1.1
dens_min = -0.1
dens_alpha = 0.5
dens_fill = NA
ydens <- axis_canvas(xy, axis = "y", coord_flip = FALSE)+
  # geom_vline(xintercept=mean(compare_all_plot$elev.shifts.med), col='red' ) +
  geom_vline(xintercept=0) +
  geom_density(data = elev.shift.all.df, aes(x = elev.shift, ..scaled..),
               alpha = dens_alpha, size = 0.5, outline.type = "full", fill = dens_fill, inherit.aes = FALSE) +
  geom_errorbarh(data = elev.shift.all.df.summary, aes(y = ybar, xmin = lower25, xmax = upper75),
                height = 0, size = 2, col = error_color) +
  geom_segment( data = elev.dens, aes(x = x, y = y, xend = xend, yend = yend ), col = error_color, size = 1, lineend = "round") +
  # geom_point(data = elev.shift.all.df.summary, aes(y = ybar, x = median ), size = 1.5, shape = 3, col = "lightgrey") +
  coord_flip(ylim = c(dens_min,dens_height), xlim = c(-125,125) ) 
ydens
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
hist( log(abun.shift.all.df$abun.shift,base=2), freq = F, las = 1)
abun.shift.all.df.summary <- log(abun.shift.all.df,base=2) %>% 
  summarize( median = quantile(abun.shift,prob = 0.5),
             upper = quantile(abun.shift,prob = 0.975),
             lower = quantile(abun.shift,prob = 0.025),
             lower25 = quantile(abun.shift,prob = 0.25),
             upper75 = quantile(abun.shift,prob = 0.75))
xdens <- axis_canvas(xy, axis = "x", coord_flip = TRUE)+
  # geom_vline(xintercept=mean(log(compare_all_plot$abun.shifts.med,base=2)), col='red' ) +
  geom_vline(xintercept=0) +
  geom_density(data = log(abun.shift.all.df,base = 2), aes(x =  abun.shift, ..scaled.. ),
               alpha = dens_alpha, size = 0.5, outline.type = "full", fill = dens_fill) +
  geom_errorbarh(data = abun.shift.all.df.summary, aes(y = ybar, xmin = lower25, xmax = upper75),
                 height = 0, size = 2, col = error_color)+
  # geom_point(data = abun.shift.all.df.summary, aes(y = ybar, x = median ), size = 1.5, shape = 3, col = "seagreen") +
  geom_segment( data = abun.dens, aes(x = x, y = y, xend = xend, yend = yend ), col = error_color, size = 1, lineend = "round") +
  coord_cartesian(ylim = c(dens_min,dens_height), xlim = c(log(1/32,base=2), log(32,base=2)) )
xdens
# an empty plot for the upper corner
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

p1 <- insert_xaxis_grob(xy, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
ggsave(file="R/Figs/abundance~peak.svg",width = 3.5, height = 3.5)

write_csv( compare_all_plot_fun, "R/output/shifts_predicted.csv")
#


# add plot for functional groups through time (fun1 using raw data means and predictions)
dplot2 <- read_csv( "R/output/funtional_groups_annual_mean.csv")
bareraw <- read_csv( "R/output/bare.csv")
dplot2$FunGroup <- factor( dplot2$FunGroup, levels = c("canopy","blade","crust","thin turf","turf","animal") )
dplot2$FunGroup
# numbers of taxa in each group
d.simple %>% 
  group_by( funct_2021) %>% 
  summarize( S = length(unique(taxon)))
dplot2$`Functional Group` <-   factor( dplot2$`Functional Group`, levels = c("canopy (12)","blade (12)","crust (13)","thin turf (28)","turf (41)","animal (10)")) 
qpred_median$`Functional Group` <- factor( qpred_median$FG, levels = c("canopy","blade","crust","thin_turf","turf","animal"),
                                           labels = c("canopy (12)","blade (12)","crust (13)","thin turf (28)","turf (41)","animal (10)") )
qpred_median$Year <- as.numeric(as.character(qpred_median$year))
# add vertical lines for ends of densities
max_pred_fg <- qpred_median %>% group_by(year) %>% summarize(yend = sum(value)) %>% filter(year %in% c(2012,2019)) %>% select(yend)
segs <- data.frame( x = c(2012,2019), xend = c(2012,2019), y = 0, yend = max_pred_fg)
fun1 <- ggplot(dplot2, aes(x = Year, y = mean)) +
  geom_area(data = qpred_median, aes(x = Year, y = value, group = `Functional Group`, fill = `Functional Group`), 
            col="black", alpha = 0.75) +
  # geom_segment( data = segs, aes(x = x,xend = xend,y = y,yend = yend)) +
  geom_bar(aes(fill = `Functional Group`), position="stack", 
           stat="identity", col='black', lwd=0.25, width = 0.5)+
  geom_smooth( data=bareraw, aes(x=Year,y=Abundance, group=1), 
               fill="black",col="yellow" ) +
  theme_classic()+
  scale_fill_manual(values = c("midnightblue", "darkred", "red","pink", "darkgrey", "#996633") %>% rev())+  #"darkgreen",
  theme( legend.position = 'top',
         legend.title = element_blank(),
         legend.key.size = unit(0.25, "cm"),
         legend.box.margin=margin(-8,-10,-3,-10) ) +
  guides( fill =  guide_legend(nrow=2,byrow=T) ) +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         legend.text=element_text(size = 6.5)) +
  ylab("Mean cover (%)") +
  coord_cartesian( ylim = c(0,180) ) 
fun1

# add plot for functional group trends
hmsc_FG


cowplot::plot_grid( fun1, hmsc_FG, p2, ncol = 3, 
                    align='hv', axis="b",
                    rel_widths = c(1.125,1.5,1.25), labels="auto" )
ggsave(file="R/Figs/fun_hmsc_shift.svg",width = 10, height = 10/3)



# see file Functional_groups_throughTime.R for prod_empir_trend_transect
a <- cowplot::plot_grid( fun1, prod_empir_trend_transect,
                    ncol = 2, 
                    align='hv', axis="tblr",
                    rel_widths = c(1, 1), labels="auto" )
a
b <- cowplot::plot_grid( hmsc_FG, p2, ncol = 2, 
                    align='hv', axis="b",
                    rel_widths = c(1,1), labels=c("c","d") )
b
c <- cowplot::plot_grid( a, b, ncol = 1, 
                    align='hv', axis="b", labels = NULL,
                    rel_heights = c(1,1))

ggsave(file="R/Figs/fun_hmsc_shift_revision.svg",width = 5.5, height = 5.5 )


#
