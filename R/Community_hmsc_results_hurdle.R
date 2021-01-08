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
library( vioplot )
library( rstan )
library(scales)
source( "R/mcmc.list2array.R")

## useful references for this code
# citation( "Hmsc" )
# https://github.com/hmsc-r/HMSC

## directories
ModelDir = paste0( here::here(), "/R/models" )
MixingDir = paste0( here::here(), "/R/mixing")



## load the model
list.files( ModelDir )
model = "model_elevxyear_hurdle_chains_4_thin_100_samples_250.Rdata"
mload <- load( paste(ModelDir,model, sep="/") )

# load the data
metacomm <- read_csv("R/output/community.csv" )

##
## MCMC convergence

# We evaluate MCMC convergence in terms of four kinds of parameters that we are especially interested in: the species niches `Beta`, the influence of traits on species niches `Gamma`, residual species associations `Omega`, and the strength of phylogenetic signal `rho`. As there are `r ns` species, the matrix `Omega` is a `r ns` x `r ns` matrix with `r ns*ns` elements. To avoid excessive computational times, we evaluate MCMC convergence for `Omega` only for a subsample of 100 randomly selected species pairs.
mpost_pa   = convertToCodaObject(models[[1]])
mpost_abun = convertToCodaObject(models[[2]])
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
    ma=ge.beta#[,1]
  } else {
    ma = cbind(ma,ge.beta)#[,1])
  }
}
png(file="Figs/MCMC_convergence_hurdle.png")
windows(4.5,6)
par(mfrow=c(2,1),las=1, mar=c(3,4,2,1))
vioplot(ma,col=rainbow_hcl(nm),names=c("presence-absence","abundance-cop"),
        ylim=c(0.9,max(ma)),main="psrf(beta)", ylab="Gelman statistic" )
abline(a=1.05,b=0)
vioplot(ma,col=rainbow_hcl(nm),names=c("presence-absence","abundance-cop"),
        ylim=c(0.9,1.1),main="psrf(beta)", ylab="Gelman statistic" )
abline(a=1.05,b=0)
dev.off()



## Assess model fit
preds = lapply( models, computePredictedValues )
MF <- mapply( evaluateModelFit, models, preds )
# windows(5,5)
for(i in 1:length(models)){
  print(lapply( MF[[i]], summary))
}




## parameter estimates
postBeta = lapply( models, getPostEstimate, parName = "Beta")
# windows(5,8)
# plotBeta(m, post = postBeta, param = "Sign", supportLevel = 0.95, mar=c(7,11,0,6))
postBeta[[1]]$mean[, c("Alaria.marginata","Hedophyllum.sessile","Polysiphonia")]
postBeta[[2]]$mean[, c("Alaria.marginata","Hedophyllum.sessile","Polysiphonia")]

cor( as.vector(postBeta[[1]]$mean), as.vector(postBeta[[2]]$mean) )
plot( as.vector(postBeta[[1]]$mean), as.vector(postBeta[[2]]$mean) )

pos.negs <- NULL
for(i in 1:length(models)){
  pos.neg <- data.frame(pos = c(postBeta[[i]]$support), neg = c(postBeta[[i]]$supportNeg))
  pos.neg[pos.neg< 0.95] <- 0
  pos.neg$neg <- -pos.neg$neg
  pos.neg$value <- pos.neg$pos + pos.neg$neg
  names(pos.neg) <- paste0( c('pos','neg','value'),i)
  if(is.null(pos.negs)){
    pos.negs=pos.neg#[,1]
  } else {
    pos.negs = cbind(pos.negs,pos.neg)#[,1])
  }
  pos.negs$parameter <- factor(c("intercept", "year1", "year2", "elev1",
                                "elev2", "elev1:year1", "elev1:year2" ),
                              levels = c("intercept", "year1", "year2", "elev1",
                                         "elev2", "elev1:year1", "elev1:year2" ), 
                              ordered = TRUE)
  pos.negs$species <- factor(rep(colnames(postBeta[[i]]$mean), each = 7), 
                            levels = colnames(models[[i]]$Y)[order(colSums(models[[i]]$Y),decreasing = TRUE)],
                            ordered = TRUE)
}


windows(12,8)
support1 <- ggplot(pos.negs, aes(y = parameter, x = species, fill = value1))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")
support2 <- ggplot(pos.negs, aes(y = parameter, x = species, fill = value2))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")
plot_grid( support1, support2, ncol=1 )


# # pull taxa with positive estimated response to temperature anomaly
# taxa.pos.anom <- as.character(pos.neg[ pos.neg$parameter=='temp.anom' & pos.neg$value>0, 'species'])

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
VP = lapply( models, computeVariancePartitioning ) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP)
VP.dfs <- NULL
for( i in 1:length(models) ){
  VP.df <- as.data.frame(VP[[i]]$vals) %>% 
    mutate(effect = factor(c("year1","year2","elev1","elev2",
                             "elev:year","elev2:year",
                             "site","transect"), 
                           levels = rev(c("year1","year2","elev1","elev2",
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

windows(8,5)
ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  facet_wrap(~model,ncol=1) +
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon", viridis(7)), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.07, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")





## associations
OmegaCor = lapply( models, computeAssociations )
supportLevel = 0.95
# choose the random variable to plot
rlevel = 1
pick <- 1
toPlot = ((OmegaCor[[pick]][[rlevel]]$support>supportLevel)
          + (OmegaCor[[pick]][[rlevel]]$support<(1-supportLevel))>0)*OmegaCor[[pick]][[rlevel]]$mean
# reorder species matrix
plotorder <- order( postBeta[[pick]]$mean[2,], decreasing = TRUE )
toPlot <- toPlot[ plotorder, plotorder]
# reorder automatically
library(lessR)
mynewcor <- corReorder( toPlot, order="hclust", nclusters=4 )
# windows(12,12)
corrplot( mynewcor, method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
           title = paste("random effect level:", models[[1]]$rLNames[rlevel]), mar=c(0,0,0.5,0), tl.cex=0.6 )





#####
# Model  predictions
# response curves across all combinations
newXData <- data.frame(shore.height = seq(60,380,by = 1),
                       anom.pine.sum.1 = 1.1706,
                       anom.pine.win = 0.8530 )
years <- as.data.frame( models[[1]]$XData %>% select(year,year1,year2) %>% mutate(year1=round(year1,7),year2=round(year2,7)) %>% distinct() )
elevs <- as.data.frame( models[[1]]$XData %>% select(elev,elev1,elev2) %>% mutate(elev1=round(elev1,5),elev2=round(elev2,5)) %>% distinct() %>% arrange(elev1) )
plot(elevs)
plot(years)

newDF <- expand.grid( shore.height = seq(60,380,by = 1),
             year=c(2012:2019),
             site="new unit",
             transect="new unit" )
newDF <- data.frame( merge(years, elevs),
             site="new unit",
             transect="new unit" )
# trim dataset to every nth observation, only particular years
newDF <- data.frame( merge(years[years$year %in% c(2012,2019),], elevs[seq(1,170, by = 3),]),
                     site="new unit",
                     transect="new unit" )

newDFsel <- newDF %>% select(year1,year2,elev1,elev2,site,transect)
newXData   <- newDFsel[,1:(ncol(newDFsel)-2)]
newDesign  <- newDFsel[,(ncol(newDFsel)-1):ncol(newDFsel)]

## predictions of individual models
predY_pa <- predict(models[[1]], XData = newXData,
                 studyDesign= newDesign,
                 ranLevels = list(site=rL_site,transect=rL), expected = TRUE) #, transect=rL, year=rL_year
predY_cop <- predict(models[[2]], XData = newXData,
                 studyDesign= newDesign,
                 ranLevels = list(site=rL_site,transect=rL), expected = TRUE) #, transect=rL, year=rL_year
## predictions of both models multiplied together
predY_abun <- Map('*', predY_pa, lapply(predY_cop,exp) )
# predY_pa[[1]][1,1] * predY_cop[[1]][1,1]
# predY_abun[[1]][1,1]


#####
## get peak elevation and "abundance" for each species in YEAR in each run
# elevation peaks for each species in each YEAR in each run
peaks  <- lapply( predY_abun, 
                  function(i) lapply( split( i, newDF$year ), 
                          function(l) apply(matrix(l,byrow = F,ncol = 47), 2, 
                                            function(z) unique(newDF$elev)[which(z==max(z))]) ) )
peaks_bind <- lapply( peaks, function(z) do.call(rbind,z) )
peaks_array <- abind::abind(peaks_bind, along=3)
# abundances of each species in each YEAR in each run
abunds  <- lapply( predY_abun, 
                  function(i) lapply( split( i, newDF$year ), 
                                      function(l) apply(matrix(l,byrow = F,ncol = 47), 2,sum) ) )
abunds_bind <- lapply( abunds, function(z) do.call(rbind,z) )
abunds_array <- abind::abind(abunds_bind, along=3)

## calculate differences between 2012 and 2019 to get shift in end member states
# elevation peak
elev.shifts.run     <- apply( peaks_array, c(2,3), function(z) z[8]-z[1] )
elev.init           <- apply( peaks_array, c(2,3), function(z) z[1] )
elev.init.med       <- apply(elev.init, 1, quantile, prob = 0.5, na.rm = TRUE)
elev.shifts.med     <- apply(elev.shifts.run, 1, quantile, prob = 0.5, na.rm = TRUE)
elev.shifts.high    <- apply(elev.shifts.run, 1, quantile, prob = 0.025, na.rm = TRUE)
elev.shifts.low     <- apply(elev.shifts.run, 1, quantile, prob = 0.975, na.rm = TRUE)
elev.shifts.summary <- data.frame(elev.init.med, elev.shifts.med, elev.shifts.low, elev.shifts.high)
summary( lm(elev.shifts.med~1) )

# abundance
abun.shifts.run     <- apply( abunds_array, c(2,3), function(z) z[8]/z[1] )
abun.init           <- apply( abunds_array, c(2,3), function(z) z[1] )
abun.init.med       <- apply(abun.init, 1, quantile, prob = 0.5, na.rm = TRUE)
abun.shifts.med     <- apply(abun.shifts.run, 1, quantile, prob = 0.5, na.rm = TRUE)
abun.shifts.high    <- apply(abun.shifts.run, 1, quantile, prob = 0.025, na.rm = TRUE)
abun.shifts.low     <- apply(abun.shifts.run, 1, quantile, prob = 0.975, na.rm = TRUE)
abun.shifts.summary <- data.frame(abun.init.med,abun.shifts.med, abun.shifts.low, abun.shifts.high)
summary( lm(log(abun.shifts.med,base=2)~1) )

## combine these results
shift.summary <- data.frame( elev.shifts.summary, abun.shifts.summary )
shift.summary$taxon <- colnames(Y)
shift.summary %>% arrange(-abun.init.med)
shift.summary$rank <- 1:47
ggplot( shift.summary, aes(x = log(abun.shifts.med,base=2), y = elev.shifts.med) ) +
  # geom_hline( yintercept=0, lty=2 ) + geom_vline( xintercept = 0, lty=2 ) +
  geom_linerange( aes( ymin = elev.shifts.low, ymax= elev.shifts.high), alpha=0.25) +
  geom_linerange( aes(xmin = log(abun.shifts.low,base=2), 
                      xmax = log(abun.shifts.high,base=2)), alpha=0.25 ) +
  geom_point() + 
  scale_x_continuous( breaks=c(4,2,0,-2,-4), 
                      labels=c('16x','4x','0','1/4x','1/16x')) +
  ylab( "Elevation shift (cm)" ) + xlab( "Abundance shift" ) +
  theme_bw() + theme( panel.grid.minor = element_blank() )
# ggsave( "R/Figs/shifts_error.pdf", width=4, height=4 )
# correlation of median responses
with( shift.summary, cor.test( elev.shifts.med, abun.shifts.med, method = 'spearman' ) )
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
elev.shift.plot <- ggplot( elev.shifts.summary, aes(x=1:47,y=elev.shifts.med)) + 
  geom_hline( yintercept=0 ) +
  geom_errorbar( aes(ymin=elev.shifts.low,ymax=elev.shifts.high) ) +
  geom_point() +
  ylab( "Elevation peak shift (cm)" ) + xlab("") +
  scale_x_continuous(breaks = c(1,10,20,30,40,47)) +
  theme( panel.grid.minor.x = element_blank() )
abun.shift.plot <- ggplot( abun.shifts.summary, aes(x=1:47,y=log(abun.shifts.med,base=2)) ) + 
  geom_hline( yintercept=0 ) +
  geom_errorbar( aes(ymin=log(abun.shifts.low,base=2),ymax=log(abun.shifts.high,base=2)) ) +
  geom_point() +
  ylab( "Percent cover shift" ) + xlab("Rank occurrence") +
  scale_y_continuous( breaks=c(6,4,2,0,-2,-4,-6), 
                      labels=c('64x','16x','4x','0','1/4x','1/16x','1/64x')) +
  scale_x_continuous( breaks = c(1,10,20,30,40,47)) +
  theme( panel.grid.minor.x = element_blank() )
plot_grid( elev.shift.plot, abun.shift.plot, ncol=1 )
ggsave( "R/Figs/shifts_2panel.pdf", width=6, height=5 )




# ###
### other ways to summarize the model results
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



# plot responses
predictions_pa$year <- as.character(factor(predictions_pa$year1, levels = unique(predictions_pa$year1), labels = unique(newDF$year)))
predictions_cop$year <- as.character(factor(predictions_pa$year1, levels = unique(predictions_pa$year1), labels = unique(newDF$year)))
predictions_abun$year <- as.character(factor(predictions_pa$year1, levels = unique(predictions_pa$year1), labels = unique(newDF$year)))
predictions_pa$elev <- as.numeric(as.character(factor(predictions_pa$elev1, levels = unique(predictions_pa$elev1), labels = unique(newDF$elev))))
predictions_cop$elev <- as.numeric(as.character(factor(predictions_pa$elev1, levels = unique(predictions_pa$elev1), labels = unique(newDF$elev))))
predictions_abun$elev <- as.numeric(as.character(factor(predictions_pa$elev1, levels = unique(predictions_pa$elev1), labels = unique(newDF$elev))))


# color scheme
as.survey <- read_csv(  "R/output/sst_anoms_survey.csv" )
as.survey$year <-  as.character(as.survey$year)
library( RColorBrewer )
anom.range <- c(-2,2)
n=9
cols <- brewer.pal(n,"RdBu")
pal <- colorRampPalette(rev(cols))

predictions_pa <- left_join(predictions_pa, as.survey)
predictions_cop <- left_join(predictions_cop, as.survey)
predictions_abun <- left_join(predictions_abun, as.survey)
# colors

cols.two <- c( rgb( 211,230,240, maxColorValue=255), rgb( 232,139,110, maxColorValue=255))

hist( comm$Fucus.distichus )
taxon <- "Hedophyllum.sessile"

a <- ggplot(predictions_pa, aes_string(x = 'elev', y = taxon, col='year'))+
  geom_smooth(aes(group=year1),se=F,lwd=1.5) +
  scale_color_manual(values=cols.two) +
  ylab("") +
  xlab("") +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b <- ggplot(predictions_cop, aes_string(x = 'elev', y = taxon, col='year'))+
  geom_smooth(aes(group=year1),se=F,lwd=1.5) +
  scale_color_manual(values=cols.two) +
  ylab("") +
  xlab("") +
  theme_classic() + 
  # theme(legend.position = "none") +
  theme(legend.position = c(1,0.25), legend.justification = c(1,0) ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
c <- ggplot(predictions_abun, aes_string(x = 'elev', y = taxon, col='year'))+
  geom_smooth(aes(group=year1),se=F,lwd=1.5) +
  # scale_color_gradient(low = "grey75", high="gray20") +
  # scale_color_gradientn(colours=pal2(100),limits=c(-2,2)) +
  scale_color_manual(values=cols.two) +
  # ylab("percent cover") +
  # xlab("elevation (cm)") +
  ylab("") +
  xlab("") +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
plot_grid(a,b,c,ncol=1, align='hv')
ggsave( paste0("R/Figs/hmsc_elev_metrics_",taxon,".svg"), height=4.5, width=1.5 )

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
predictions_abund <- left_join(predictions_abund,newDF)
save( predictions_abund, file = paste("R/output/hmsc_pred",model, sep="_") )
# ###


load( file = paste("R/output/hmsc_pred",model, sep="_") )
# load( file = paste("R/output/hmsc_pred_model_5_chains_4_thin_100_samples_1000.Rdata") )






##### Counter-factual predicitons
## predictions over on factor only
# get predictions over time
Gradient <- constructGradient(models[[1]], focalVariable = "year1", non.focalVariables = list("elev1" = list(1)),
                              ngrid = 8) #length(unique(models[[1]]$XData$year1)))

Gradient$XDataNew$year1
XData_choose <- models[[1]]$XData %>% 
  select(year,year1,year2) %>% mutate_all(round,7) %>%  distinct()
XData_choose_elev <- models[[1]]$XData %>% 
  select(elev,elev1,elev2) %>% mutate_all(round,7) %>%  distinct()
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

species_occur <- left_join(apply(scaled_occur, c(1,2), median) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "occurrence") %>% gather(key = Species, value = median, -year1, -metric),
                           apply(scaled_occur, c(1,2), quantile, prob = 0.25) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "occurrence") %>% gather(key = Species, value = quant_0.25, -year1, -metric)) %>%
  left_join(apply(scaled_occur, c(1,2), quantile, prob = 0.75) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "occurrence") %>% gather(key = Species, value = quant_0.75, -year1, -metric))

species_occur <- species_occur %>%
  gather(key = quant, value = value, median:quant_0.75) %>%
  mutate(value = scale(value)) %>%
  spread(key = quant, value = value)

species_con_bmass <- left_join(apply(log_con_bmass, c(1,2), median) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "conditional cover") %>% gather(key = Species, value = median, -year1, -metric),
                               apply(log_con_bmass, c(1,2), quantile, prob = 0.25) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "conditional cover") %>% gather(key = Species, value = quant_0.25, -year1, -metric)) %>%
  left_join(apply(log_con_bmass, c(1,2), quantile, prob = 0.75) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "conditional cover") %>% gather(key = Species, value = quant_0.75, -year1, -metric))

species_con_bmass <- species_con_bmass %>%
  gather(key = quant, value = value, median:quant_0.75) %>%
  mutate(value = scale(value)) %>%
  spread(key = quant, value = value)

species_bmass <- left_join(apply(log_biomass, c(1,2), median) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "cover") %>% gather(key = Species, value = median, -year1, -metric),
                           apply(log_biomass, c(1,2), quantile, prob = 0.25) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "cover") %>% gather(key = Species, value = quant_0.25, -year1, -metric)) %>%
  left_join(apply(log_biomass, c(1,2), quantile, prob = 0.75) %>%  as.data.frame() %>% mutate(year1 = Gradient$XDataNew$year1, metric = "cover") %>% gather(key = Species, value = quant_0.75, -year1, -metric))

species_bmass <- species_bmass %>%
  gather(key = quant, value = value, median:quant_0.75) %>%
  mutate(value = scale(value)) %>%
  spread(key = quant, value = value)

species_temporal.df <- rbind(species_occur, species_con_bmass, species_bmass)
species_temporal.df <- left_join( mutate(species_temporal.df,year1=round(year1,7)), years )

library(broom)
sp_scaled_trends <- species_temporal.df %>%
  ungroup() %>%
  dplyr::select(-quant_0.25, -quant_0.75) %>%
  nest_by(Species, metric) %>%
  mutate(fitYear = list(lm(median ~ year1, data = data))) %>%
  summarise(tidy(fitYear)) %>%
  filter(term != "(Intercept)") %>%
  mutate(estimate_sig = estimate * as.numeric(p.value<0.05))

order_occur <- sp_scaled_trends %>%
  filter(metric == "cover") %>%
  arrange(estimate_sig) 

sp_scaled_trends$Species <- factor(sp_scaled_trends$Species, levels = order_occur$Species, ordered = TRUE)
sp_scaled_trends$metric <- factor(sp_scaled_trends$metric, levels = c('occurrence','conditional cover','cover'), ordered = FALSE)


ggplot(sp_scaled_trends, aes(x = metric, y = Species, fill = estimate_sig))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue3", high = "red", name = "scaled\nchange/year")+
  ylab("")+
  xlab("") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave("R/Figs/hmsc_scale_change_heatplot.svg", height = 6.2, width = 5)

# boxplot of estimates
fill_cols <- c('mintcream','mediumseagreen','mediumspringgreen')
a <- ggplot(sp_scaled_trends, aes(x = metric, y = estimate, fill=metric)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot(width=0.5) +
  xlab("Metric") + ylab("Estimate") + 
  theme_bw() + 
  theme( legend.position = 'none',panel.grid = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual( values = fill_cols)
a
ggsave("R/Figs/hmsc_scale_change_boxplot.png", height = 4, width = 4)

# add functional groups
sp_scaled_trends$taxon <- sp_scaled_trends$Species
sp_scaled_trends_fun <- left_join( sp_scaled_trends, taxon.key  )
sp_scaled_trends_fun$funct_Sep2020[ is.na(sp_scaled_trends_fun$funct_Sep2020)] <- "animal"
sp_scaled_trends_fun$funct_Sep2020 <- factor(sp_scaled_trends_fun$funct_Sep2020,
                                             levels = c('canopy','turf','thin_turf','crust','blade','animal'))
b <- ggplot(sp_scaled_trends_fun, aes(x = funct_Sep2020, y = estimate, fill=metric)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot(width=0.5) +
  xlab("Functional Group") + ylab("") + 
  theme_bw() +
  theme( legend.position = 'none',panel.grid = element_blank()) +
  theme( panel.grid = element_blank())+
  scale_fill_manual( values = fill_cols)
b
plot_grid(a,NULL,b,nrow=1,rel_widths = c(0.5,-0.01,1),
          align = 'hv', labels=c("A","","B") )
ggsave("R/Figs/hmsc_scale_change_boxplot_combo.png", height = 3, width = 6)

# combine datasets for boxplots
sp_scaled_trends$funct_Sep2020 <- "pooled"
sp_scaled_trends_comb <- bind_rows( sp_scaled_trends, sp_scaled_trends_fun )
# remove hyphen from functional group names
sp_scaled_trends_comb$funct_Sep2020 <- gsub("_"," ",as.character(sp_scaled_trends_comb$funct_Sep2020))
sp_scaled_trends_comb$funct_Sep2020 <- factor(sp_scaled_trends_comb$funct_Sep2020,
                                             levels = c('pooled','canopy','turf','thin turf','crust','blade','animal'))
# # colors
# fill_cols <- c('mintcream','mediumseagreen','mediumspringgreen')
sp_scaled_trends_comb$colors <- "mintcream" 
sp_scaled_trends_comb$colors[sp_scaled_trends_comb$metric == 'conditional cover'] <- "mediumseagreen" 
sp_scaled_trends_comb$colors[sp_scaled_trends_comb$metric == 'cover'] <- "mediumspringgreen" 
# fill_cols <- c('mintcream','mediumseagreen','mediumspringgreen')
# line types, colors
sp_scaled_trends_comb$line.type <- 1
sp_scaled_trends_comb$line.type[ sp_scaled_trends_comb$funct_Sep2020 == "pooled"] <- 3
sp_scaled_trends_comb$a <- 0.75
sp_scaled_trends_comb$a[ sp_scaled_trends_comb$funct_Sep2020 == "pooled"] <- 1



# get the number of datapoints for each group
reps <- sp_scaled_trends_comb %>% 
  group_by(funct_Sep2020) %>% 
  summarize(n=length(estimate)/3) %>% 
  mutate(estimate = -22.5, metric = 'conditional cover')

ggplot(sp_scaled_trends_comb, aes(x = funct_Sep2020, y = estimate, fill=metric)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_text(data = reps, aes(label = paste0('(',n,')')), size = 3) +
  xlab("Functional Group") + ylab("Estimate") + 
  ylim(c(-23,13)) +
  theme_bw() +
  theme( legend.position="top",
         panel.grid = element_blank(),
         legend.title = element_text(size=10),
         legend.text = element_text(size=8),
         legend.key.size = unit(0.5, "cm"),
         legend.key = element_rect(colour = NA, fill = NA)) +
  scale_fill_manual( values = fill_cols)
ggsave("R/Figs/hmsc_scale_change_boxplot_comb_single.png", height = 3.5, width = 4)



## plot trends for all species
species_temporal.df$Species <- factor(species_temporal.df$Species, levels = order_occur$Species, ordered = TRUE)
species_temporal.df$Species <- factor(species_temporal.df$Species, levels = order_occur$Species, labels = gsub("[.]","\n",order_occur$Species), ordered = TRUE)
species_temporal.df$metric  <- factor(species_temporal.df$metric, levels = c("occurrence","conditional cover","cover"))

sp.trends <- ggplot(species_temporal.df, aes(x = year1, y = median, color = metric,  fill = metric, group = metric))+ #fill = metric,
  geom_ribbon(aes(ymin = quant_0.25, ymax= quant_0.75), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  facet_wrap(~Species, scales = "free", ncol = 6)+
  # scale_color_brewer(type = "qual", palette = "Dark2", name = "")+
  # scale_fill_brewer(type = "qual", palette = "Dark2", name = "", guide = FALSE)+
  scale_color_manual( values = c('gray25','seagreen','springgreen') ) +
  scale_fill_manual( values = c('gray50','seagreen','springgreen3') ) +
  xlab("")+
  ylab("")+
  theme_classic() + theme(legend.position = "top")
ggsave(sp.trends, "R/Figs/hmsc_sp_trends.svg", height = 11*1.5, width = 8.5*1.5)
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

point_colors <- c("slateblue","firebrick", "goldenrod" )
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
species <-list( c("Fucus.distichus","Elachista.fucicola"),c("Palmaria.hecatensis","Palmaria.mollis") ,c("Polysiphonia","Lithothamnion.phymatodeum"),"Pyropia",c("Alaria.marginata","Hedophyllum.sessile"),c( "Mazzaella.parvula", "Mazzaella.oregona", "Mazzaella.splendens" ),"Cladophora.columbiana","Corallina",c("Colpomenia.bullosa","Colpomenia.peregrina"),c("Mastocarpus","Petrocelis"),c("Phyllospadix.sp."),c("Gloiopeltis.furcata"),c("Barnacles","Mytilus.sp."),c("Costaria.costata","Osmundea.spectabilis","Nemalion.helminthoides"),c("Scytosiphon.lomentaria","Lomentaria.hakodatensis","Salishia.firma"),c("Erythrotrichia.carnea","Ectocarpus.commensalis","Elachista.fucicola"))

species_temporal.df$taxon <- gsub("\n", ".", species_temporal.df$Species)

alpha_choose = 0.1
alpha_choose2 = 1
size_choose = 2.5
size_choose2 = 5
for(i in 1:length(species)){
  abun.time.plot <- species_temporal.df %>% 
  filter(taxon %in% species[[i]] ) %>% ungroup() %>% 
  filter( metric == "conditional cover" ) %>% 
  mutate( N = median )
  
  ggplot( abun.time.plot, aes( x=year, y=exp(N), group=taxon, col=taxon ) ) + 
  # geom_ribbon( aes(x=year, ymin=N_high/90, ymax=N_low/90), fill="grey70", alpha=0.5 ) +
  # facet_wrap(~taxon, scales="free_y",ncol=1) +
  geom_path(size=size_choose, alpha=alpha_choose2) +
  # geom_point( data=filter(ogd.mean.pa, taxon %in% species ), aes(y=N+0.01111111), size=5, alpha = 0.1 ) +
  geom_point( data=filter(ogd.mean, taxon %in% species[[i]] ), aes(y=N+0.01111111), size=size_choose2, 
              alpha = alpha_choose, position = position_dodge(width = 0.25) ) +
  geom_path( data=filter(ogd.mean.pa.mean, taxon %in% species[[i]] ), aes(y=N+0.01111111), lwd=0.75, alpha = 1 ) +
  # stat_summary( data=filter(ogd.mean.pa, taxon %in% species ), aes(y=N), fun = "median", geom="line", size = 0.5 ) +
  scale_color_manual(values=point_colors[1:length(species[[i]])]) +
  ylab("Percent Cover") +
  theme_classic() +
  theme(legend.position="top") +
  # ylim( c(0,18) ) +
  # scale_y_sqrt(labels = comma, breaks = c(0, 1.01111111, 10, 25, 50, 75, 100)) +
  scale_y_log10(labels = comma) +
  # scale_y_continuous(labels = comma) +
  guides(color=guide_legend(nrow=2,ncol=2,byrow=F))  
  
  ggsave( paste0("R/Figs/temporal_trends_select_taxa/",paste0(species[[i]],collapse="_"),"_model+data.svg"), 
        width=3,height=3.5)
}
#















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
comm_final$taxon <- factor(comm_final$taxon, 
                           levels = colnames(models[[1]]$Y)[order(colSums(models[[2]]$Y),decreasing = T)], ordered = TRUE)

inverts <- c("Barnacles","Mytilus.sp.","Anemone","Bryozoan","Tunicata/Porifera","Pollicipes.polymerus","Tube.worms","Hydroid" )
comm_final$alga <- "alga"
comm_final$alga[ comm_final$taxon %in% inverts ] <- "invert"
comm_final$elev <- comm_final$shore.height
predictions_abund$alga <- "alga"
predictions_abund$alga[ predictions_abund$taxon %in% inverts ] <- "invert"
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
FunGroups$funct_Sep2020 <- gsub("large_brown","canopy",FunGroups$funct_Sep2020)
FunGroups$funct_Sep2020 <- gsub("seagrass","canopy",FunGroups$funct_Sep2020)
FunGroups$taxon[ FunGroups$taxon == "Bossiella.articulate" ] <- "Bossiella_articulate"

# Check which species names don't match #These should all be animals and other non-algal fields
colnames(comm)[colnames(comm) %in% FunGroups$taxon == "FALSE"] 
# fix naming discrepancies
# names(comm)[colnames(comm)=="Bossiella_articulate"] <- "Bossiella.articulate"

# Creat2 matrix of functional groups summed
taxon<-data.frame(taxon = colnames(comm))
taxon.key<-left_join(taxon, FunGroups, by = "taxon") # animals should be missing

predictions_abund_trait <- left_join( predictions_abund, taxon.key )
sort(unique(predictions_abund_trait$taxon))
# add trait for animals
predictions_abund_trait$funct_Sep2020[ is.na(predictions_abund_trait$funct_Sep2020)] <- "animal"
predictions_abund_trait$funct <- predictions_abund_trait$funct_Sep2020 




# only pull the 10 most common taxa
taxa <- colnames(models[[1]]$Y)
top6 <- taxa[c(1:6,9,14,15)]
noxshift <- taxa[c(3,9,11,15,20,25,26,47)]
upX <- taxa[c(45,42,44,16,37,35)]
downX <- taxa[c(43,36,5,24,28,33,30)]
upY <- taxa[c(15,42,39,47,44,46,20,22,9)]
customXY <- taxa[c(1,4,5,6,7, 
                   14,15,13,27,
                   30,34,37)]
# 10 rarest taxa
bot10 <- taxa[(length(taxa)-8):length(taxa)]

taxa2plot <- customXY

# windows(6,4)
ggplot( filter(predictions_abund,taxon %in% taxa2plot & year %in% c(2012,2019)), 
        aes(x = elev, y = N,
            fill=factor(year), col=factor(year) ))+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, col="gray75")+
  facet_wrap(~taxon, scales = "free_y", ncol=4)+
  theme_classic() +#+
  # scale_fill_viridis_d() +
  # scale_color_manual(values=c("black","black"))+
  scale_fill_manual(values=c("gold1", "darkslategray4"))+
  geom_point( data = filter( comm_final, taxon %in% taxa2plot, year %in% c(2012,2019)), pch=21, alpha = 0.5 ) +
  geom_line(size = 0.5 ) +
  scale_color_manual(values = c("darkslategray","black")) +
  scale_y_sqrt() + ylab("Percent cover") + xlab("Shore height (cm)") +
  theme( strip.text.x = element_text(size = 7) ) +
  coord_cartesian(ylim=c(0,100))
ggsave("R/Figs/hmsc_response_curves_hurdle.svg", width = 7, height = 5)

# just show Fucus
fuc <- "Fucus.distichus"
fuc <- "Polysiphonia"
fuc <- "Mazzaella.parvula"
fuc <- "Microcladia.borealis"
fuc <- "Palmaria.hecatensis"
fuc <- "Farlowia.mollis"
fuc <- "Cladophora.columbiana"
fuc <- "coralline.crust"
fuc <- "Alaria.marginata"
fuc <- "Egregia.menziesii"
fuc <- "Hedophyllum.sessile"
fuc <- "Lithothamnion.phymatodeum"
fuc <- "Neopolyporolithon.reclinatum"
fuc <- "Mytilus.sp."
fuc <- "Barnacles"
fuc <- "Elachista.fucicola"
fuc <- "Phyllospadix.sp."
fuc <- "Acrosiphonia"
fuc <- "Pyropia"
fuc <- "Petrocelis"
fuc <- "Ralfsioid"
# windows(5,4)
ggplot( filter(predictions_abund,taxon %in% fuc ), 
        aes(x = elev, y = N ))+
  facet_wrap(~year)+
  geom_ribbon(aes(ymin = N_low, ymax = N_high), alpha = 0.5, fill="gray75")+
  geom_line(size = 0.5)+
  theme_classic() +
  ggtitle(fuc) +
  geom_point( data = filter( comm_final, taxon %in% fuc), pch=21 ) +
  scale_y_sqrt(breaks=c(1,10,50,100,200)) +
  ylab("Percent cover") + xlab("Shore height (cm)") +
  coord_cartesian(ylim = c(-0, 100)) 
#
compare_all %>% arrange(shift.y)
compare_all %>% arrange(shift.x)
#


# Find the predicted peak for each instance
peaks <- predictions_abund_trait %>% 
  group_by( year, taxon, funct ) %>% 
  summarize( peak = mean(elev[which(N==max(N))]), peaksd = mean(elev[which(N==max(N))],na.rm=T) ) 
peaks$peak
ggplot( filter(peaks, taxon %in% upX), aes(x=year,y=peak) ) + facet_wrap(~taxon) +
  geom_point()
ggplot( peaks, aes(x=year,y=peak) ) + #facet_wrap(~taxon) +
  # geom_path(aes(group=taxon), alpha=0.5) + 
  geom_smooth(method="lm") + geom_smooth(aes(group=taxon),col='black',se=F,lwd=0.5,alpha=0.5)
summary(lm(peak~1+year,data=peaks))

# get difference between peaks for 2012 and 2019
peak_shift <- peaks %>% 
  group_by( taxon, funct ) %>% 
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
trait.p <- ggplot( peak_shift, aes(x=reorder(funct, shift, FUN = median), y=shift/100) ) + 
  geom_hline (yintercept=0, lty=2 ) +
  geom_boxplot()  + geom_point() +
  xlab("Functional group") + ylab("Peak shift (meters)") +
  # scale_fill_manual(values=c("whitesmoke","dodgerblue")) +
  theme_classic() +
  theme( axis.text.x = element_text(angle=45,hjust=1,vjust=1) ) +
  theme(legend.position = "none")

peak_shift %>% arrange(-shift)
peak_shift %>% arrange(shift)
peak_shift %>% filter(shift==0)
peak_shift %>% filter(shift >-10 & shift < 10)
shift_increase <- peak_shift %>% arrange(-shift)
choose <- shift_increase$taxon[1:6]
ggplot( filter(peaks, taxon %in% choose), aes(x=year,y=peak) ) + facet_wrap(~taxon) +
  geom_point()


summary(lm( shift~1, peak_shift))
peak_compare <- peak_compare %>% 
  ungroup() %>% 
  mutate(peak2=peak-min(peak))
summary(lm( shift~peak2, peak_compare))



# need to add lines showing the realm of possible shifts
xs <- range( XData$elev )
y1 <- c(0,diff(xs))
y2 <- c(-diff(xs),0)
df.bound <- data.frame( x1=xs[1],x2=xs[2],y1,y2 )
df.poly <- data.frame( x=rep(xs,each=2), y=c(0,diff(xs),0,-diff(xs)) )



# (a <- ggplot( peak_compare, aes(x=peak,y=shift)) +
#     geom_polygon( data=df.poly, aes(x=x,y=y), fill='whitesmoke', col='slategray', lty=2) +
#     geom_hline( yintercept = 0, lty=2 ) +
#      geom_smooth(method='lm', se=T, col='black') +
#     geom_point( size=3, pch=1, col='slateblue' ) +
#   ylab("peak elevationshift (cm)") + xlab("initial peak elevation (cm)") +
#   theme_classic() )
summary(lm( elev.shifts.med ~ elev.init.med, data=shift.summary ))
cor.test( x=shift.summary$elev.shifts.med, y = shift.summary$elev.init.med )
(a <- ggplot( shift.summary, aes(x=elev.init.med,y=elev.shifts.med)) +
    geom_polygon( data=df.poly, aes(x=x,y=y), fill='whitesmoke', col='slategray', lty=2) +
    geom_hline( yintercept = 0, lty=2 ) +
    geom_smooth(method='lm', se=T, col='black') +
    geom_point( size=3, pch=1, col='slateblue' ) +
    ylab("Peak elevationshift (cm)") + xlab("Initial peak elevation (cm)") +
    theme_classic() )


# shift in abundance ~ initial peak elevation
# (b <- ggplot( compare_all, aes(x=peak,y=log(shift.y,base=2))) + 
#   geom_hline( yintercept=0, lty=2 ) +
#   geom_point(size=3, pch=1, col='slateblue') + 
#   ylab("abundance shift") + xlab("initial peak elevation (cm)") +
#   scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
#                       labels=c('10x','2x','0','1/2x','1/10x')) +
#   theme_classic() )
(b <- ggplot( shift.summary, aes(x=elev.init.med,y=log(abun.shifts.med,base=2))) + 
    geom_hline( yintercept=0, lty=2 ) +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("Cover shift") + xlab("Initial peak elevation (cm)") +
    scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
                        labels=c('10x','2x','0','1/2x','1/10x')) +
    theme_classic() )

# peak shift ~ initial abundance
# (c <- ggplot( compare_all, aes(x=integral/30,y=shift.x)) + 
#     geom_hline( yintercept=0, lty=2 ) +
#     geom_point(size=3, pch=1, col='slateblue') + 
#     ylab("peak elevationshift (cm)") + xlab("appox. initial cover (%)") +
#     # scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
#     #                     labels=c('10x','2x','0','1/2x','1/10x')) +
#     scale_x_continuous(trans = "log2") +
#     theme_classic() )
(c <- ggplot( shift.summary, aes(x=abun.init.med/30,y=elev.shifts.med)) + 
    geom_hline( yintercept=0, lty=2 ) +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("Peak elevation shift (cm)") + xlab("Appox. initial cover (%)") +
    # scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
    #                     labels=c('10x','2x','0','1/2x','1/10x')) +
    scale_x_continuous(trans = "log2") +
    theme_classic() )

# abund shift ~ initial abundance
# summary(lm( log(shift.y,base=2)~log(integral,base=2), compare_all ))
# (d <- ggplot( compare_all, aes(x=integral/30,y=log(shift.y,base=2))) + 
#     geom_hline( yintercept=0, lty=2 ) +
#     geom_smooth(method = 'lm', se = T, col='black' ) +
#     geom_point(size=3, pch=1, col='slateblue') + 
#     ylab("Abundance shift") + xlab("appox. initial cover (%)") +
#     scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
#                         labels=c('10x','2x','0','1/2x','1/10x')) +
#     scale_x_continuous(trans = "log2") +
#     theme_classic() )
summary(lm( log(abun.shifts.med,base=2)~log(abun.init.med,base=2), shift.summary ))
cor.test( x=log(shift.summary$abun.shifts.med,base=2), y = log(shift.summary$abun.init.med,base=2) )
(d <- ggplot( shift.summary, aes(x=abun.init.med/30,y=log(abun.shifts.med,base=2))) + 
    geom_hline( yintercept=0, lty=2 ) +
    geom_smooth(method = 'lm', se = T, col='black' ) +
    geom_point(size=3, pch=1, col='slateblue') + 
    ylab("Cover shift") + xlab("Appox. initial cover (%)") +
    scale_y_continuous( breaks=c(sqrt(10),1,0,-1,-sqrt(10)), 
                        labels=c('10x','2x','0','1/2x','1/10x')) +
    scale_x_continuous(trans = "log2") +
    theme_classic() )

cowplot::plot_grid( a, c, b, d, ncol=2, align = 'hv', labels = "AUTO" )
ggsave(file="R/Figs/abundance+peak_shift_intial.svg",width = 6, height = 6)

head(shift.summary %>% arrange(-elev.init.med))
head(shift.summary %>% arrange(elev.shifts.med),20)
head(shift.summary %>% arrange(-abun.init.med))
head(shift.summary %>% arrange(abun.shifts.med))
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
(bpa <- ggplot( shift.summary, aes(y=elev.shifts.med,x=1) ) +
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
(bpb <- ggplot( shift.summary, aes(y=log(abun.shifts.med,base=2),x=1) ) +
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
compare_all_plot <- shift.summary
compare_all_plot$labels <- factor( compare_all_plot$taxon, labels=1:47 )
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

cor.test( compare_all_plot_fun$elev.shifts.med, log(compare_all_plot_fun$abun.shifts.med,base =2 ))
(xy <- ggplot( compare_all_plot_fun, aes(x=log(abun.shifts.med,base=2),y=elev.shifts.med)) + 
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point( aes(fill=group), pch=21,size=3) +
    # geom_text_repel(aes(label=labels),
                    # box.padding = 0, point.padding = 0, size=2) +
    # geom_text(aes(label=labels),size=3, nudge_y = 10) +
    theme_classic() +
    scale_x_continuous(breaks=c(log(50,base=2),log(10,base=2),log(5,base=2),log(2,base=2),0,log(0.5,base=2),log(1/5,base=2),log(1/10,base=2),log(1/50,base=2)),
                       labels=c('50x','10x','5x','2x','0','1/2x','1/5x','1/10x','1/50x'),
                       position="bottom") +
    scale_fill_manual( values=(c("white","darkred", "red","pink", "darkgrey", "#996633")), guide='none' ) +
    xlab("Abundance shift") + ylab("Elevation shift (cm)"))

# a densities not violin plots
ydens <- axis_canvas(xy, axis = "y", coord_flip = TRUE)+
  geom_vline(xintercept=mean(compare_all_plot$elev.shifts.med), col='red' ) +
  geom_vline(xintercept=0) +
  geom_density(data = compare_all_plot, aes(x = elev.shifts.med),
               alpha = 0.7, size = 0.5, outline.type = "full") +
  coord_flip()

# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
xdens <- axis_canvas(xy, axis = "x")+
  geom_vline(xintercept=mean(log(compare_all_plot$abun.shifts.med,base=2)), col='red' ) +
  geom_vline(xintercept=0) +
  geom_density(data = compare_all_plot, aes(x = log(abun.shifts.med,base=2) ),
               alpha = 0.7, size = 0.5, outline.type = "full")

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

p1 <- insert_xaxis_grob(xy, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
ggsave(file="R/Figs/abundance~peak.svg",width = 3.5, height = 3.5)

write_csv( compare_all_plot_fun, "R/output/shifts_predicted.csv")
#

 

#
































