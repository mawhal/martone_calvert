# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses HMSC to model communities change over time and space

## useful references for this code
citation( "Hmsc" )
# https://github.com/hmsc-r/HMSC

## directories
localDir = "."
# dataDir = file.path(localDir, "data")
ModelDir = file.path( localDir, "R Code and Analysis/models" )
MixingDir = file.path(localDir, "R Code and Analysis/mixing")

## load libraries
library( tidyverse )
library( Hmsc )
library( corrplot )
library( viridis )


## load the model
list.files( ModelDir )
model = 2
mload <- load( paste(ModelDir,list.files( ModelDir ), sep="/")[model] )
m


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
windows(5,8)
plotBeta(m, post = postBeta, param = "Sign", supportLevel = 0.95, mar=c(7,11,0,6))

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "shore.elevation", "shore.elev.sq", "temp.anom.sum", "temp.anom.win" ),
                            levels = c("intercept", "shore.elevation", "shore.elev.sq",  "temp.anom.sum", "temp.anom.win" ), 
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 5), levels = colnames(Y)[order(colSums(Y),decreasing = TRUE)], ordered = TRUE)

windows(12,4)
ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")

# pull taxa with positive estimated response to temperature anomaly
taxa.pos.anom <- as.character(pos.neg[ pos.neg$parameter=='temp.anom' & pos.neg$value>0, 'species'])


## variance partitioning
VP = computeVariancePartitioning(m) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
plotVariancePartitioning(m, VP = VP)

VP.df <- as.data.frame(VP$vals) %>% 
  mutate(effect = factor(c("elevation","elev.square","temp.anom.sum", "temp.anom.win","transect","year","site"), 
                         levels = rev(c("elevation","elev.square","temp.anom.sum", "temp.anom.win","transect","year","site")), 
                         ordered = TRUE)) %>% 
  gather(key = species, value = variance, -effect) %>% 
  group_by(species) %>% 
  mutate(tempR2 = variance[effect == "temp.anom.sum"])

hold <- VP.df %>% filter(effect == "temp.anomaly") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, levels = colnames(Y)[order(colSums(Y),decreasing = TRUE)], ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(Y))

windows(8,5)
ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon", viridis(5)), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")





## associations
OmegaCor = computeAssociations(m)
supportLevel = 0.95
# choose the random variable to plot
rlevel = 2
toPlot = ((OmegaCor[[rlevel]]$support>supportLevel) 
          + (OmegaCor[[rlevel]]$support<(1-supportLevel))>0)*OmegaCor[[rlevel]]$mean
# reorder species matrix
plotorder <- order( postBeta$mean[5,], decreasing = TRUE )
toPlot <- toPlot[ plotorder, plotorder]
# reorder automatically
library(lessR)
mynewcor <- corReorder( toPlot, order="hclust", nclusters=4 )
# windows(12,12)
corrplot( toPlot, method = "color", 
         col = colorRampPalette(c("blue","white","red"))(200),
           title = paste("random effect level:", m$rLNames[rlevel]), mar=c(0,0,0.5,0), tl.cex=0.6 )

