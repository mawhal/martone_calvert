# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses community trajectory analysis on Martone seaweed community data

# load libraries
library( tidyverse )
library( vegclust )   # CTA, length, angles, directionality, projection, distances between segments/trajectories
library( adespatial ) # dynamic-based beta diversity, local contributions to beta diversity
library( smacof )     # stress minimization using majorization (smacof) -- new approaches for MDS
library( vegan )
library( RColorBrewer )
library( BiodiversityR )
#
# custom color scheme (three shades of three colors) for plotting
lighten <- function(color, factor=1.5){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}
darken <- function(color, factor=2){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}
base <- c('goldenrod', '#424F9F', 'black')
cols = as.vector(matrix( c(sapply( base, darken ), sapply( base, lighten ), base), ncol=3, byrow = TRUE ))


# To specify community dynamics, we need three data items:
#   
# A set of community states (i.e. coordinates in a space Ω), described using a distance matrix d;
# A vector specifying the site (i.e. sampling unit) corresponding to each community state;
# A vector specifying the survey (i.e. time point) corresponding to the sampling of each community state.
# 

## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )


## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
muse <- am[ am$Year != "2011", ]
# remove Meay Channel
## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
muse <- muse[ muse$Site != "Meay Channel", ]
# Only use Mid-shore transects for now
# muse <- muse[ muse$Zone == "MID", ]
# muse <- muse[ muse$Site == "North Beach", ]
muse <- droplevels(muse)

# restrict rows of d to ones with UID's in the metadata
duse <- ad[ ad$UID %in% muse$UID, ]
dm <- left_join( duse, muse )


# for now, restrict community analysis to algae only
d <- dm %>% 
  # filter( non.alga.flag %in% c("Algae"  )
  filter( motile_sessile == "sessile" )

# calculate mean abundance per transect in each year
dmean <- d %>% 
  group_by( Year, Site, Zone, taxon_lumped ) %>%
  summarise( Abundance=mean(Abundance) )

# spread out
d.comm.mean <- dmean %>%
  spread( taxon_lumped, Abundance, fill=0 )

# remove UID column from community data
meta <- d.comm.mean[ ,1:3 ]
comm <- as.matrix(d.comm.mean[,-c(1:3)])

# define site as the particular trasect
Zone <- factor( meta$Zone, levels=c("LOW","MID","HIGH"), labels=c("Low","Mid","High") )
Site <- factor( meta$Site, labels=c("5","N","W") )
site <- paste( Site, Zone, sep="." )
site <- factor( site, levels= c("5.Low","5.Mid","5.High","N.Low","N.Mid","N.High","W.Low","W.Mid","W.High"),
                labels= c("5.L","5.M","5.H","N.L","N.M","N.H","W.L","W.M","W.H"))

# define year
year <- meta$Year
# define 2016 for plotting later
year2 <- as.character(year)
year2[year2!=2016] <- ""

# x<-radfit(cs)
# plot(x)

# which species are exceedingly rare?
sort(colSums(comm))
hist(colSums(comm),breaks = seq(0,1600,2))
# which species are most abundant through the time series?
dominance <- apply( comm, 1, function(z) names(z)[order(z,decreasing = T)]  )
dominance[1:3,]

# quick beta diversity by group
z <- betadiver(comm,"z")
mod <- betadisper(z, site)
plot(mod)
plot(mod, axes=c(3,1) )
boxplot(mod)
anova(mod)
(mod3B <- betadisper(z, site, type = "median", bias.adjust=TRUE))
anova(mod3B)
permutest(mod3B, permutations = 99)
plot(mod3B)
boxplot(mod3B)

adonis2( comm~site*year, by='margin'  )

# Trajectory analysis
# calculate a distance 
D_man <- vegdist( comm, method="manhattan", transform = function(x) log(x+1) ) # uses log abundances, so will downplay importance of rarest taxa
D_bray <- vegdist( comm, method="bray" )
D_bray_root <- vegdist( comm^0.5, method="bray" )
D_bray_root2 <- vegdist( comm^0.25, method="bray" )
D_use <- D_bray_root


# display trjectories in PCoA - uses cmdscale (MDS/PCoA)
# windows(5,5)
par(mar=c(5,4,1,1)+0.1,lty=1)
# new scheme for colors and tide heights
cols2 <- c("goldenrod","#646FB6","black")
fills <- c("goldenrod","#646FB6","white")
ltys  <- c(3,2,1)
x <- trajectoryPCoA( D_use,  site, year,
                     traj.colors = cols2[as.numeric(Site)],
                     axes=c(1,2), length=0, lwd=ltys[as.numeric(Zone)] )
sites <- site
surveys <- year
traj.colors = cols2[as.numeric(Site)]
siteIDs = unique(sites)
nsite = length(siteIDs)
selection = 1:nsite
selIDs = siteIDs[selection]
D2 = as.dist(as.matrix(D_use))#[sites %in% selIDs, sites %in% selIDs])
cmd_D2 <- cmdscale(D2, eig = TRUE, add = TRUE, k = nrow(as.matrix(D2)) - 
                     1)
cmd_D2 <- add.spec.scores(cmd_D2, comm^0.5, 
                          method='pcoa.scores', Rscale=T, scaling=3, multi=0.1)
axes=c(1,2)
x <- cmd_D2$points[, axes[1]]
y <- cmd_D2$points[, axes[2]]
spec1 <- cmd_D2$cproj[, axes[1]]
spec2 <- cmd_D2$cproj[, axes[2]]

# reduce number of taxa in species scores
sel1 <- names(spec1[abs(spec1)>2.5*sd(spec1, na.rm=T)])
sel2 <- names(spec2[abs(spec2)>2.5*sd(spec2, na.rm=T)])

spec.sel <- unique(c(sel1,sel2))
spec1.sel <- spec1[names(spec1)%in%spec.sel]
spec2.sel <- spec2[names(spec2)%in%spec.sel]
specs <- data.frame( x=spec1.sel, y=spec2.sel, taxon=names(spec1.sel) )
specs$taxon <- vegan::make.cepnames(specs$taxon)
windows(4,4)
par(mar=c(5,4,2,2)+0.1,xpd=T)
plot( y~x, data=specs, type = "p", asp = 1,
     xlab = paste0("PCoA ", axes[1],
                   " (", round(100 * cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)),
                   "%)"), 
     ylab = paste0("PCoA ", axes[2], 
                   " (", round(100 * cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)), "%)"),
     xlim=c(-0.3,0.3), ylim=c(-0.2, 0.25))
text( y~x, data=specs, labels=specs$taxon,  cex=0.8,
      pos=c(1,1,3,2,3,1,2,1,4,4,3,2,2,1,4))
# sitesred = sites[sites %in% selIDs]
# surveysred = surveys[sites %in% selIDs]
# 
# for (i in 1:length(selIDs)) {
#   ind_surv = which(sitesred == selIDs[i])
#   # define zone and site for each instance
#   li = ltys[unique(as.numeric(Zone[ind_surv]))]
#   ci = cols2[unique(as.numeric(Site[ind_surv]))]
#     ind_surv = ind_surv[order(surveysred[sitesred ==
#                                            selIDs[i]])]
#   for (t in 1:(length(ind_surv) - 1)) {
#     niini = ind_surv[t]
#     nifin = ind_surv[t + 1]
#       arrows(x[niini], y[niini], x[nifin], y[nifin],
#              col = ci, lty=li, lwd=2, length=0 )
# 
#   }
# }
# points( x,y, pch=16, bg='white', col=cols2[as.numeric(Site)],cex=0.8 )
# points( x[1:9],y[1:9], pch=21, bg=fills[as.numeric(Site)], cex=1.5 )
# # points( x[which(surveys==2016)],y[which(surveys==2016)], pch=25, bg=fills[as.numeric(Site)], cex=1.5 )
# points( x[(length(x)-8):length(x)],y[(length(x)-8):length(x)], pch=22, bg=fills[as.numeric(Site)], cex=1.5 )
# # text( x,y,labels = year-2012, adj =c(0.4,0.4), col='black', cex=0.6 )
# legend("topleft",    bty="o", legend = c("Fifth","North","West"), title = "Shore",
#        col = cols2, lwd=2 )
# legend("bottomleft", bty="o", legend = c("High","Mid","Low"), title="Zone",
#        col = 'black', lty=c(1,2,3),lwd=2 )
# legend("topright",   bty="o", legend = c("2012","other","2019"), title="Year",
#        pch=c(21,16,22),  col = 'black', bg='white', lty=0 )



## wrap into a single data.frame
trajdf <- data.frame( x, y, transect=site, site=Site, zone=Zone,  year )
trajdf$site <- as.character( trajdf$site )
trajdf$site[ trajdf$site=="5" ] <- "Fifth Beach"
trajdf$site[ trajdf$site=="N" ] <- "North Beach"
trajdf$site[ trajdf$site=="W" ] <- "West Beach"
trajdf$site <- factor( trajdf$site )
trajdf$site <- relevel(trajdf$site, ref="West Beach")

windows(4,6)
ggplot( trajdf, aes(x=x,y=y,col=zone, group=transect)) + facet_wrap(~site,ncol=1) +
  geom_point(size=1) + geom_path() +
  geom_point( data=filter(trajdf,year==2012), aes(x=x,y=y,fill=zone), size=3, pch=21, col='black' ) +
  geom_point( data=filter(trajdf,year==2019), aes(x=x,y=y,fill=zone), size=3, pch=22, col='black' ) +
  # geom_point( data=specs, aes(x,y, group=1), col="black",size=3) +
  xlab( paste0("PCoA ", axes[1], " (", round(100 * cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)),"%)")) +
  ylab( paste0("PCoA ", axes[2], " (", round(100 * cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)),"%)"))+
  scale_color_manual(values=c("black","grey50","grey75")) + 
  scale_fill_manual(values=c("black","grey50","grey75")) + 
  theme_minimal()  

ggsave( 'R Code and Analysis/Figs/PCoA_4root_ggplot.pdf', width=4, height=6 )




# centered trajectories
x<-trajectoryPCoA( centerTrajectories(D_use, as.numeric(factor(site)) ),  as.numeric(factor(site)), year,
                traj.colors = cols[as.numeric(site)], 
                axes=c(1,2), length=0, lwd=1, lty=0 )
points( x$points[1:9,1:2], pch=21, bg=cols[as.numeric(site)], cex=2 )
points( x$points[(nrow(x$points)-8):nrow(x$points),1:2], pch=22, bg=cols[as.numeric(site)], cex=2 )
arrows( x$points[1:9,1], x$points[1:9,2], 
        x$points[(nrow(x$points)-8):nrow(x$points),1],
        x$points[(nrow(x$points)-8):nrow(x$points),2],
        col=cols[as.numeric(site)], length = 0.1, lwd=2 )
# text( x$points[,1:2],labels = year2, pos = 1, col=rep(cols,each=length(unique(year))) )
legend("bottomright", bty="n", legend = levels(site), col = cols, lwd=3 )


## try to run a PCoA from which commcies scores can be extracted
# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix of commcies
comm.bray <- vegdist(comm)
comm.b.pcoa <- cmdscale(comm.bray, k=(nrow(comm)-1), eig=TRUE)
# Plot of the sites
dev.new(title="PCoA on commcies - Percentage difference")
ordiplot(scores(comm.b.pcoa, choices=c(1,2)), type="t", main="PCoA with commcies weighted averages")
abline(h=0, lty=3)
abline(v=0, lty=3)
# Add weighted average projection of commcies
comm.wa <- wascores(comm.b.pcoa$points[,1:2], comm)
text(comm.wa, rownames(comm.wa), cex=0.7, col="red")

# PCoA and projection of commcies vectors using function pcoa() ####
comm.h <- decostand(comm, "hellinger")
comm.h.pcoa <- pcoa(dist(comm.h))
# Biplots
dev.new(title="PCoA with commcies vectors", width=14, height=8)
par(mfrow=c(1,2))
# First biplot: Hellinger-transformed commcies data
biplot.pcoa(comm.h.pcoa, comm.h, dir.axis1=-1)
abline(h=0, lty=3)
abline(v=0, lty=3)
# Second biplot: standardized Hellinger-transformed commcies data
comm.std <- scale(comm.h)
biplot.pcoa(comm.h.pcoa, comm.std, dir.axis1=-1)
abline(h=0, lty=3)
abline(v=0, lty=3)


# can also use MDS to represent trajectories
# stress plot
mMDS  <- mds( D_use, ndim=8 )
mMDS
mMDS2 <- mds( centerTrajectories(D_man, as.numeric(factor(site)) ) ) 
mMDS2
use <- mMDS$conf
trajectoryPlot( use,  as.numeric(factor(site)), year,
                traj.colors = cols, 
                axes=c(1,2), length=0.1, lwd=2 )
text( use,labels = year2, pos = 1, col=rep(cols,each=length(unique(year))) )
legend("topright", bty="n", legend = unique(site), col = cols, lwd=2)




## Trajectory statistics
trajectoryLengths(        D_use, site, year ) 
trajectoryAngles(         D_use, site, year )
trajectoryAngles(         D_use, site, year, all=TRUE ) # high degree of angle homogeneity
trajectoryDirectionality( D_use, site, year ) # despite longer segements in LOW, MID and HIGH often are more directional
plot( trajectoryDirectionality( D_use, site, year ), cex=3, pch=21, bg=cols ) # similar directionality among sites
trajectoryProjection(     D_use, 1,2:6 )
trajectoryConvergence(    D_use, site, year, symmetric = FALSE )
trajectoryPlot( x, site, year, axes=1:2  )

Ds = segmentDistances( D_use, site, year )$Dseg
Ds
mMDS = mds(Ds)
mMDS

xret = mMDS$conf
par(mar=c(4,4,1,1))
plot(xret, xlab="axis 1", ylab = "axis 2", asp=1, pch=21,
     bg=rep(cols, each=7))#, 
     # xlim=c(-1.5,1), ylim=c(-1,1.5))
text(xret, labels=rep(paste0("s",1:7),9), pos=1)
legend("topleft", pt.bg=cols, pch=21, bty="n", legend=paste0("trajectory",1:9))


trajectoryDistances( D_use, site, year, distance.type = "Hausdorff")

                    