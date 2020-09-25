##Functional groups through time
##This code imports the community matrix, collapses by functional groups and makes some plots
##Code written by Sam Starko
library(tidyverse)
library(Rmisc)

#Import functional groups by species
FunGroups<-read_csv("./Data/taxa/Algae_functional_groups.csv")
FunGroups$taxon<-gsub(" ",".",FunGroups$taxon)


#Import community matrix
comm<-read_csv("./R Code and Analysis/output from r/community.csv")
#Check which species names don't match #These should all be animals and other non-algal fields
colnames(comm)[colnames(comm) %in% FunGroups$taxon == "FALSE"] 
#Remove non-algal fields from community matrix
comm2<-comm[,colnames(comm) %in% FunGroups$taxon]

#Creat matrix of functional groups summed
taxon<-data.frame(taxon = colnames(comm2))
taxon.key<-left_join(taxon, FunGroups, by = "taxon") 
colnames(comm2)<-taxon.key$funct_Sep2020
comm3<-comm2[,which(colnames(comm2) != "NA")]
comm.sum <- t(rowsum(t(comm3), group = colnames(comm3), na.rm = T)) %>% as_tibble()
comm.sum$Year<-comm$Year
comm.sum$Site<-comm$Site
comm.sum$Zone<-comm$Zone

blade<-summarySE(data = comm.sum, measurevar = "blade", groupvars = c(as.character("Year"), "Site", "Zone"))
colnames(blade)[5]<-"mean"
blade$FunGroup<-"blade"
crust<-summarySE(data = comm.sum, measurevar = "crust", groupvars = c(as.character("Year"), "Site", "Zone"))
colnames(crust)[5]<-"mean"
crust$FunGroup<-"crust"
thin_turf<-summarySE(data = comm.sum, measurevar = "thin_turf", groupvars = c(as.character("Year"), "Site", "Zone"))
colnames(thin_turf)[5]<-"mean"
thin_turf$FunGroup<-"thin_turf"
kelp<-summarySE(data = comm.sum, measurevar = "kelp", groupvars = c(as.character("Year"), "Site", "Zone"))
colnames(kelp)[5]<-"mean"
kelp$FunGroup<-"kelp"
seagrass<-summarySE(data = comm.sum, measurevar = "seagrass", groupvars = c(as.character("Year"), "Site", "Zone"))
colnames(seagrass)[5]<-"mean"
seagrass$FunGroup<-"seagrass"
turf<-summarySE(data = comm.sum, measurevar = "turf", groupvars = c(as.character("Year"), "Site", "Zone"))
colnames(turf)[5]<-"mean"
turf$FunGroup<-"turf"

d<-rbind(blade, crust, thin_turf, kelp, seagrass, turf) %>% as_tibble()
d$Zone<-factor(d$Zone, levels = c("LOW", "MID", "HIGH"))
d$FunGroup<-factor(d$FunGroup, levels = rev(c("turf", "thin_turf", "blade", "crust", "seagrass", "kelp")))

ggplot(d, aes(x = Year, y = mean, fill = FunGroup))+
  geom_bar(position="stack", stat="identity")+
  theme(plot.title = element_text(size = 5))+
  theme_cowplot()+
  facet_wrap(~Site + Zone)+
  scale_fill_manual(values = c("darkred", "red","pink", "darkgrey", "darkgreen", "#996633") %>% rev())
  

#filter(d, Site == "Fifth Beach"&Zone=="LOW")
