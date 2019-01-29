#download and start readxl 
install.packages('readxl')
library("readxl")

#Summon and view all of the sheets within the 2013_Hakai_R document 


WBlow2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="WB_low_R")
WBlow2017

WBmid2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet= "WB_mid_R")
WBmid2017

WBhigh2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet= "WB_high_R")
WBhigh2017

NBlow2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet= "NB_low_R")
NBlow2017

NBmid2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="NB_mid_R")
NBmid2017

NBhigh2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="NB_high_R")
NBhigh2017

FBlow2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="FB_low_R")
FBlow2017

FBmid2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="FB_mid_R")
FBmid2017

FBhigh2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="FB_high_R")
FBhigh2017

MClow2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="MC_low_R")
MClow2017

MCmid2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="MC_mid_R")
MCmid2017

MChigh2017 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2017/2017_Hakai_R.xlsx", sheet="MC_high_R")
MChigh2017

#Combine and stack all of the data according to column names 

data2017 <-rbind(WBlow2017, WBmid2017, WBhigh2017, NBlow2017, NBmid2017, NBhigh2017, FBlow2017, FBmid2017, FBhigh2017, MClow2017, MCmid2017, MChigh2017)
View(data2017)
