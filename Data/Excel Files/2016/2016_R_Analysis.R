#download and start readxl 
install.packages('readxl')
library("readxl")

#Summon and view all of the sheets within the 2013_Hakai_R document 


WBlow2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="WB_low_R")
WBlow2016

WBmid2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet= "WB_mid_R")
WBmid2016

WBhigh2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet= "WB_high_R")
WBhigh2016

NBlow2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet= "NB_low_R")
NBlow2016

NBmid2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="NB_mid_R")
NBmid2016

NBhigh2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="NB_high_R")
NBhigh2016

FBlow2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="FB_low_R")
FBlow2016

FBmid2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="FB_mid_R")
FBmid2016

FBhigh2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="FB_high_R")
FBhigh2016

MClow2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="MC_low_R")
MClow2016

MCmid2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="MC_mid_R")
MCmid2016

MChigh2016 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2016/2016_Hakai_R.xlsx", sheet="MC_high_R")
MChigh2016

#Combine and stack all of the data according to column names 

data2016 <-rbind(WBlow2016, WBmid2016, WBhigh2016, NBlow2016, NBmid2016, NBhigh2016, FBlow2016, FBmid2016, FBhigh2016, MClow2016, MCmid2016, MChigh2016)
View(data2016)
