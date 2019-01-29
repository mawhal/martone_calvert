#download and start readxl 
install.packages('readxl')
library("readxl")

#Summon and view all of the sheets within the 2013_Hakai_R document 


WBlow2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet="WB_low_R")
WBlow2014

WBmid2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet= "WB_mid_R")
WBmid2014

WBhigh2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet= "WB_high_R")
WBhigh2014

NBlow2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet= "NB_low_R")
NBlow2014

NBmid2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet="NB_mid_R")
NBmid2014

NBhigh2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet="NB_high_R")
NBhigh2014

FBlow2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet="FB_low_R")
FBlow2014

FBmid2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet="FB_mid_R")
FBmid2014

FBhigh2014 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2014/2014_Hakai_R.xlsx", sheet="FB_high_R")
FBhigh2014

#Combine and stack all of the data according to column names 

data2014 <-rbind(WBlow2014, WBmid2014, WBhigh2014, NBlow2014, NBmid2014, NBhigh2014, FBlow2014, FBmid2014, FBhigh2014)
View(data2014)
