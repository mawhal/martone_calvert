#download and start readxl 
install.packages('readxl')
library("readxl")

#Summon and view all of the sheets within the 2012_Hakai_R document 


WBlow2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls")
WBlow2012

WBmid2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet= "WB_mid_R")
WBmid2012

WBhigh2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet= "WB_high_R")
WBhigh2012

NBlow2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet= "NB_low_R")
NBlow2012

NBmid2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet="NB_mid_R")
NBmid2012

NBhigh2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet="NB_high_R")
NBhigh2012

FBlow2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet="FB_low_R")
FBlow2012

FBmid2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet="FB_mid_R")
FBmid2012

FBhigh2012 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2012/2012_Hakai_R.xls", sheet="FB_high_R")
FBhigh2012

#Combine and stack all of the data according to column names 

data2012 <-rbind(WBlow2012, WBmid2012, WBhigh2012, NBlow2012, NBmid2012, NBhigh2012, FBlow2012, FBmid2012, FBhigh2012)
View(data2012)
