#download and start readxl 
install.packages('readxl')
library("readxl")

#Summon and view all of the sheets within the 2013_Hakai_R document 


WBlow2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx")
WBlow2013

WBmid2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet= "WB_mid_R")
WBmid2013

WBhigh2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet= "WB_high_R")
WBhigh2013

NBlow2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet= "NB_low_R")
NBlow2013

NBmid2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet="NB_mid_R")
NBmid2013

NBhigh2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet="NB_high_R")
NBhigh2013

FBlow2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet="FB_low_R")
FBlow2013

FBmid2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet="FB_mid_R")
FBmid2013

FBhigh2013 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2013/2013_Hakai_R.xlsx", sheet="FB_high_R")
FBhigh2013

#Combine and stack all of the data according to column names 

data2013 <-rbind(WBlow2013, WBmid2013, WBhigh2013, NBlow2013, NBmid2013, NBhigh2013, FBlow2013, FBmid2013, FBhigh2013)
View(data2013)
