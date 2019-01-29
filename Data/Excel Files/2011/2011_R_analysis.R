#download and start readxl 
install.packages('readxl')
library("readxl")

#Summon and view all of the sheets within the 2011_Hakai_R document 

WBlow2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="WB_low_R")
WBlow2011

WBmid2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="WB_mid_R")
WBmid2011

WBhigh2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="WB_high_R")
WBhigh2011

NBlow2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="NB_low_R")
NBlow2011

NBmid2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="NB_mid_R")
NBmid2011

NBhigh2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="NB_high_R")
NBhigh2011

FBlow2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="FB_low_R")
FBlow2011

FBmid2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="FB_mid_R")
FBmid2011

FBhigh2011 <- read_excel("~/Documents/University/2016-2017/Martone Data/JJ HAKAI/2011/2011_Hakai_R.xlsx", sheet="FB_high_R")
FBhigh2011

#Combine and stack all of the data according to column names 

data2011 <-rbind(WBlow2011, WBmid2011, WBhigh2011, NBlow2011, NBmid2011, NBhigh2011, FBlow2011, FBmid2011, FBhigh2011)
View(data2011)
