###This script is to compare iButton data from each site
#Load in required packages
library(tidyverse)
library(lubridate)

#Load in and reformat iButton data
FB.high<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FB_High_N_BlackPVCCap.csv", skip=14)
FB.high$`Date/Time`<-as.POSIXct(FB.high$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
FB.low1<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FB_Low_S_BlackPVCCap.csv", skip=14)
FB.low1$`Date/Time`<-as.POSIXct(FB.low1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
FB.low2<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FB_Low_S_BlackPVCCap2.csv", skip=14)
FB.low2$`Date/Time`<-as.POSIXct(FB.low2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
FB.low3<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FBB_Low_N_BlackPVCCap.csv", skip=14)
FB.low3$`Date/Time`<-as.POSIXct(FB.low3$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.high1<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_High_N_BlackPVCCap.csv", skip=14)
NB.high1$`Date/Time`<-as.POSIXct(NB.high1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.high2<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_High_N_BlackPVCCap2.csv", skip=14)
NB.high2$`Date/Time`<-as.POSIXct(NB.high2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.low1<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_Low_N_BlackPVCCap.csv", skip=14)
NB.low1$`Date/Time`<-as.POSIXct(NB.low1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.low2<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_Low_N_BlackPVCCap2.csv", skip=14)
NB.low2$`Date/Time`<-as.POSIXct(NB.low2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.low<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_Low_N_BlackPVCCapE.csv", skip=14)
WB.low$`Date/Time`<-as.POSIXct(WB.low$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.high1<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_Low_N_BlackPVCCapE.csv", skip=14)
WB.high1$`Date/Time`<-as.POSIXct(WB.high1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.high2<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_Low_N_BlackPVCCapW.csv", skip=14)
WB.high2$`Date/Time`<-as.POSIXct(WB.high2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.high3<-read_csv("./Data/Environmental Data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_Low_N_BlackPVCCapW2.csv", skip=14)
WB.high3$`Date/Time`<-as.POSIXct(WB.high3$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")

##FB.low1 and NB.low2 have fewer datapoints than they should. Consider this when combining

FB.LOW<-FB.low1
FB.LOW$Value2<-FB.low2$Value

lines(FB.high$Value~FB.high$`Date/Time`,cex=0.4, col="blue")
lines(NB.high1$Value~FB.high$`Date/Time`,cex=0.4, col="red")
points(WB.high3$Value~FB.high$`Date/Time`,cex=0.4, col="lightgreen")
plot(FB.high$Value~WB.high3$Value)
abline(0,1)
