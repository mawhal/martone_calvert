###This script is to compare iButton data from each site
#Load in required packages
library(tidyverse)
library(lubridate)

#Load in and reformat iButton data
FB.high<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FB_High_N_BlackPVCCap.csv", skip=14)
FB.high$`Date/Time`<-as.POSIXct(FB.high$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
FB.low1<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FB_Low_S_BlackPVCCap.csv", skip=14)
FB.low1$`Date/Time`<-as.POSIXct(FB.low1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
FB.low2<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FB_Low_S_BlackPVCCap2.csv", skip=14)
FB.low2$`Date/Time`<-as.POSIXct(FB.low2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
FB.low3<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_FBB_Low_N_BlackPVCCap.csv", skip=14)
FB.low3$`Date/Time`<-as.POSIXct(FB.low3$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.high1<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_High_N_BlackPVCCap.csv", skip=14)
NB.high1$`Date/Time`<-as.POSIXct(NB.high1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.high2<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_High_N_BlackPVCCap2.csv", skip=14)
NB.high2$`Date/Time`<-as.POSIXct(NB.high2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.low1<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_Low_N_BlackPVCCap.csv", skip=14)
NB.low1$`Date/Time`<-as.POSIXct(NB.low1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
NB.low2<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_NBB_Low_N_BlackPVCCap2.csv", skip=14)
NB.low2$`Date/Time`<-as.POSIXct(NB.low2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.low<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_Low_N_BlackPVCCapE.csv", skip=14)
WB.low$`Date/Time`<-as.POSIXct(WB.low$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.high1<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_High_N_BlackPVCCapW.csv", skip=14)
WB.high1$`Date/Time`<-as.POSIXct(WB.high1$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.high2<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_High_N_BlackPVCCapW.csv", skip=14)
WB.high2$`Date/Time`<-as.POSIXct(WB.high2$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")
WB.high3<-read_csv("./Data/environmetal_data/iButtons/ibuttons/2016 ibuttons/2016 ibuttons/2016_Hakai_WBB_High_N_BlackPVCCapW.csv", skip=14)
WB.high3$`Date/Time`<-as.POSIXct(WB.high3$`Date/Time`,format="%m/%d/%Y %H:%M:%OS")

##FB.low1 and NB.low2 have fewer datapoints than they should. Consider this when combining

FB.LOW<-FB.low1
FB.LOW$Value2<-FB.low2$Value

lines(FB.high$Value~FB.high$`Date/Time`,cex=0.4, col="blue")
lines(NB.high1$Value~FB.high$`Date/Time`,cex=0.4, col="red")
points(WB.high3$Value~FB.high$`Date/Time`,cex=0.4, col="lightgreen")
plot(FB.high$Value~WB.high3$Value)
abline(0,1)

##For high zone use FB.high, NB.high1, WB.high1 (all 6 am start)
##For low zone, use FB.low2 (4am), NB.low1 (6am), WB.low (6am), 

FB.high.summ<-FB.high %>%
  group_by(month=floor_date(`Date/Time`, "month")) %>%
  summarize(above_20=sum(Value>19.9&Value<24.9), above_25=sum(Value>24.9&Value<29.9), above_30=sum(Value>29.9), below_5=sum(Value<5&Value>0), below_0=sum(Value<0), site = "FifthBeach.high")

NB.high1.summ<-NB.high1 %>%
  group_by(month=floor_date(`Date/Time`, "month")) %>%
  summarize(above_20=sum(Value>19.9&Value<24.9), above_25=sum(Value>24.9&Value<29.9), above_30=sum(Value>29.9), below_5=sum(Value<5&Value>0), below_0=sum(Value<0), site = "NorthBeach.high")

WB.high1.summ<-WB.high1 %>%
  group_by(month=floor_date(`Date/Time`, "month")) %>%
  summarize(above_20=sum(Value>19.9&Value<24.9), above_25=sum(Value>24.9&Value<29.9), above_30=sum(Value>29.9), below_5=sum(Value<5&Value>0), below_0=sum(Value<0), site = "WestBeach.high")

FB.low1.summ<-FB.low1 %>%
  group_by(month=floor_date(`Date/Time`, "month")) %>%
  summarize(above_20=sum(Value>19.9&Value<24.9), above_25=sum(Value>24.9&Value<29.9), above_30=sum(Value>29.9), below_5=sum(Value<5&Value>0), below_0=sum(Value<0), site = "FifthBeach.low")

NB.low1.summ<-NB.low1 %>%
  group_by(month=floor_date(`Date/Time`, "month")) %>%
  summarize(above_20=sum(Value>19.9&Value<24.9), above_25=sum(Value>24.9&Value<29.9), above_30=sum(Value>29.9), below_5=sum(Value<5&Value>0), below_0=sum(Value<0), site = "NorthBeach.low")

WB.low.summ<-WB.low %>%
  group_by(month=floor_date(`Date/Time`, "month")) %>%
  summarize(above_20=sum(Value>19.9&Value<24.9), above_25=sum(Value>24.9&Value<29.9), above_30=sum(Value>29.9), below_5=sum(Value<5&Value>0), below_0=sum(Value<0), site = "WestBeach.low")

iBut<-rbind(FB.high.summ, NB.high1.summ,WB.high1.summ,FB.low1.summ, NB.low1.summ,WB.low.summ )
iBut_sum<-iBut %>% 
  group_by(site) %>%
  summarize(above_20 = sum(above_20),above_25 = sum(above_25), above_30 = sum(above_30), below_5 = -1*sum(below_5), below_0 = -1*sum(below_0) )
  


iBut_gg<-cbind(gather(iBut_sum, key = "temp", value = "count", 2:6))
iBut_gg$site<-factor(iBut_gg$site, levels=c("FifthBeach.low", "NorthBeach.low", "WestBeach.low", "FifthBeach.high", "NorthBeach.high", "WestBeach.high"))
iBut_gg$temp<-factor(iBut_gg$temp, levels=c("above_30","above_25","above_20",  "below_0", "below_5"))

ggplot(iBut_gg, aes(x = site, y = count), position="stack")+
  geom_col(aes(fill=temp))+
  ylim(c(-200,200))+
  theme_cowplot(15)+
  scale_fill_manual(values = c("darkred", "red","pink", "navy", "blue"))
  





   