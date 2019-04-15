# Leaflet

Leaflet 是一个 R 地图应用包，可以用于病原流行病学分析。

## 1. 安装

```R
# https://cloud.tencent.com/developer/article/1092335
library(plyr)
library(mapdata)
library(leaflet)
library(maptools)
library(ggplot2)

#导入中国各省会城市地理信息数据：
province_city <- read.csv("c:/rstudy/chinaprovincecity.csv")
province_city$size<-round(runif(34,5,10),2)
province_city$type<-factor(sample(LETTERS[1:5],34,replace=TRUE))
co<-substr(rainbow(34),1,7)
province_city<-data.frame(province_city,co)

leaflet(province_city)%>%addTiles()%>%setView(lng=116.38,lat=39.9,zoom=3)%>%addMarkers(lng=~jd,lat=~wd,popup=~city)
```

## 2. 使用
