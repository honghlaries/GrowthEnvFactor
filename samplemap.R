## ======================cleaning==========================================
library(ggplot2)
library(maptools)
library(grid)
library(ggmap)
library(dplyr)
library(tidyr)
rm(list=ls())

## ===========================data readln==================================
##========================================================================
cname <- c("Longitude","Latitude")

coo.land <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Land.kml", ignoreAltitude = T)[[1]])
colnames(coo.land) <- cname

coo.sm1 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh1.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm1) <- cname

coo.sm2 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh2.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm2) <- cname

coo.sm3 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh3.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm3) <- cname

coo.sm4 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh4.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm4) <- cname

coo.sm5 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh5.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm5) <- cname

coo.sm6 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh6.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm6) <- cname

coo.sm7 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh7.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm7) <- cname

coo.sw1 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall1.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw1) <- cname

coo.sw2 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall2.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw2) <- cname

coo.sw3 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall3.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw3) <- cname

coo.sw4 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall4.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw4) <- cname

coo.site <- read.csv("Data/map/sitecor.csv")
coo.site$col <- as.character(coo.site$col)


##================= mapping============================
##============================================================
ggplot()+
  geom_polygon(data = coo.land,
               aes(x = Longitude, y = Latitude, fill = "#B7B7B7"))+
  geom_polygon(data = coo.sm1,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm2,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm3,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm4,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm5,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm6,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm7,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sw1,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw2,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw3,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw4,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_point(data = coo.site %>% filter(group != "NV"), 
             aes(x = Longitude, y = Latitude, col = col ,size = 10, shape = group))+
  coord_equal(xlim = c(121.2475, 121.2725), ylim = c(32.4475, 32.4725))+
  scale_x_continuous(breaks = 121.25 + 0:4 * 0.005, 
                     labels = paste("121°",c(15 + 0:4 * 0.3),"'E",sep=""))+
  geom_vline(xintercept = 121.25 + 0:4 * 0.005, alpha = 0.5)+
  scale_y_continuous(breaks = 32.45 + 0:4 * 0.005, 
                     labels = paste("32°",c(27 + 0:4 * 0.3),"'W",sep=""))+
  geom_hline(yintercept = 32.45 + 0:4 * 0.005, alpha = 0.5) +
  scale_size_identity()+
  scale_fill_identity() +
  scale_colour_identity()+
  scale_shape_manual(values = c(16,17,15))+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(),
        panel.background = element_rect(fill = "#B9D3EE"),
        panel.grid = element_blank(),
        legend.position = "none") 
ggsave(filename = "map/line_1.png", dpi = 600, width = 6, height = 6)


ggplot()+
  geom_polygon(data = coo.land,
               aes(x = Longitude, y = Latitude, fill = "#B7B7B7"))+
  geom_polygon(data = coo.sm1,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm2,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm3,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm4,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm5,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm6,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm7,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sw1,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw2,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw3,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw4,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_point(data = coo.site %>% filter(group != "NV"), 
             aes(x = Longitude, y = Latitude, col = col ,size = 10, shape = group))+
  coord_equal(xlim = c(121.41, 121.425), ylim = c(32.3625, 32.3775))+
  scale_x_continuous(breaks = 121.4125 + 0:2 * 0.005, 
                     labels = paste("121°",(0.4125 + 0:2 * 0.005) * 60 ,"'E",sep=""))+
  geom_vline(xintercept = 121.4125 + 0:2 * 0.005, alpha = 0.5) +
  scale_y_continuous(breaks = 32.365 + 0:2 * 0.005, 
                     labels = paste("32°",(0.365 + 0:2 * 0.005) * 60,"'W",sep=""))+
  geom_hline(yintercept = 32.365 + 0:2 * 0.005, alpha = 0.5) +
  scale_size_identity()+
  scale_fill_identity() +
  scale_colour_identity()+
  scale_shape_manual(values = c(16,17,15))+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(),
        panel.background = element_rect(fill = "#B9D3EE"),
        panel.grid = element_blank(),
        legend.position = "none") 
ggsave(filename = "map/line_2.png", dpi = 600, width = 6, height = 6)


ggplot()+
  geom_polygon(data = coo.land,
               aes(x = Longitude, y = Latitude, fill = "#B7B7B7"))+
  geom_polygon(data = coo.sm1,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm2,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm3,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm4,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm5,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm6,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm7,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sw1,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_polygon(data = coo.sw2,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_polygon(data = coo.sw3,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_polygon(data = coo.sw4,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_point(data = coo.site %>% filter(group != "NV"), 
             aes(x = Longitude, y = Latitude, col = col ,size = 2, shape = group))+
  coord_equal(ratio = 1, xlim = c(121.23, 121.45), ylim = c(32.35, 32.48))+
  scale_x_continuous(breaks = 121.235 + 0:6 * 0.035, 
                     labels = paste("121°",(0.235 + 0:6 * 0.035) * 60 ,"'E",sep=""))+
  #geom_vline(xintercept = 121.235 + 0:6 * 0.035, alpha = 0.5) +
  scale_y_continuous(breaks = 32.355 + 0:3 * 0.040, 
                     labels = paste("32°",(0.355 + 0:3 * 0.040) * 60,"'W",sep=""))+
  #geom_hline(yintercept = 32.355 + 0:3 * 0.040, alpha = 0.5) +
  scale_size_identity()+
  scale_fill_identity() +
  scale_colour_identity()+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(),
        panel.background = element_rect(fill = "#B9D3EE"),
        panel.grid = element_blank(),
        legend.position = "bottom")  
ggsave(filename = "map/base.png", dpi = 600, width = 12, height = 8)