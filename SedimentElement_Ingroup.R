## ---------------------------------------------------
## init
##
## ---------------------------------------------------
rm(list = ls())
library(dplyr);library(tidyr);library(ggplot2);library(car)

## ---------------------------------------------------
## Functions 
##
## ---------------------------------------------------
datareadln <- function() {
  read.csv("./Data/Result_Sediment.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    dplyr::filter(group != "NV") %>% 
    dplyr::select(N:Pb,SiteID,SplMonth,group,dist) 
}

## ---------------------------------------------------
## comparing the difference inside groups 
##
## ---------------------------------------------------

dat <- datareadln() %>%
  gather(trait, value, N:Pb)


####               ####
####  Gather Plot  ####
####               ####
ggplot(aes(x = SiteID,y = value), data = dat) +
  geom_point(aes(col = group, shape = SplMonth), size = 1) + 
  ylab("Element Content (mg/kg)") + xlab("Site") +
  facet_wrap(~trait, scales = "free") + 
  theme_bw() 
ggsave("sediment/plot/gather_distribution.png",width = 10, height = 8, dpi = 600)

ggplot(aes(x = dist,y = value), data = dat) +
  geom_point(aes(col = group, shape = SplMonth), size = 2) + 
  geom_smooth(aes(x = dist,y = value), data = dat %>% filter(group != "EA"), method = "lm", col = "black") +
  ylab("Element Content (mg/kg)") + xlab("Distance from the edge of community (m)") +
  facet_wrap(~trait, scales = "free") + 
  theme_bw() 
ggsave("sediment/plot/gather_dist.png",width = 10, height = 8, dpi = 600)

####          ####
####    CL    ####
####          ####
dat1 <- dat %>%
  filter(group == "CL")
tags <- unique(dat1$trait)
sites <- unique(dat1$SiteID)
rst <- NULL
for(i in 1: length(tags)) {
  dat2 <- dat1 %>% filter(trait == tags[i])
  aov <- Anova(lm(value~SiteID,data = dat2))
  p <- aov$`Pr(>F)`[1]
  if (p >= 0.05) p <- format(p,digit = 3) else 
    if (p >= 0.01) p <- paste(format(p,digit = 3),"*",sep = " ") else 
      if (p >= 0.001) p <- paste(format(p,digit = 3),"**",sep = " ") else 
         p <- paste("<0.001","***",sep = " ")
  rst <- rbind(rst, data.frame(element = tags[i],p = p))
}
cal <- dat1 %>% group_by(SiteID,trait) %>%
  summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                          format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                          sep = " ± ")) %>%
  spread(SiteID, value) %>%
  rename(element = trait)
rst.tot<- inner_join(cal,rst)

####          ####
####    EA    ####
####          ####
dat1 <- dat %>%
  filter(group == "EA")
tags <- unique(dat1$trait)
sites <- unique(dat1$SiteID)
rst <- NULL
for(i in 1: length(tags)) {
  dat2 <- dat1 %>% filter(trait == tags[i])
  aov <- Anova(lm(value~SiteID,data = dat2))
  p <- aov$`Pr(>F)`[1]
  if (p >= 0.05) p <- format(p,digit = 3) else 
    if (p >= 0.01) p <- paste(format(p,digit = 3),"*",sep = " ") else 
      if (p >= 0.001) p <- paste(format(p,digit = 3),"**",sep = " ") else 
        p <- paste("<0.001","***",sep = " ")
  rst <- rbind(rst, data.frame(element = tags[i],p = p))
}
cal <- dat1 %>% group_by(SiteID,trait) %>%
  summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                          format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                          sep = " ± ")) %>%
  spread(SiteID, value) %>%
  rename(element = trait)
rst.tot<- inner_join(rst.tot,inner_join(cal,rst),by = c("element" = "element"))

####          ####
####    WE    ####
####          ####
dat1 <- dat %>%
  filter(group == "WE")
tags <- unique(dat1$trait)
sites <- unique(dat1$SiteID)
rst <- NULL
for(i in 1: length(tags)) {
  dat2 <- dat1 %>% filter(trait == tags[i])
  aov <- Anova(lm(value~SiteID,data = dat2))
  p <- aov$`Pr(>F)`[1]
  if (p >= 0.05) p <- format(p,digit = 3) else 
    if (p >= 0.01) p <- paste(format(p,digit = 3),"*",sep = " ") else 
      if (p >= 0.001) p <- paste(format(p,digit = 3),"**",sep = " ") else 
        p <- paste("<0.001","***",sep = " ")
  rst <- rbind(rst, data.frame(element = tags[i],p = p))
}
cal <- dat1 %>% group_by(SiteID,trait) %>%
  summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                          format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                          sep = " ± ")) %>%
  spread(SiteID, value) %>%
  rename(element = trait)
rst.tot<- inner_join(rst.tot,inner_join(cal,rst),by = c("element" = "element"))

write.csv(rst.tot,"sediment/log/varingroup.csv")

