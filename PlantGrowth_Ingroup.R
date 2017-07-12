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
  read.csv("./Data/Result_PlantTillerGrowthSummary.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_PlantSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    dplyr::filter(group != "NV") %>% 
    dplyr::select(SiteID,SplMonth,group,QudNo,
                  Height = Height_mean,
                  BasalDiameter = BsDmr_mean,
                  LeafThickness = LvThk_mean,
                  AboveGroundBiomass = TlrWg_mean,
                  LeafBiomass = LvWgG_mean,
                  SeedRate = Seed_rate,
                  dist) 
}

## ---------------------------------------------------
## comparing the difference inside groups 
##
## ---------------------------------------------------

dat <- datareadln() %>%
  gather(trait, value, Height:SeedRate)
dat[dat$trait == "SeedRate" & dat$value == 0,]$value <- NA

####               ####
####  Gather Plot  ####
####               ####
ggplot(aes(x = SiteID,y = value), data = dat) +
  geom_point(aes(col = group, shape = SplMonth), size = 1) + 
  ylab("") + xlab("Site") +
  facet_wrap(~trait, scales = "free") + 
  theme_bw() 
ggsave("growth/plot/gather_distribution.png",width = 10, height = 8, dpi = 600)

ggplot(aes(x = dist,y = value), data = dat) +
  geom_point(aes(col = group, shape = SplMonth), size = 2) + 
  geom_smooth(aes(x = dist,y = value), 
              data = dat %>% filter(group != "EA"), method = "lm", col = "black") +
  ylab("") + xlab("Distance from the edge of community (m)") +
  facet_wrap(~ trait + SplMonth, scales = "free",ncol = 6) + 
  theme_bw() + theme(legend.position = "bottom")
ggsave("growth/plot/gather_dist.png",width = 10, height = 8, dpi = 600)

####          ####
####    CL    ####
####          ####
dat1 <- dat %>% 
  filter(!is.na(value)) %>%
  filter(group == "CL")
tags <- unique(dat1$trait)
sites <- unique(dat1$SiteID)
months <- unique(dat1$SplMonth)
rst <- NULL
for(i in 1: length(tags)) { 
  for(j in 1:length(months)) {
    dat2 <- dat1 %>% filter(trait == tags[i],SplMonth == months[j])
    if (length(dat2$SiteID) > 3) {
      aov <- Anova(lm(value~SiteID,data = dat2))
      p <- aov$`Pr(>F)`[1]
      if (p >= 0.05) p <- format(p,digit = 3) else 
        if (p >= 0.01) p <- "*" else 
          if (p >= 0.001) p <- "**" else 
            p <- "***"
      rst <- rbind(rst, data.frame(trait = tags[i],month = months[j],p = p))
    }
  }
}
cal <- dat1 %>% group_by(SiteID,SplMonth,trait) %>%
  summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                          format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                          sep = " ± ")) %>%
  spread(SiteID, value) %>%
  rename(month = SplMonth)
rst.tot<- inner_join(cal,rst)

####          ####
####    EA    ####
####          ####
dat1 <- dat %>% 
  filter(!is.na(value)) %>%
  filter(group == "EA")
tags <- unique(dat1$trait)
sites <- unique(dat1$SiteID)
months <- unique(dat1$SplMonth)
rst <- NULL
for(i in 1: length(tags)) { 
  for(j in 1:length(months)) {
    dat2 <- dat1 %>% filter(trait == tags[i],SplMonth == months[j])
    if (length(dat2$SiteID) > 3) {
      aov <- Anova(lm(value~SiteID,data = dat2))
      p <- aov$`Pr(>F)`[1]
      if (p >= 0.05) p <- format(p,digit = 3) else 
        if (p >= 0.01) p <- "*" else 
          if (p >= 0.001) p <- "**" else 
            p <- "***"
      rst <- rbind(rst, data.frame(trait = tags[i],month = months[j],p = p))
    }
  }
}
cal <- dat1 %>% group_by(SiteID,SplMonth,trait) %>%
  summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                          format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                          sep = " ± ")) %>%
  spread(SiteID, value) %>%
  rename(month = SplMonth)
rst.tot<- inner_join(rst.tot,inner_join(cal,rst),by = c("trait" = "trait","month" = "month"))

####          ####
####    WE    ####
####          ####
dat1 <- dat %>% 
  filter(!is.na(value)) %>%
  filter(group == "WE")
tags <- unique(dat1$trait)
sites <- unique(dat1$SiteID)
months <- unique(dat1$SplMonth)
rst <- NULL
for(i in 1: length(tags)) { 
  for(j in 1:length(months)) {
    dat2 <- dat1 %>% filter(trait == tags[i],SplMonth == months[j])
    if (length(dat2$SiteID) > 3) {
      aov <- Anova(lm(value~SiteID,data = dat2))
      p <- aov$`Pr(>F)`[1]
      if (p >= 0.05) p <- format(p,digit = 3) else 
        if (p >= 0.01) p <- "*" else 
          if (p >= 0.001) p <- "**" else 
            p <- "***"
      rst <- rbind(rst, data.frame(trait = tags[i],month = months[j],p = p))
    }
  }
}
cal <- dat1 %>% group_by(SiteID,SplMonth,trait) %>%
  summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                          format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                          sep = " ± ")) %>%
  spread(SiteID, value) %>%
  rename(month = SplMonth)
rst.tot<- inner_join(rst.tot,inner_join(cal,rst),by = c("trait" = "trait","month" = "month"))

write.csv(rst.tot,"growth/log/varingroup.csv")