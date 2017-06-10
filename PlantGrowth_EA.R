## ---------------------------------------------------
## init
##
## ---------------------------------------------------
rm(list = ls())
library(dplyr);library(tidyr);library(ggplot2)

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
                  Density = Number,
                  Height = Height_mean,
                  BasalDiameter = BsDmr_mean,
                  LeafThickness = LvThk_mean,
                  LeafCount = LvCt_mean,
                  GreenLeafCount = LvCtG_mean,
                  FallenLeafCount = LvCtY_mean,
                  AboveGroundBiomass = TlrWg_mean,
                  LeafBiomass = LvWgG_mean,
                  StemBiomass = StWg_mean,
                  SeedRate = Seed_rate) 
}

singlePlot <- function(dat,tag) {
  dat <- dat %>% group_by() %>%  filter(trait == tag) %>% filter(!is.na(mean)) %>%
    mutate(SplMonth = factor(as.character(SplMonth)))
  ggplot(data = dat) + 
    geom_bar(aes(x = SplMonth, y = mean, fill = SiteID),stat = "identity", position = position_dodge(0.9)) +
    geom_errorbar(aes(x = SplMonth, ymin = mean - se, ymax = mean + se, fill = SiteID),
                  position = position_dodge(0.9), width = 0.3) +
    xlab("") + ylab(tag) +
    scale_fill_discrete("Location",breaks = c("EA1","EA2"),labels = c("Seaward(EA1)","Landward(EA2)")) 
}

seasonal.t.test <- function(dat,tag) {
  library(dplyr)
  dat <- dat %>% filter(trait == tag) %>% filter(!is.na(value)) %>%
    mutate(SplMonth = factor(as.character(SplMonth)))
  t.sum <- NULL
  for (i in 1: length(unique(dat$SplMonth))) {
    subdat <- dat %>% filter(SplMonth == unique(dat$SplMonth)[i])
    t <- t.test(value~SiteID, data = subdat, var.equal = TRUE)
    conf <- format(t$conf.int, digit = 3)
    p <- format(t$p.value,digit = 3)
    t.sum <- rbind(t.sum, data.frame(
                   tag = tag,
                   month = unique(dat$SplMonth)[i],
                   conf = paste("[", conf[1], ",", conf[2], "]", sep = ""),
                   p = p))
  }
  t.sum
}

doTandPlot <- function(dat,tag) {
  library(dplyr);library(ggplot2)
  singlePlot(dat%>%
               dplyr::group_by(SiteID,SplMonth,trait) %>%
               dplyr::summarise(mean = mean(value,na.rm = T),se = (sd(value,na.rm = T)/sum(!is.na(value)))),
             tag) 
  ggsave(filename = paste("growth/plot/bar_",tag,"_EA.png",sep = ""), width = 6, height = 4)
  seasonal.t.test(dat,tag)
}

shapiroTest <- function(dat) {
  print(stats::shapiro.test(dat$value))
  par(mar = c(4, 4, 0.5, 0.5), mfrow = c(1, 2))
  qqnorm(dat$value)
  plot(density(dat$value))
  print("pass? 1: Yes; 0: No"); if(!scan(n = 1)) return("test of normality failed")
  
  SiteIDlv <- unique(dat$SiteID)
  for(i in 1:length(SiteIDlv)) {
    subdat <- dat %>% filter(SiteID == SiteIDlv[i])
    print("------SHAPIRO TEST RESULT-------")
    print(paste(SiteIDlv[i],":"))
    stats::shapiro.test(subdat$value) %>% print()
    print("--------------------------------")
    print("                                ")
  }
  SplMonthlv <- unique(dat$SplMonth)
  for(i in 1:length(SplMonthlv)) {
    subdat <- dat %>% filter(SplMonth == SplMonthlv[i])
    print("------SHAPIRO TEST RESULT-------")
    print(paste(SplMonthlv[i],":"))
    stats::shapiro.test(subdat$value) %>% print()
    print("--------------------------------")
    print("                                ")
  }
  par(mar = c(4, 4, 0.5, 0.5), mfrow = c(2, 3))
  plot(density(dat$value[dat$SiteID == "EA1"]), main = "EA1")
  plot(density(dat$value[dat$SiteID == "EA2"]), main = "EA2")
  plot(1)
  if("Apr" %in% unique(dat$SplMonth)) plot(density(dat$value[dat$SplMonth == "Apr"]), main = "Apr") else plot(1)
  if("Jul" %in% unique(dat$SplMonth)) plot(density(dat$value[dat$SplMonth == "Jul"]), main = "Jul") else plot(1)
  plot(density(dat$value[dat$SplMonth == "Nov"]), main = "Nov")
  par(mfrow = c(1,1))
  print("pass? 1: Yes; 0: No"); if(!scan(n = 1)) return("test of normality failed")
  return("test of normality passed!")
}

leveneTest <- function(dat){
  library(car)
  plot(value ~ as.factor(SiteID), data = dat)
  print(with(car::leveneTest(value,group = SiteID), data = dat))
  print("pass? 1: Yes; 0: No"); if(!scan(n = 1)) return("homogeneity test of variances failed")
  
  plot(value ~ as.factor(SplMonth), data = dat)
  print(with(car::leveneTest(value,group = SplMonth), data = dat))
  print("pass? 1: Yes; 0: No"); if(!scan(n = 1)) return("homogeneity test of variances failed")
  
  plot(value ~ as.factor(paste(SiteID,SplMonth,sep = ".")), data = dat)
  print(with(car::leveneTest(value,group = paste(SiteID,SplMonth,sep = ".")), data = dat))
  print("pass? 1: Yes; 0: No"); if(!scan(n = 1)) return("homogeneity test of variances failed")
  return("homogeneity test of variances passed!")
}

doadjustedTest <- function(dat, tag) {
  
  aov <- aov(value ~ SplMonth * SiteID, data = dat)
  print(summary(aov))
  scan(n = 0)
  
  library(multcomp)
  library(multcompView)
  hsd <- TukeyHSD(aov)
  #plot(hsd)
  print(multcompLetters4(aov,hsd))
  scan(n = 0)
  
  seasonal.t.test(dat,tag)
}
## ---------------------------------------------------
## comparing before/after the dyke 
##
## ---------------------------------------------------
dat <- datareadln() %>% 
  dplyr::filter(group == "EA") %>% 
  dplyr::select(SiteID, SplMonth, Density:SeedRate)  %>% 
  dplyr::mutate(L_S_Ratio = LeafBiomass / StemBiomass) %>%
  tidyr::gather(trait, value, Density:L_S_Ratio) 
  

ggplot(data = dat%>%
         dplyr::group_by(SiteID,SplMonth,trait) %>%
         dplyr::summarise(mean = mean(value,na.rm = T),se = (sd(value,na.rm = T)/sum(!is.na(value))))) + 
  geom_bar(aes(x = SplMonth, y = mean, fill = SiteID),stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(x = SplMonth, ymin = mean - se, ymax = mean + se, fill = SiteID),
                position = position_dodge(0.9), width = 0.3) +
  xlab("") + ylab("") +
  scale_fill_discrete("Location",breaks = c("EA1","EA2"),labels = c("Seaward(EA1)","Landward(EA2)")) + 
  facet_wrap(facets = ~ trait, nrow = 4, scales = "free") 
ggsave(filename = "growth/plot/bar_growth_EA.png")

t.summary <- NULL
tagList <- unique(dat$trait)
for(i in 1: length(tagList)){
  t.summary <- rbind(t.summary,doTandPlot(dat,tagList[i]))
} %>% write.csv("growth/log/t.test_EA.csv",row.names = F)

read.csv("growth/log/t.test_EA.csv") %>% print(row.names = F)
## ---------------------------------------------------
## scaling and re-analysis
##
## ---------------------------------------------------

# Density
tag <- "Density"

dat1 <- dat %>% filter(trait == tag) %>% mutate(value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Height
tag <- "Height"

dat1 <- dat %>% filter(trait == tag) %>% mutate(value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Basal diameter
tag <- "BasalDiameter"

dat1 <- dat %>% filter(trait == tag) %>% mutate(value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# LeafThickness
tag <- "LeafThickness"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
    value = log(1+log(value)))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# LeafCount
tag <- "LeafCount"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# GreenLeafCount
tag <- "GreenLeafCount"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# FallenLeafCount
tag <- "FallenLeafCount"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

aov <- aov(value ~ SiteID, data = dat1)
print(summary(aov))

seasonal.t.test(dat,tag)

# AboveGroundBiomass
tag <- "AboveGroundBiomass"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = log(value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# LeafBiomass
tag <- "LeafBiomass"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = log(log(value)))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# StemBiomass
tag <- "StemBiomass"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = log(value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# StemBiomass
tag <- "StemBiomass"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = log(value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# SeedRate
tag <- "SeedRate"

dat1 <- dat %>% filter(SplMonth == "Nov" & trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

aov <- aov(value ~ SiteID, data = dat1)
print(summary(aov))

seasonal.t.test(dat,tag)

# L_S_Ratio
tag <- "L_S_Ratio"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

