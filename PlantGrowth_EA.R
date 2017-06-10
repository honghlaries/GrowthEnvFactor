## clean 
rm(list = ls())
library(dplyr);library(tidyr);library(ggplot2)

## Functions 
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
    t <- t.test(value~SiteID, data = subdat)
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
} %>% write.csv("growth/log/t.test_EA.csv")

## ---------------------------------------------------
## scaling and re-analysis
##
## ---------------------------------------------------

# Density
dat1 <- dat %>% filter(trait == "Density") %>% mutate(value = log(value))
seasonal.t.test(dat1,"Density")

