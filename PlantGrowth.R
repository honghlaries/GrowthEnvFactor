## clean ----
rm(list = ls())
library(dplyr);library(tidyr)

## Functions ----
datareadln <- function() {
  read.csv("./Data/Result_PlantTillerGrowthSummary.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_PlantSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    dplyr::filter(group != "NV") %>% 
    dplyr::select(SiteID,SplMonth,group,Number,
                  Height = Height_mean,
                  BsDmr = BsDmr_mean,
                  LvThk = LvThk_mean,
                  LvCt = LvCt_mean,
                  LvCtG = LvCtG_mean,
                  LvCtY = LvCtY_mean,
                  TlrWg = TlrWg_mean,
                  LvWgG = LvWgG_mean,
                  StWg = StWg_mean,
                  SeedRate = Seed_rate) 
}

meanseCal <- function(dat) { ## calulate and output mean and se
  se <- function(v) {
    if(length(v) < 2) {NA} else {sd(v, na.rm = T)/sqrt(length(v)-sum(is.na(v)))}
  }
  dat %>%
    dplyr::group_by(SiteID, group, SplMonth) %>%
    dplyr::summarise(Density_avg = mean(Number,na.rm = T) * 4,Density_se = se(Number) * 4,
                     Height_avg = mean(Height,na.rm = T),Height_se = se(Height),
                     BsDmr_avg = mean(BsDmr,na.rm = T),BsDmr_se = se(BsDmr),
                     LvThk_avg = mean(LvThk,na.rm = T),LvThk_se = se(LvThk),
                     LvCt_avg = mean(LvCt,na.rm = T),LvCt_se = se(LvCt),
                     LvCtG_avg = mean(LvCtG,na.rm = T),LvCtG_se = se(LvCtG),
                     LvCtY_avg = mean(LvCtY,na.rm = T),LvCtY_se = se(LvCtY),
                     TlrWg_avg = mean(TlrWg,na.rm = T),TlrWg_se = se(TlrWg),
                     LvWg_avg = mean(LvWgG,na.rm = T),LvWg_se = se(LvWgG),
                     StWg_avg = mean(StWg,na.rm = T),StWg_se = se(StWg),
                     SeedRate_avg = mean(SeedRate,na.rm = T),SeedRate_se = se(SeedRate)) %>%
    write.csv("growth/log/GrowthTraits.csv")
}


## Basic Stat information ----
meanseCal(datareadln())



## Two-way anova ----
tmp.model <- lm(Number ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Height_mean ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(BsDmr_mean ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(LvThk_mean ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(LvCtG_mean ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(LvCtY_mean ~ SiteID, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(LvCt_mean ~ SiteID, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(TlrWg_mean ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(LvWgG_mean ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(StWg_mean ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Seed_rate ~ SiteID, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID"), group=TRUE, console=TRUE, alpha = 0.05)

## barplot ----
ggplot(data =  datareadln() %>%
         mutate(Density = 4* Number) %>%
         select(group,SplMonth,SiteID,Density,contains("_mean"),Seed_rate) %>%
         gather(Traits,Value,Density:Seed_rate) %>%
         group_by(group,SplMonth,SiteID,Traits) %>%
         summarise(n = n(),
                   mean = mean(Value,na.rm = T),
                   se = sd(Value,na.rm = T)/sqrt(n()-sum(is.na(Value)))) %>%
         filter(SplMonth == "Nov") %>%
         filter(Traits == "Seed_rate"| Traits == "LvCtY_mean" ))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Traits~SplMonth, scales = "free") + 
  #scale_x_continuous("") +
  scale_y_continuous("") +
  scale_x_discrete("SiteID") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404", "#013ADF"))+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 7, angle = 30),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
ggsave("./PlantGrowth/bar_PlantGrowth_SexReproduce.png",width = 5, height = 4,scale = 1)

ggplot(data =  datareadln() %>%
         mutate(Density = 4* Number) %>%
         select(group,SplMonth,SiteID,Density,contains("_mean"),Seed_rate) %>%
         gather(Traits,Value,Density:Seed_rate) %>%
         group_by(group,SplMonth,SiteID,Traits) %>%
         summarise(n = n(),
                   mean = mean(Value,na.rm = T),
                   se = sd(Value,na.rm = T)/sqrt(n()-sum(is.na(Value)))) %>%
         filter(Traits == "Density"))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Traits~SplMonth, scales = "free") + 
  #scale_x_continuous("") +
  scale_y_continuous("") +
  scale_x_discrete("SiteID") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404", "#013ADF"))+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 7, angle = 30),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
ggsave("./PlantGrowth/bar_PlantGrowth_Density.png",width = 5, height = 2.5,scale = 1)

ggplot(data =  datareadln() %>%
         mutate(Density = 4* Number) %>%
         select(group,SplMonth,SiteID,Density,contains("_mean"),Seed_rate) %>%
         gather(Traits,Value,Density:Seed_rate) %>%
         group_by(group,SplMonth,SiteID,Traits) %>%
         summarise(n = n(),
                   mean = mean(Value,na.rm = T),
                   se = sd(Value,na.rm = T)/sqrt(n()-sum(is.na(Value)))) %>%
         filter(Traits == "Height_mean" | Traits == "BsDmr_mean" | Traits == "LvThk_mean" ))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Traits~SplMonth, scales = "free") + 
  #scale_x_continuous("") +
  scale_y_continuous("") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404", "#013ADF"))+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 7, angle = 30),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
ggsave("./PlantGrowth/bar_PlantGrowth_Structure.png",width = 6, height = 5,scale = 1)

ggplot(data =  datareadln() %>%
         mutate(Density = 4* Number) %>%
         select(group,SplMonth,SiteID,Density,contains("_mean"),Seed_rate) %>%
         gather(Traits,Value,Density:Seed_rate) %>%
         group_by(group,SplMonth,SiteID,Traits) %>%
         summarise(n = n(),
                   mean = mean(Value,na.rm = T),
                   se = sd(Value,na.rm = T)/sqrt(n()-sum(is.na(Value)))) %>%
         filter(Traits == "LvCt_mean" | Traits == "LvCtG_mean" | Traits == "LvCtY_mean" ))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Traits~SplMonth, scales = "free") + 
  #scale_x_continuous("") +
  scale_y_continuous("") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404",  "#013ADF"))+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 7, angle = 30),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
ggsave("./PlantGrowth/bar_PlantGrowth_LeafNum.png",width = 6, height = 5,scale = 1)

ggplot(data =  datareadln() %>%
         mutate(Density = 4* Number) %>%
         select(group,SplMonth,SiteID,Density,contains("_mean"),Seed_rate) %>%
         gather(Traits,Value,Density:Seed_rate) %>%
         group_by(group,SplMonth,SiteID,Traits) %>%
         summarise(n = n(),
                   mean = mean(Value,na.rm = T),
                   se = sd(Value,na.rm = T)/sqrt(n()-sum(is.na(Value)))) %>%
         filter(Traits == "TlrWg_mean" | Traits == "LvWgG_mean" | Traits == "StWg_mean" ))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Traits~SplMonth, scales = "free") + 
  #scale_x_continuous("") +
  scale_y_continuous("") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404",  "#013ADF"))+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 7, angle = 30),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
ggsave("./PlantGrowth/bar_PlantGrowth_Biomass.png",width = 6, height = 5,scale = 1)





