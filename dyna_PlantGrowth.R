## clean ----
rm(list = ls())

## package loading ----
library(dplyr)
library(tidyr)
library(agricolae)
library(ggplot2)

## functions ----
datareadln <- function() {
  library(dplyr)
  library(tidyr)
    read.csv("./Data/Result_PlantTillerGrowthSummary.csv") %>%
    inner_join(read.csv("./Data/meta_PlantSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    filter(group != "NV") %>%
    return()
}



## calculation ----
#progressing
source("calMeanSe.R")

calMeanSe(dat = datareadln(), group = c("SiteID", "group", "SplMonth"), 
          trait.input = "Number", trait.output = "Density")

datareadln() %>% 
  group_by(SiteID, group, SplMonth) %>%
  summarise_(x = interp(~ mean(var), var = as.name("Height_mean"))) %>%
  group_by() %>%
  select(x)
tmp = "Height"
            

datareadln() %>% 
  group_by(SiteID, group, SplMonth) %>%
  summarise(n = n(),
            Density_avg = mean(Number) * 4,
            Density_se = sd(Number)/sqrt(n()) * 4,
            Height_avg = mean(Height_mean),
            Height_se = sd(Height_mean)/sqrt(n()),
            BsDmr_avg = mean(BsDmr_mean),
            BsDmr_se = sd(BsDmr_mean)/sqrt(n()),
            LvThk_avg = mean(LvThk_mean),
            LvThk_se = sd(LvThk_mean)/sqrt(n()),
            LvCt_avg = mean(LvCt_mean),
            LvCt_se = sd(LvCt_mean)/sqrt(n()),
            LvCtG_avg = mean(LvCtG_mean),
            LvCtG_se = sd(LvCtG_mean)/sqrt(n()),
            LvCtY_avg = mean(LvCtY_mean),
            LvCtY_se = sd(LvCtY_mean)/sqrt(n()),
            TlrWg_avg = mean(TlrWg_mean),
            TlrWg_se = sd(TlrWg_mean)/sqrt(n()),
            LvWgG_avg = mean(LvWgG_mean),
            LvWgG_se = sd(LvWgG_mean)/sqrt(n()),
            StWg_avg = mean(StWg_mean),
            StWg_se = sd(StWg_mean)/sqrt(n()),
            Seed_rate_avg = mean(Seed_rate),
            Seed_rate_se = sd(Seed_rate)/sqrt(n())
            ) %>% 
  write.csv("./PlantGrowth/summary_PlantGrowthSummary.csv")

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





