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
  read.csv("./Data/Result_Sediment.csv") %>%
    inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    return()
}

## calculation ----
datareadln() %>% 
  group_by(SiteID, group, SplMonth) %>%
  summarise(n = n(),
            N_avg = mean(N),N_se = sd(N)/sqrt(n()),
            C_avg = mean(C),C_se = sd(C)/sqrt(n()),
            S_avg = mean(S),S_se = sd(S)/sqrt(n()),
            P_avg = mean(P),P_se = sd(P)/sqrt(n()),
            orgC_avg = mean(orgC),orgC_se = sd(orgC)/sqrt(n()),
            Al_avg = mean(Al),Al_se = sd(Al)/sqrt(n()),
            Cr_avg = mean(Cr),Cr_se = sd(Cr)/sqrt(n()),
            Mn_avg = mean(Mn),Mn_se = sd(Mn)/sqrt(n()),
            Fe_avg = mean(Fe),Fe_se = sd(Fe)/sqrt(n()),
            Ni_avg = mean(Ni),Ni_se = sd(Ni)/sqrt(n()),
            Cu_avg = mean(Cu),Cu_se = sd(Cu)/sqrt(n()),
            Zn_avg = mean(Zn),Zn_se = sd(Zn)/sqrt(n()),
            As_avg = mean(As),As_se = sd(As)/sqrt(n()),
            Cd_avg = mean(Cd),Cd_se = sd(Cd)/sqrt(n()),
            Pb_avg = mean(Pb),Pb_se = sd(Pb)/sqrt(n())
  ) %>%
  write.csv("./SedimentElement/summary_SedimentElement.csv")

## One-way anova:SiteID ----
tmp.fit <- aov(N ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(C ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(S ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(P ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(orgC ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Al ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Cr ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Mn ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Fe ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Ni ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Cu ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Zn ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(As ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Cd ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Pb ~SiteID, data = datareadln())
HSD.test(tmp.fit, "SiteID", group=TRUE, console=TRUE, alpha = 0.05)

## One-way anova:SplMonth ----
tmp.fit <- aov(N ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(C ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(S ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(P ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(orgC ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Al ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Cr ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Mn ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Fe ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Ni ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Cu ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Zn ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(As ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Cd ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

tmp.fit <- aov(Pb ~SplMonth, data = datareadln())
HSD.test(tmp.fit, "SplMonth", group=TRUE, console=TRUE, alpha = 0.05)

## Two-way anova ----
tmp.model <- lm(N ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(C ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(S ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(P ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(orgC ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Al ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Cr ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Mn ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Fe ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Ni ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Cu ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Zn ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(As ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Cd ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

tmp.model <- lm(Pb ~ SiteID + SplMonth + SiteID * SplMonth, data = datareadln())
anova(tmp.model)
HSD.test(tmp.model, c("SiteID", "SplMonth"), group=TRUE, console=TRUE, alpha = 0.05)

## barplot ----
ggplot(data =  datareadln() %>%
         select(group,SplMonth,SiteID,N:Pb) %>%
         gather(Element,Content,N:Pb) %>%
         group_by(group,SplMonth,SiteID,Element) %>%
         summarise(n = n(),
                   mean = mean(Content,na.rm = T),
                   se = sd(Content,na.rm = T)/sqrt(n()-sum(is.na(Content)))) %>%
         filter(Element == "C" | Element == "orgC" | Element == "N" | Element == "S" | Element == "P" ))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Element~SplMonth, scales = "free") + 
  scale_y_continuous("Element Content (mg/kg)") +
  scale_x_discrete("SiteID") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404", "#013ADF", "grey50"))+
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
ggsave("./SedimentElement/bar_SedimentElement_CNSP.png",width = 6, height = 5,scale = 1)

ggplot(data = datareadln() %>%
         select(group,SplMonth,SiteID,N:Pb) %>%
         gather(Element,Content,N:Pb) %>%
         group_by(group,SplMonth,SiteID,Element) %>%
         summarise(n = n(),
                   mean = mean(Content,na.rm = T),
                   se = sd(Content,na.rm = T)/sqrt(n()-sum(is.na(Content)))) %>%
         filter(Element == "Cr" | Element == "Ni" | Element == "As" | Element == "Pb" | Element == "Cd" ))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Element~SplMonth, scales = "free") + 
  scale_y_continuous("Element Content (mg/kg)") +
  scale_x_discrete("SiteID") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404", "#013ADF", "grey50"))+
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
ggsave("./SedimentElement/bar_SedimentElement_CrNiAsPbCd.png",width = 6, height = 5,scale = 1)

ggplot(data = datareadln() %>%
         select(group,SplMonth,SiteID,N:Pb) %>%
         gather(Element,Content,N:Pb) %>%
         group_by(group,SplMonth,SiteID,Element) %>%
         summarise(n = n(),
                   mean = mean(Content,na.rm = T),
                   se = sd(Content,na.rm = T)/sqrt(n()-sum(is.na(Content)))) %>%
         filter(Element == "Al" | Element == "Fe" | Element == "Mn" | Element == "Cu" | Element == "Zn" ))+
  geom_bar(aes(x = SiteID, y = mean, fill = group),
           position = "dodge",stat = "identity", alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(x = SiteID, ymax = mean+se, ymin = mean-se, group = SplMonth),
                stat = "identity", col = "black", 
                position = position_dodge(width=0.9),
                width = 0.25, size = 0.3) + 
  facet_grid(Element~SplMonth, scales = "free") + 
  scale_y_continuous("Element Content (mg/kg)") +
  scale_x_discrete("SiteID") +
  scale_fill_manual("Location",
                    values = c("#B45F04", "#31B404", "#013ADF", "grey50"))+
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
ggsave("./SedimentElement/bar_SedimentElement_AlFeMnCuZn.png",width = 6, height = 5,scale = 1)




