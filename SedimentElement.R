## clean ----
rm(list = ls())
library(dplyr,tidyr)

## functions ----
datareadln <- function() {
  read.csv("./Data/Result_Sediment.csv") %>%
    inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    filter(group != "NV") %>%
    return()
}

meanseCal <- function(dat) { ## calulate and output mean and se ----
  se <- function(v) {
    if(length(v) < 2) {NA} else {sd(v, na.rm = T)/sqrt(length(v)-sum(is.na(v)))}
  }
  dat %>%
    group_by(SiteID, group, SplMonth) %>%
    summarise(N_avg = mean(N,na.rm = T),N_se = se(N),
              C_avg = mean(C,na.rm = T),C_se = se(C),
              S_avg = mean(S,na.rm = T),S_se = se(S),
              P_avg = mean(P,na.rm = T),P_se = se(P),
              orgC_avg = mean(orgC,na.rm = T),orgC_se = se(orgC),
              Al_avg = mean(Al,na.rm = T),Al_se = se(Al),
              Cr_avg = mean(Cr,na.rm = T),Cr_se = se(Cr),
              Mn_avg = mean(Mn,na.rm = T),Mn_se = se(Mn),
              Fe_avg = mean(Fe,na.rm = T),Fe_se = se(Fe),
              Ni_avg = mean(Ni,na.rm = T),Ni_se = se(Ni),
              Cu_avg = mean(Cu,na.rm = T),Cu_se = se(Cu),
              Zn_avg = mean(Zn,na.rm = T),Zn_se = se(Zn),
              As_avg = mean(As,na.rm = T),As_se = se(As),
              Cd_avg = mean(Cd,na.rm = T),Cd_se = se(Cd),
              Pb_avg = mean(Pb,na.rm = T),Pb_se = se(Pb)) %>%
    write.csv("sediment/log/SedimentElement.csv")
}

## STAT begin
meanseCal(datareadln())
library(lme4)
library(car)

modelCompare <- function(dat, formulaList, log = TRUE) {
  if (!file.exists("sediment/log/SedimentElementModel.csv")) {
    write.table(t(c("Df","AIC","BIC","Loglik","deviance","Chisq","ChiDf","Pr(>Chisq)","mod",	"ref","result")),
                file = "sediment/log/SedimentElementModel.csv", sep = ",", row.names = F, col.names = F)
  }
  mod.champion <- lmer(formulaList[[1]], data = dat)
  for (i in 2:length(formulaList)) {
    mod.challenger <- lmer(formulaList[[i]], data = dat)
    comp <- anova(mod.champion, mod.challenger)
    mod <- as.character(formula(mod.champion@call))
    modp <- as.character(formula(mod.challenger@call))
    if(log) write.table(cbind(comp,
                            mod = c(paste(mod[2],"~",mod[3]),paste(modp[2],"~",modp[3])),
                            ref = c("champion", "challenger"),
                            result = if (comp$AIC[2] < comp$AIC[1]) c("", "Win") else c("Win", "")) ,
                      file = "sediment/log/SedimentElementModel.csv", 
                      append = T, sep = ",", row.names = F, col.names = F)
    if (comp$AIC[2] < comp$AIC[1]) mod.champion <- mod.challenger
  }
  mod.champion
}

formulaList <- c(N ~ group + (1|group:SiteID) + (1|SiteID:SplMonth),
                 N ~ group + SplMonth + (1|group:SiteID) + (1|SiteID:SplMonth),
                 N ~ group * SplMonth + (1|group:SiteID) + (1|SiteID:SplMonth))
modelCompare(dat = datareadln(), formulaList = formulaList)

a <- "N";b <- "group + (1|group:SiteID) + (1|SiteID:SplMonth)"
as.formula(paste(a,"~",b)) 

mod.challenger <- lmer(as.formula(paste(a,"~",b)), data = dat)



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




