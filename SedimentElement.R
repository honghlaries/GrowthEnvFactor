## clean ----
rm(list = ls())
library(dplyr);library(tidyr)

## Functions ----
datareadln <- function() {
  read.csv("./Data/Result_Sediment.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    dplyr::filter(group != "NV") %>% 
    dplyr::select(N:Pb,SiteID,SplMonth,group) %>%
    return()
}

meanseCal <- function(dat) { ## calulate and output mean and se
  se <- function(v) {
    if(length(v) < 2) {NA} else {sd(v, na.rm = T)/sqrt(length(v)-sum(is.na(v)))}
  }
  dat %>%
    dplyr::group_by(SiteID, group, SplMonth) %>%
    dplyr::summarise(N_avg = mean(N,na.rm = T),N_se = se(N),
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

shapiroTest <- function(dat, tag, grouplv, SplMonthlv) {
  dat <- dat %>% filter(elem == tag) %>% filter(!is.na(resp))
  shapiro <- NULL
  for(i in 1:length(grouplv)) {
    for(j in 1:length(SplMonthlv)) {
      subdat <- dat %>% filter(group == grouplv[i]) %>% filter(SplMonth == SplMonthlv[j])
      result <- stats::shapiro.test(subdat$resp)
      shapiro <- rbind(shapiro, data.frame(tag = tag,
                                      group = grouplv[i],
                                      SplMonth = SplMonthlv[j],
                                      W = result$statistic,
                                      Pvalue = result$p.value))
    }
  }
  shapiro
}

modelCompare <- function(dat, tag, fact, log = TRUE) {
  library(lme4)
  dat <- dat %>% filter(elem == tag) %>% filter(!is.na(resp))
  formulaList <- NULL; modlog <- NULL
  for (j in 1:length(fact)) {
    formulaList <- c(formulaList, as.formula(paste("resp","~",fact[j])))
  }
  mod.champion <- lmer(formulaList[[1]], data = dat)
  for (i in 2:length(formulaList)) {
    mod.challenger <- lmer(formulaList[[i]], data = dat)
    comp <- anova(mod.champion, mod.challenger)
    mod <- as.character(formula(mod.champion@call))[3]
    modp <- as.character(formula(mod.challenger@call))[3]
    modlog <- rbind(modlog, cbind(comp,
                                  model = c(paste(tag,"~",mod),paste(tag,"~",modp)),
                                  ref = c("champion", "challenger"),
                                  result = if (comp$AIC[2] < comp$AIC[1]) c("", "Win") else c("Win", ""),
                                  time = date())) 
    if (comp$AIC[2] < comp$AIC[1]) mod.champion <- mod.challenger
  }
  if (!file.exists("sediment/log/SedimentElementModel.csv")) {
    write.table(t(c("Df","AIC","BIC","Loglik","deviance","Chisq","ChiDf","Pr(>Chisq)","model","ref","result","time")),
                file = "sediment/log/SedimentElementModel.csv", sep = ",", row.names = F, col.names = F)
  } 
  if(log) {
    write.table(modlog, append = T, sep = ",", row.names = F, col.names = F,
                file = "sediment/log/SedimentElementModel.csv")
  }
  mod.champion
}

plotFEsim2 <- function (fe, modinf = NULL, level = 0.95, stat = "mean", sd = TRUE,  
                        sigmaScale = NULL, oddsRatio = FALSE, glv, gcode) {
  ## plotFEsim2, modified from merTools::plotFEsim
  if (!missing(sigmaScale)) {
    fe[, "sd"] <- fe[, "sd"]/sigmaScale
    fe[, stat] <- fe[, stat]/sigmaScale
  }
  fe[, "sd"] <- fe[, "sd"] * qnorm(1 - ((1 - level)/2))
  fe[, "ymax"] <- fe[, stat] + fe[, "sd"]
  fe[, "ymin"] <- fe[, stat] - fe[, "sd"]
  hlineInt <- 0
  if (oddsRatio == TRUE) {
    fe[, "ymax"] <- exp(fe[, "ymax"])
    fe[, stat] <- exp(fe[, stat])
    fe[, "ymin"] <- exp(fe[, "ymin"])
    hlineInt <- 1
  }
  xvar <- "term"
  fe$term <- as.character(fe$term)
  fe$term <- factor(fe$term, levels = fe[order(fe[, stat]), 1])
  
  gcol <- rep("black",length(levels(fe$term)))
  gsize <- rep(1,length(levels(fe$term)))
  gltp <- rep(2,length(levels(fe$term)))
  for (m in 1: length(glv)) {
    gcol[levels(fe$term) == glv[m]] <- gcode[m]
    gsize[levels(fe$term) == glv[m]] <- 4
  }
  llv <- c("Nov","Jul","EA","NV","WE")
  for (m in 1: length(llv)) {
    gltp[levels(fe$term) == llv[m]] <- 1
  }
  p <- ggplot(aes_string(x = xvar, y = stat, ymax = "ymax", ymin = "ymin"), data = fe) + 
    geom_hline(yintercept = hlineInt, color = I("red")) +
    geom_point(aes(col = term, size = term)) + 
    labs(x = "Group", y = "Fixed Effect") + 
    scale_color_manual("Group", breaks = levels(fe$term), values = gcol) + 
    scale_size_manual("Group", breaks = levels(fe$term), values = gsize) + 
    coord_flip() + 
    theme_HL
  if (sd) {
    p <- p + geom_errorbar(aes(linetype = term), width = 0.2) +
      scale_linetype_manual("Group", breaks = levels(fe$term), values = gltp)
  }
  p
}

plotREsim2facet <- function (data, level = 0.95, stat = "median", sd = TRUE, sigmaScale = NULL, 
                        oddsRatio = FALSE, labs = FALSE, taglv = NA, ncol) {
  ## plotREsim2, modified from merTools::plotREsim
  if (!missing(sigmaScale)) {
    data[, "sd"] <- data[, "sd"]/sigmaScale
    data[, stat] <- data[, stat]/sigmaScale
  }
  data[, "sd"] <- data[, "sd"] * qnorm(1 - ((1 - level)/2))
  data[, "ymax"] <- data[, stat] + data[, "sd"]
  data[, "ymin"] <- data[, stat] - data[, "sd"]
  data[, "sig"] <- data[, "ymin"] > 0 | data[, "ymax"] < 0
  hlineInt <- 0
  if (oddsRatio == TRUE) {
    data[, "ymax"] <- exp(data[, "ymax"])
    data[, stat] <- exp(data[, stat])
    data[, "ymin"] <- exp(data[, "ymin"])
    hlineInt <- 1
  }
  #data <- data[order(data[, "groupFctr"], data[, "term"], data[,stat]), ]
  #rownames(data) <- 1:nrow(data)
  data[, "xvar"] <- factor(paste(data$groupFctr, data$groupID, 
                                 sep = ""), levels = unique(paste(data$groupFctr, data$groupID, 
                                                                  sep = "")), ordered = TRUE)
  if (labs == TRUE) {
    xlabs.tmp <- element_text(face = "bold", angle = 90, 
                              vjust = 0.5)
  }
  else {
    data[, "xvar"] <- as.numeric(data[, "xvar"])
    xlabs.tmp <- element_blank()
  }
  if (!is.na(taglv)) { 
    tmp <- NULL
    for (k in 1: length(taglv)) {
      tmp <- rbind(tmp, dplyr::filter(data, tag == taglv[k]))
    }
    data <- tmp
    data$tag <- factor(data$tag, levels = taglv)
  }
  p <- ggplot(data, aes_string(x = "xvar", y = stat, ymax = "ymax", 
                               ymin = "ymin")) +
    geom_hline(yintercept = hlineInt, color = I("red"), size = I(1.1)) +
    geom_point(color = "gray75", alpha = 1/(nrow(data)^0.33),  size = I(0.5)) +
    geom_point(data = subset(data, sig ==  TRUE), size = I(3)) + 
    labs(x = "Group", y = "Effect Range") + 
    theme_bw() #+ theme(panel.grid.major = element_blank(), 
  #        panel.grid.minor = element_blank(), axis.text.x = xlabs.tmp, 
  #        axis.ticks.x = element_blank())
  if (sd) {
    p <- p + geom_pointrange(alpha = 1/(nrow(data)^0.33)) + 
      geom_pointrange(data = subset(data, sig == TRUE), 
                      alpha = 0.25)
  }
  p + facet_wrap(~tag + groupFctr, ncol = ncol, scales = "free")
}

plotFEsim2facet <- function (fe, modinf = NULL, level = 0.95, stat = "mean", sd = TRUE,  
                        sigmaScale = NULL, oddsRatio = FALSE, taglv = NA, glv, gcode, ncol) {
  ## plotFEsim2, modified from merTools::plotFEsim
  if (!missing(sigmaScale)) {
    fe[, "sd"] <- fe[, "sd"]/sigmaScale
    fe[, stat] <- fe[, stat]/sigmaScale
  }
  fe[, "sd"] <- fe[, "sd"] * qnorm(1 - ((1 - level)/2))
  fe[, "ymax"] <- fe[, stat] + fe[, "sd"]
  fe[, "ymin"] <- fe[, stat] - fe[, "sd"]
  hlineInt <- 0
  if (oddsRatio == TRUE) {
    fe[, "ymax"] <- exp(fe[, "ymax"])
    fe[, stat] <- exp(fe[, stat])
    fe[, "ymin"] <- exp(fe[, "ymin"])
    hlineInt <- 1
  }
  xvar <- "term"
  fe$term <- as.character(fe$term)
  fe$term <- factor(fe$term, levels = fe[order(fe[, stat]), 1])
  if (!is.na(taglv)) { 
    tmp <- NULL
    for (k in 1: length(taglv)) {
      tmp <- rbind(tmp, dplyr::filter(fe, tag == taglv[k]))
    }
    fe <- tmp
    fe$tag <- factor(fe$tag, levels = taglv)
  }
  gcol <- rep("black",length(levels(fe$term)))
  gsize <- rep(1,length(levels(fe$term)))
  gltp <- rep(2,length(levels(fe$term)))
  for (m in 1: length(glv)) {
    gcol[levels(fe$term) == glv[m]] <- gcode[m]
    gsize[levels(fe$term) == glv[m]] <- 4
  }
  llv <- c("Nov","Jul","EA","NV","WE")
  for (m in 1: length(llv)) {
    gltp[levels(fe$term) == llv[m]] <- 1
  }
  p <- ggplot(aes_string(x = xvar, y = stat, ymax = "ymax", ymin = "ymin"), data = fe) + 
    geom_hline(yintercept = hlineInt, color = I("red")) +
    geom_point(aes(col = term, size = term)) + 
    labs(x = "Group", y = "Fixed Effect") + 
    scale_color_manual("Group", breaks = levels(fe$term), values = gcol) + 
    scale_size_manual("Group", breaks = levels(fe$term), values = gsize) + 
    facet_wrap(~ tag, ncol = ncol, scales = "free") +
    # coord_flip() + 
    theme_HL
  if (sd) {
    p <- p + geom_errorbar(aes(linetype = term), width = 0.2) +
      scale_linetype_manual("Group", breaks = levels(fe$term), values = gltp)
  }
  if(!missing(modinf)) {
    p <- p + geom_text(x = 2, y = 5, label = "abc")
  }
  p
}

multiElementMod <- function(dat, 
                            fact = c("group + (1|SiteID)", 
                                           "group + (1|SiteID:SplMonth)",
                                           "group + (1|SiteID) + (1|SiteID:SplMonth)", 
                                           "group + SplMonth + (1|SiteID)",
                                           "group + SplMonth + (1|SiteID:SplMonth)",
                                           "group + SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                                           "group * SplMonth + (1|SiteID)",
                                           "group * SplMonth + (1|SiteID:SplMonth)",
                                           "group * SplMonth + (1|SiteID) + (1|SiteID:SplMonth)"),
                            SplMonthlv = c("Apr","Jul","Nov"), grouplv = c("EA","CL","WE"), 
                            glv = c("EA","WE"), gcode = c("#31B404","#013ADF"),
                            tag = c("N","C","orgC","S","P","Al","Fe","Mn","Cu","Zn","Ni","Cr","Pb","As","Cd"),
                            theme_HL = theme_bw() + 
                              theme(legend.position = "none",
                                    axis.text = element_text(angle = 30)),
                            archiplot = TRUE, log = TRUE) {
  library(merTools);library(lsmeans);library(multcompView);library(car)
  fe.g <- NULL; re.g <- NULL; modinf.g <- NULL; shapiro.g <- NULL; posthoc.g <- NULL; modavo.g <- NULL; 
  resid.g <- NULL;
  for (i in 1:length(tag)) {
    
    #shapiro
    shapiroTest(dat = dat, tag = tag[i], grouplv = grouplv, SplMonthlv = SplMonthlv) -> shapiro
    print(shapiro); 
    shapiro.g <- rbind(shapiro.g, cbind(shapiro, tag[i]))
    
    #mod choose
    mod <- modelCompare(dat = dat, tag = tag[i], fact = fact,log = log)
    paste(tag[i],"~",as.character(formula(mod@call))[3]) -> modinf
    print(modinf); 
    modinf.g <- rbind(modinf.g, cbind(modinf, tag[i]))
    
    #anova
    nyma <- anova(mod)
    nymb <- Anova(mod)
    cbind(nyma[,c(1,2,3)],Chisq = nymb[,1], Fvalue = nyma[,4],Pr = nymb[,3], 
          FixedFact = row.names(nyma), tag = tag[i]) -> modavo; print(modavo); 
    modavo.g <- rbind(modavo.g, cbind(modavo, tag[i]))
    
    #posthoc
    lsmeans::cld(lsmeans(object = mod, adjust = "tukey",
                         specs = pairwise ~ group + SplMonth),
                 Letters = LETTERS) -> posthoc; print(posthoc); 
    posthoc.g <- rbind(posthoc.g,cbind(posthoc, tag[i]))
    
    # fe
    fe <- FEsim(mod)
    fe.g <- rbind(fe.g, cbind(tag = tag[i],fe))
    ggsave(plot = plotFEsim2(fe%>% filter(term != "(Intercept)") %>%
                               mutate(term = gsub("group","",term)) %>% 
                               mutate(term = gsub("SplMonth","",term)),
                             glv = glv, gcode = gcode), 
           paste("sediment/plot/",tag[i],"_Fixeff.png",sep=""),
           width = 6, height = 4)
    #re
    re <- REsim(mod)
    re.g <- rbind(re.g, cbind(tag = tag[i],re))
    ggsave(plot = merTools::plotREsim(re), 
           paste("sediment/plot/",tag[i],"_Raneff.png",sep=""),
           width = 6, height = 4)
    #resid
    png(paste("sediment/plot/",tag[i],"_resid.png",sep=""),
        width = 30, height = 20, units = "cm", res = 600)
    plot(mod)
    dev.off()
    stats::shapiro.test(resid(mod)) -> resid
    resid.g <- rbind(resid.g, 
                     data.frame(tag = tag[i],
                                W = resid$statistic,
                                Pvalue = resid$p.value))
  }
  ## output
  if (log) {write.csv(fe.g, "sediment/log/FixedEff.csv")
    write.csv(re.g, "sediment/log/RandomEff.csv")
    write.csv(modinf.g, "sediment/log/ModelChoice.csv") 
    write.csv(shapiro.g, "sediment/log/ShapiroRawData.csv") 
    write.csv(posthoc.g, "sediment/log/Posthoc.csv") 
    write.csv(modavo.g, "sediment/log/ModelAnova.csv") 
    write.csv(resid.g, "sediment/log/ShapiroResid.csv")
  }
  "DONE!"
}

## Basic Stat information ----
meanseCal(datareadln())

## LMM fiting and ploting ----
multiElementMod(dat = datareadln() %>% gather(key = elem, value = resp, N:Pb))


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




