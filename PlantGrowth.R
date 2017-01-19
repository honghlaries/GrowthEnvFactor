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
    dplyr::select(SiteID,SplMonth,group,
                  Density = Number *4,
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
    dplyr::summarise(Density_avg = mean(Density,na.rm = T),Density_se = se(Density),
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

shapiroTest <- function(dat, tag, grouplv, SplMonthlv) {
  dat <- dat %>% filter(traits == tag) %>% filter(!is.na(resp))
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
  dat <- dat %>% filter(traits == tag) %>% filter(!is.na(resp))
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
  if (!file.exists("growth/log/GrowthTraitsModel.csv")) {
    write.table(t(c("Df","AIC","BIC","Loglik","deviance","Chisq","ChiDf","Pr(>Chisq)","model","ref","result","time")),
                file = "growth/log/GrowthTraitsModel.csv", sep = ",", row.names = F, col.names = F)
  } 
  if(log) {
    write.table(modlog, append = T, sep = ",", row.names = F, col.names = F,
                file = "growth/log/GrowthTraitsModel.csv")
  }
  mod.champion
}

plotFEsim2 <- function (fe, modinf = NULL, level = 0.95, stat = "mean", sd = TRUE,  
                        sigmaScale = NULL, oddsRatio = FALSE, glv, gcode, theme) {
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
    theme
  if (sd) {
    p <- p + geom_errorbar(aes(linetype = term), width = 0.2) +
      scale_linetype_manual("Group", breaks = levels(fe$term), values = gltp)
  }
  p
}

plotFEsim2facet <- function (fe, modinf = NULL, level = 0.95, stat = "mean", sd = TRUE,  
                             sigmaScale = NULL, oddsRatio = FALSE, taglv = NA, glv, gcode, ncol, theme) {
  ## plotFEsim2facet, modified from merTools::plotFEsim
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
  fe$term <- factor(fe$term, levels = c("WE:Jul","WE:Nov","WE","Jul","Nov","EA","EA:Jul","EA:Nov"))
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
  gltp <- rep(0.5,length(levels(fe$term)))
  for (m in 1: length(glv)) {
    gcol[levels(fe$term) == glv[m]] <- gcode[m]
    gsize[levels(fe$term) == glv[m]] <- 2
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
    theme
  if (sd) {
    p <- p + geom_errorbar(aes(alpha = term), width = 0.2) +
      scale_alpha_manual("Group", breaks = levels(fe$term), values = gltp)
  }
  if(!missing(modinf)) {
    p <- p + geom_text(x = 2, y = 5, label = "abc")
  }
  p
}

multiGrowthMod <- function(dat, fact, SplMonthlv, grouplv, glv, gcode, 
                           tag = NULL, suffix = "", archiplot = TRUE, log = TRUE) {
  library(merTools);library(lsmeans);library(multcompView);library(car);library(ggplot2);library(lattice)
  fe.g <- NULL; re.g <- NULL; modinf.g <- NULL; shapiro.g <- NULL; posthoc.g <- NULL; modavo.g <- NULL; 
  resid.g <- NULL;
  theme_HL <- theme_bw() + theme(legend.position = "none",axis.text = element_text(angle = 30))
  if(missing(tag)) tag <- unique(dat$traits)
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
    modfact <- strsplit(as.character(formula(mod@call))[3],"+", fixed = TRUE)[[1]]
    if (grepl('*', modfact[1], fixed = TRUE)) {
      modfact <- strsplit(modfact[1],"*", fixed = TRUE)[[1]]
    }
    modfxfact <-  modfact[!grepl('(', modfact, fixed = TRUE)]
    modfxfm <- if(length(modfxfact) == 1) {modfxfact[1]} else {paste(modfxfact[1],modfxfact[2],sep = "+")}
    posthoc <-lsmeans::cld(lsmeans(object = mod, adjust = "tukey",
                                   specs = as.formula(paste("pairwise~",modfxfm))),
                           Letters = LETTERS) 
    print(posthoc) 
    if(length(modfxfact) == 1) {
      if (modfxfact == "group ") {
        posthoc <- cbind(posthoc, SplMonth = NA)
      } else {
        if (modfxfact == "SplMonth ") {
          posthoc <- cbind(posthoc, group = NA)
        } else {stop("Kidding me?")}
      }
    }
    posthoc.g <- rbind(posthoc.g,cbind(posthoc, tag[i]))
    
    # fe
    fe <- FEsim(mod)
    fe.g <- rbind(fe.g, cbind(tag = tag[i],fe))
    ggsave(plot = plotFEsim2(fe%>% filter(term != "(Intercept)") %>%
                               mutate(term = gsub("group","",term)) %>% 
                               mutate(term = gsub("SplMonth","",term)),
                             glv = glv, gcode = gcode, theme = theme_HL), 
           paste("growth/plot/",tag[i],"_Fixeff.png",sep=""),
           width = 6, height = 4)
    #re
    re <- REsim(mod)
    re.g <- rbind(re.g, cbind(tag = tag[i],re))
    ggsave(plot = merTools::plotREsim(re), 
           paste("growth/plot/",tag[i],"_Raneff.png",sep=""),
           width = 6, height = 4)
    #resid
    png(paste("growth/plot/",tag[i],"_diag_resid.png",sep=""), 
        width = 30, height = 20, units = "cm", res = 600)
    print(plot(mod, type = c("p", "smooth")))
    dev.off()
    png(paste("growth/plot/",tag[i],"_diag_residQQ.png",sep=""), 
        width = 30, height = 20, units = "cm", res = 600)
    print(qqmath(mod, id = 0.05))
    dev.off()
    #png(paste("growth/plot/",tag[i],"_diag_zeta.png",sep=""), 
    #    width = 30, height = 20, units = "cm", res = 600)
    #something using xyplot
    #dev.off()
    #png(paste("growth/plot/",tag[i],"_diag_dens.png",sep=""), 
    #    width = 30, height = 20, units = "cm", res = 600)
    #something using densityplot
    #dev.off()
    #png(paste("growth/plot/",tag[i],"_diag_pair.png",sep=""), 
    #    width = 30, height = 20, units = "cm", res = 600)
    #something using splom
    #dev.off()
    stats::shapiro.test(resid(mod)) -> resid
    resid.g <- rbind(resid.g, 
                     data.frame(tag = tag[i],
                                W = resid$statistic,
                                Pvalue = resid$p.value))
  }
  ## output
  if (log) {write.csv(x = fe.g, file = paste("growth/log/FixedEff",suffix,".csv",sep = ""), row.names = F)
    write.csv(x = re.g, file = paste("growth/log/RandomEff",suffix,".csv",sep = ""), row.names = F)
    write.csv(x = modinf.g, file = paste("growth/log/ModelChoice",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = shapiro.g, file = paste("growth/log/ShapiroRawData",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = posthoc.g, file = paste("growth/log/Posthoc",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = modavo.g, file = paste("growth/log/ModelAnova",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = resid.g, file = paste("growth/log/ShapiroResid",suffix,".csv",sep = ""), row.names = F)
  }
  "DONE!"
}



## Basic Stat information ----
meanseCal(datareadln())

## LMM fiting and ploting ----
multiGrowthMod(dat = datareadln() %>% gather(key = traits, value = resp, Density:SeedRate),
                fact = c("group + (1|SiteID)",
                         "SplMonth + (1|SiteID)",
                         "group + (1|SiteID:SplMonth)",
                         "SplMonth + (1|SiteID:SplMonth)",
                         "group + (1|SiteID) + (1|SiteID:SplMonth)",
                         "SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                         "group + SplMonth + (1|SiteID)",
                         "group + SplMonth + (1|SiteID:SplMonth)",
                         "group + SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                         "group * SplMonth + (1|SiteID)",
                         "group * SplMonth + (1|SiteID:SplMonth)",
                         "group * SplMonth + (1|SiteID) + (1|SiteID:SplMonth)"),
                SplMonthlv = c("Apr","Jul","Nov"), grouplv = c("EA","CL","WE"), 
                glv = c("EA","WE"), gcode = c("#31B404","#013ADF"),
               tag = c("Density","Height","BsDmr","TlrWg"), suffix = "_ajn")

multiGrowthMod(dat = datareadln() %>% gather(key = traits, value = resp, Density:SeedRate),
               fact = c("group + (1|SiteID)",
                        "SplMonth + (1|SiteID)",
                        "group + (1|SiteID:SplMonth)",
                        "SplMonth + (1|SiteID:SplMonth)",
                        "group + (1|SiteID) + (1|SiteID:SplMonth)",
                        "SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                        "group + SplMonth + (1|SiteID)",
                        "group + SplMonth + (1|SiteID:SplMonth)",
                        "group + SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                        "group * SplMonth + (1|SiteID)",
                        "group * SplMonth + (1|SiteID:SplMonth)",
                        "group * SplMonth + (1|SiteID) + (1|SiteID:SplMonth)"),
               SplMonthlv = c("Apr","Jul","Nov"), grouplv = c("EA","CL","WE"), 
               glv = c("EA","WE"), gcode = c("#31B404","#013ADF"),
               tag = c("LvThk","LvCt","LvWgG","StWg"), suffix = "_jn")

multiGrowthMod(dat = datareadln() %>% gather(key = traits, value = resp, Density:SeedRate),
               fact = c("group + (1|SiteID)",
                        "SplMonth + (1|SiteID)",
                        "group + (1|SiteID:SplMonth)",
                        "SplMonth + (1|SiteID:SplMonth)",
                        "group + (1|SiteID) + (1|SiteID:SplMonth)",
                        "SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                        "group + SplMonth + (1|SiteID)",
                        "group + SplMonth + (1|SiteID:SplMonth)",
                        "group + SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                        "group * SplMonth + (1|SiteID)",
                        "group * SplMonth + (1|SiteID:SplMonth)",
                        "group * SplMonth + (1|SiteID) + (1|SiteID:SplMonth)"),
               SplMonthlv = c("Apr","Jul","Nov"), grouplv = c("EA","CL","WE"), 
               glv = c("EA","WE"), gcode = c("#31B404","#013ADF"),
               tag = c("LvCtY","SeedRate"), suffix = "_n")

## Gather ploting ----
ggsave(plot = plotFEsim2facet(read.csv("growth/log/FixedEff_ajn.csv") %>% 
                                filter(term != "(Intercept)") %>% 
                                mutate(term = gsub("group","",term)) %>% 
                                mutate(term = gsub("SplMonth","",term)),
                              glv = c("EA","WE"), gcode = c("#31B404","#013ADF"), ncol = 2,
                              theme = theme_bw() + theme(legend.position = "none",
                                                         axis.text = element_text(size= 4,angle = 30),
                                                         axis.title = element_text(size= 6),
                                                         strip.text = element_text(size= 6))), 
       "growth/plot/Fixeff_ajn.png",
       width = 6, height = 4, dpi = 600)

ggsave(plot = plotFEsim2facet(read.csv("growth/log/FixedEff_jn.csv") %>% 
                                filter(term != "(Intercept)") %>% 
                                mutate(term = gsub("group","",term)) %>% 
                                mutate(term = gsub("SplMonth","",term)),
                              glv = c("EA","WE"), gcode = c("#31B404","#013ADF"), ncol = 2,
                              theme = theme_bw() + theme(legend.position = "none",
                                                         axis.text = element_text(size= 4,angle = 30),
                                                         axis.title = element_text(size= 6),
                                                         strip.text = element_text(size= 6))), 
       "growth/plot/Fixeff_jn.png",
       width = 6, height = 4, dpi = 600)

ggsave(plot = plotFEsim2facet(read.csv("growth/log/FixedEff_n.csv") %>% 
                                filter(term != "(Intercept)") %>% 
                                mutate(term = gsub("group","",term)) %>% 
                                mutate(term = gsub("SplMonth","",term)),
                              glv = c("EA","WE"), gcode = c("#31B404","#013ADF"), ncol = 2,
                              theme = theme_bw() + theme(legend.position = "none",
                                                         axis.text = element_text(size= 4,angle = 30),
                                                         axis.title = element_text(size= 6),
                                                         strip.text = element_text(size= 6))), 
       "growth/plot/Fixeff_n.png",
       width = 6, height = 3, dpi = 600)


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





