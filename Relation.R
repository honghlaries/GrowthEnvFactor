## clean ----
rm(list = ls())

## package loading ----
library(MASS)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

## functions ----
datareadln <- function() {
  library(dplyr)
  library(tidyr)
  dat <- inner_join(x = read.csv("./Data/Result_Sediment.csv") %>%
                      inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
                      inner_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
                      inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>% 
                      mutate(PhyStress = group=="WE") %>%
                      select(QudNo,N,C,S,P,orgC,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,SplMonth,PhyStress,group,SiteID),
                    y = read.csv("./Data/Result_PlantTillerGrowthSummary.csv") %>%
                      mutate(density = 4 * Number)%>%
                      inner_join(read.csv("./Data/meta_PlantSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
                      inner_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
                      select(QudNo,density,Seed_rate,LvCtY_mean,LvCtY_max,LvCtG_mean,LvCtG_max,LvCt_mean,LvCt_max,
                             Height_mean,Height_max,BsDmr_mean,BsDmr_max,LvThk_mean,LvThk_max,
                             TlrWg_mean,TlrWg_max,LvWgG_mean,LvWgG_max,StWg_mean,StWg_max),
                    by = c("QudNo" = "QudNo"))
}

datTran <- function(dat) {
  dat1 <- dat%>% 
    filter(group == "CL") %>%
    group_by(SplMonth) %>%
    summarise(Hmean = mean(Height_mean),Hsd = sd(Height_mean),
              Dmean = mean(BsDmr_mean),Dsd = sd(BsDmr_mean),
              Bmean = mean(TlrWg_mean),Bsd = sd(TlrWg_mean),
              denmean = mean(Height_mean),densd = sd(Height_mean))
  dat %>%
    inner_join(dat1, by = c("SplMonth" = "SplMonth"))%>%
    mutate(Height = log(Height_mean / Hmean), 
           Diameter = log(BsDmr_mean / Dmean), 
           AbgBiomass = log(TlrWg_mean / Bmean),
           RametDensity = log(density / denmean))%>%     
    #mutate(Height = (Height_mean - Hmean)/Hsd, 
    #       Diameter = (BsDmr_mean - Dmean)/Dsd, 
    #       AbgBiomass = (TlrWg_mean - Bmean)/Bsd,
    #       RametDensity = (density - denmean)/densd)%>% 
    select(SiteID,SplMonth,PhyStress,group,
           N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,
           Height,Diameter,AbgBiomass,RametDensity)
}

StepRDA <- function(Growth,Env,Tag,SplMonth) {
  library(vegan)
  library(ggplot2)
  
  null <- rda(Growth~1,Env,scale=TRUE) 
  full <- rda(Growth~.,Env,scale=TRUE) 
  mod <- step(null, scope = formula(full), test = "perm")
  
  perm <- anova.cca(mod,permutations = how(nperm=9999))
  effload <- mod$CCA$v[,1:2]; effloadx <-effload - effload; efftag <- row.names(effload); effload <- data.frame(rbind(effloadx,effload),efftag)
  envload <- mod$CCA$biplot[,1:2]; envloadx <- envload - envload; envtag <- row.names(envload); envload <- data.frame(rbind(envloadx,envload),envtag)
  sampload <- mod$CCA$wa[,1:2]; sampload <- data.frame(sampload,Tag,SplMonth)
  sampload.site <- sampload %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("Tag" = "SiteID"))
  sampload.group <- sampload %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("Tag" = "SiteID")) %>%
    group_by(SplMonth, group) %>%
    dplyr::summarise(RDA1.avg = mean(RDA1,na.rm = T),
              RDA1.se = sd(RDA1,na.rm = T)/sqrt(n()-sum(is.na(RDA1))),
              RDA1.sd = sd(RDA1,na.rm = T),
              RDA2.avg = mean(RDA2,na.rm = T),
              RDA2.se = sd(RDA2,na.rm = T)/sqrt(n()-sum(is.na(RDA2))),
              RDA2.sd = sd(RDA2,na.rm = T))
  loadplot <- ggplot() +
    geom_point(aes(x = RDA1,y = RDA2, col = group, shape = SplMonth), size = 1,
               data = sampload.site) +
    geom_path(aes(x = RDA1,y = RDA2),group = envload$envtag, size = 0.7,
              data = envload, col = "black") +
    geom_path(aes(x = RDA1,y = RDA2),group = effload$efftag, size = 0.7,
              data = effload, col = "blue") +
    geom_errorbar(aes(x = RDA1.avg, y = RDA2.avg, ymax = RDA2.avg + RDA2.sd, ymin = RDA2.avg - RDA2.sd), 
                      col = "grey50", size = 0.7, data = sampload.group) +
    geom_errorbarh(aes(y = RDA2.avg, x = RDA1.avg, xmax = RDA1.avg + RDA1.sd, xmin = RDA1.avg - RDA1.sd), 
                   col = "grey50", size = 0.7, data = sampload.group) +
    geom_point(aes(x = RDA1.avg,y = RDA2.avg, col = group, shape = SplMonth), size = 3, 
               data = sampload.group) +
    geom_label(aes(x = RDA1,y = RDA2, label = envtag), size = 3.5, data = envload[7:12,], col = "black") + 
    geom_label(aes(x = RDA1,y = RDA2, label = efftag), size = 3.5, data = effload[5:8,], col = "blue", 
               nudge_y =c(0.05,0.02,-0.07,-0.05),nudge_x =c(0,-0.1,-0,0)) +
    scale_x_continuous(name = "RDA1", limits = c(-1.1,1.1)) +
    scale_y_continuous(name = "RDA2", limits = c(-1.1,1.1)) +
    scale_shape("Month") + scale_color_hue("Location") +  
    theme_bw()
  list(mod = mod, plot = loadplot, perm = perm)
}

RDAlmfit <- function(Growth,selG,Env,selE,group,SplMonth) {
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  Growth <- Growth %>% select(one_of(selG))
  Env <- Env %>% select(one_of(selE))
  
  dat.p <- cbind(Growth,Env,group,SplMonth) %>%
    tidyr::gather_("Etag","Env",selE) %>%
    tidyr::gather_("Gtag","Growth",selG)
  
  dat.smooth.gather <- NULL
  dat.smooth.group <- NULL
  for(i in 1:length(selE)) { 
    for(j in 1:length(selG)) {
      dat.tmp <- dat.p %>% filter(Etag == selE[i],Gtag == selG[j])
      if(cor.test(dat.tmp$Growth,dat.tmp$Env,method="pearson")$p.value < 0.05) {
        dat.smooth.gather <- rbind(dat.smooth.gather,dat.tmp)
      }
      for(k in 1:length(group)) {
        dat.tmp <- dat.p %>% filter(Etag == selE[i],Gtag == selG[j],group == group[k]) 
        if(cor.test(dat.tmp$Growth,dat.tmp$Env,method="pearson")$p.value < 0.05) {
          dat.smooth.group <- rbind(dat.smooth.gather,dat.tmp)
        }   
      }
    }
  }

  cor <- NULL 
  for(i in 1:length(colnames(Growth))) {
    for(j in 1:length(colnames(Env))) {
      Etag <- colnames(Env)[j]
      Gtag <- colnames(Growth)[i]
      relation <- cor.test(Growth[,i],Env[,j],method="pearson")
      cor <- rbind(cor,data.frame(Etag,Gtag,
                                  paste("Pearson's Cor = ",
                                        format(relation$estimate, digit = 2),
                                        "\n","p = ",
                                        format(relation$p.value, digit = 2), 
                                        sep = ""),
                                  x = (max(Env[,j]) + min(Env[,j]))/2, y = 0.85*max(Growth[,i])+0.15*min(Growth[,i])))
    }
  }
  colnames(cor) <- c("Etag","Gtag","label", "x", "y")
  
  ggplot() +
    geom_point(aes(x = Env, y = Growth, col = group, shape = SplMonth), data = dat.p) +
    geom_smooth(aes(x = Env, y = Growth),col = "black", method = rlm, se = F,
                data = dat.smooth.gather) +
    geom_smooth(aes(x = Env, y = Growth,  col = group), method = rlm, se = F,
                data = dat.smooth.group) +
    facet_grid(Gtag~Etag, scales = "free") +
    geom_text(aes(label = label, x= x, y = y),data = cor, size = 2.5) 
}

## DCA ----
result.dca <-datareadln()%>% 
  #select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,
  #       Height,Diameter,AbgBiomass,RametDensity)%>% 
  #filter(SplMonth != "Apr") %>%
  select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,
         density,Height_mean,BsDmr_mean,TlrWg_mean)%>% 
  decorana() 
anova.cca(result.dca)

## RDA ----
library(ggplot2)
library(vegan)
library(MASS)
data <- datareadln() %>% datTran()
Env <-  data %>% select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb)
Growth <-data %>%  select(Height,Diameter,AbgBiomass,RametDensity)
Tag <- data$SiteID
SplMonth <- data$SplMonth
group <- data$group

Result.rda <- StepRDA(Growth,Env,Tag,SplMonth)
Result.rda$mod
Result.rda$perm
Result.rda$plot


qplot(data = data, x = N, y = AbgBiomass) +
  geom_smooth(method = rlm, linetype = 1, col = "black",fill= "red") +
  geom_smooth(method = lm, linetype = 2, col = "black", fill ="blue") +
  theme_bw()

qplot(data = data, x = As, y = Diameter) +
  geom_smooth(method = rlm, linetype = 1, col = "black",fill= "red") +
  geom_smooth(method = lm, linetype = 2, col = "black", fill ="blue") +
  theme_bw()

selG <- c("Height", "Diameter", "AbgBiomass", "RametDensity")
selE <- c("N", "S", "Cu", "Mn", "Ni", "As")

RDAlmfit(Growth,selG,Env,selE,group,SplMonth)