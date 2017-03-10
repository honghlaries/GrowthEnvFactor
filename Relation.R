## clean ----
rm(list = ls())

## package loading ----
library(MASS)
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

StepRDA <- function(Growth,Env,Tag,SplMonth,path = "./GrowthvsElement/",filename) {
  library(vegan)
  library(ggplot2)
  
  null <- rda(Growth~1,Env,scale=TRUE) 
  full <- rda(Growth~.,Env,scale=TRUE) 
  mod <- step(null, scope = formula(full), test = "perm")
  
  perm <- anova.cca(mod,permutations = how(nperm=9999))
  print(perm)
  write.csv(format(perm),paste(path, "output_", filename, ".csv", sep = ""))
  
  plot(mod)
  effload <- mod$CCA$v[,1:2]; effloadx <-effload - effload; efftag <- row.names(effload); effload <- data.frame(rbind(effloadx,effload),efftag)
  envload <- mod$CCA$biplot[,1:2]; envloadx <- envload - envload; envtag <- row.names(envload); envload <- data.frame(rbind(envloadx,envload),envtag)
  sampload <- mod$CCA$wa[,1:2]; sampload <- data.frame(sampload,Tag,SplMonth)
  sampload.site <- sampload %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("Tag" = "SiteID"))
  sampload.group <- sampload %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("Tag" = "SiteID")) %>%
    select(-col) %>%
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
  ggsave(plot = loadplot,filename = paste(path, "RDAloading_", filename, ".eps", sep = ""))
  ggsave(plot = loadplot,filename = paste(path, "RDAloading_", filename, ".tiff", sep = ""))
  ggsave(plot = loadplot,filename = paste(path, "RDAloading_", filename, ".png", sep = ""), dpi = 1200)
  mod
}

Rlmfit <- function(Growth,selG,Env,selE,group,SplMonth,path = "./GrowthvsElement/",filename) {
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  Growth <- Growth %>% select(one_of(selG))
  Env <- Env %>% select(one_of(selE))
  
  dat <- NULL
  for(i in 1:length(colnames(Growth))) {
    for(j in 1:length(colnames(Env))) {
      Etag <- colnames(Env)[j]
      Gtag <- colnames(Growth)[i]
      relation <- cor.test(Growth[,i],Env[,j],method="pearson")
      dat <- rbind(dat,data.frame(Etag,Gtag,Env[,j],Growth[,i],group,SplMonth,Sign = (relation$p.value < 0.05)))
    }
  }
  colnames(dat) <- c("Etag","Gtag","Env","Growth","group","SplMonth","Sign") 
  
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
  
  smatrix <- ggplot(data = dat, aes(x = Env, y = Growth)) +
    geom_point(aes(col = group, shape = SplMonth, alpha = Sign)) +
    geom_smooth(aes(linetype = !Sign, col = group),method = rlm, fill= "grey50",se = F, size = 1) +
    geom_smooth(aes(linetype = !Sign),method = rlm, col = "black",fill= "grey50") +
    facet_grid(Gtag~Etag, scales = "free") +
    geom_text(aes(label = label, x= x, y = y),data = cor, size = 2.5) +
    scale_x_continuous(name = "Element Content (mg/kg)") +
    scale_y_continuous(name = "Growth Promotion Effect") +
    scale_color_discrete(name = "Location") +
    scale_shape(name = "Month") +
    scale_alpha_discrete(range = c(0.00,1),breaks = c(0,1), labels = (c(">= 0.05", "< 0.05")), name = "Significance") +
    scale_linetype_discrete(name = "Significance",breaks = c(0,1), labels = (c("< 0.05", ">= 0.05"))) +
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
  ggsave(plot = smatrix, filename = paste(path,"Matrix_",filename,".png", sep = ""), dpi = 1200)
  return(smatrix)
}

## DCA ----
result.dca <-datareadln()%>% 
  select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,
         Height,Diameter,AbgBiomass,RametDensity)%>% 
  #filter(SplMonth != "Apr") %>%
  #select(density,Seed_rate,LvCtY_mean,LvCtY_max,
  #       Height_mean,Height_max,BsDmr_mean,BsDmr_max,LvThk_mean,LvThk_max,
  #       TlrWg_mean,TlrWg_max,LvWgG_mean,LvWgG_max,StWg_mean,StWg_max)%>% 
  decorana() 
anova.cca(result.dca)

## RDA ----
library(ggplot2)
library(vegan)
library(MASS)
data <- datareadln()
Env <-  data %>% select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb)
Growth <-data %>%  select(Height,Diameter,AbgBiomass,RametDensity)
Tag <- data$SiteID
SplMonth <- data$SplMonth
group <- data$group

Result.rda <- StepRDA(Growth,Env,Tag,SplMonth, filename = "LogTrans")

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

Rlmfit(Growth,selG,Env,selE,group,SplMonth,filename = "selelement")

## arc ====
#Apr-all
data <- datareadln()%>% 
  filter(SplMonth == "Apr") 

Env <-  data %>%
  select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,PhyStress)
Growth <-data %>%
  select(density,Height_mean,BsDmr_mean,TlrWg_mean)

Result.rda$null <- rda(Growth~1,Env,scale=TRUE) 
Result.rda$full <- rda(Growth~.,Env,scale=TRUE) 
Result.rda$model<- step(Result.rda$null, scope = formula(Result.rda$full), test = "perm")
anova.cca(Result.rda$model,permutations = how(nperm=9999))
plot(Result.rda$model)

qplot(data = data, x = orgC, y = TlrWg_mean, col = group)
qplot(data = data, x = Al, y = density, col = group)

#Apr-ref
data <- datareadln()%>% 
  filter(SplMonth == "Apr") %>%
  filter(group != "EA")

Env <-  data %>%
  select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,PhyStress)
Growth <-data %>%
  select(density,Height_mean,BsDmr_mean,TlrWg_mean)

Result.rda$null <- rda(Growth~1,Env,scale=TRUE) 
Result.rda$full <- rda(Growth~.,Env,scale=TRUE) 
Result.rda$model<- step(Result.rda$null, scope = formula(Result.rda$full),test = "perm")
anova.cca(Result.rda$model,permutations = how(nperm=9999))
plot(Result.rda$model)

#Jul-all
data <- datareadln()%>% 
  filter(SplMonth == "Jul") 

Env <-  data %>%
  select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,PhyStress)
Growth <-data %>%
  select(#density,
         LvCtG_mean,
         Height_mean,BsDmr_mean,LvThk_mean,
         TlrWg_mean,StWg_mean,LvWgG_mean)

Result.rda <- NULL
Result.rda$null <- rda(Growth~1,Env,scale=TRUE) 
Result.rda$full <- rda(Growth~.,Env,scale=TRUE) 
Result.rda$model<- step(Result.rda$null, scope = formula(Result.rda$full), test = "perm")
anova.cca(Result.rda$model,permutations = how(nperm=9999))
plot(Result.rda$model)

Result.glm <- NULL
Result.glm$null <- lm(density~1,data,scale=TRUE) 
Result.glm$full <- lm(density~N + S +P +Al +Cr +Mn +Fe +Ni +Cu +Zn +As +Cd +Pb + PhyStress, data ,scale=TRUE)
Result.glm$model<- step(Result.glm$null, scope = formula(Result.glm$full))

qplot(data = data, x = orgC, y = TlrWg_mean, col = group)
qplot(data = data, x = As, y = density, col = group)

#Jul-ref
data <- (datareadln()%>% 
  filter(SplMonth == "Jul") %>%
  filter(group != "EA"))[-1:-3,]

Env <-  data %>%
  select(N,C,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,PhyStress)
Growth <-data %>%
  select(density,LvCtG_mean,
         Height_mean,BsDmr_mean,LvThk_mean,
         TlrWg_mean,StWg_mean,LvWgG_mean)

Result.rda$null <- rda(Growth~1,Env,scale=TRUE) 
Result.rda$full <- rda(Growth~.,Env,scale=TRUE) 
Result.rda$model<- step(Result.rda$null, scope = formula(Result.rda$full), test = "perm")
anova.cca(Result.rda$model,permutations = how(nperm=9999))
plot(Result.rda$model)

qplot(data = data, x = orgC, y = TlrWg_mean, col = group)
qplot(data = data, x = N, y = density, col = group)


#Nov-all
data <- datareadln()%>% 
  filter(SplMonth == "Nov") 

Env <-  data %>%
  select(N,S,P,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,PhyStress)
Growth <-data %>%
  select(density,LvCtG_mean,LvCtY_mean, Seed_rate,
         Height_mean,BsDmr_mean,LvThk_mean,
         TlrWg_mean,StWg_mean,LvWgG_mean)

Result.rda$null <- rda(Growth~1,Env,scale=TRUE) 
Result.rda$full <- rda(Growth~.,Env,scale=TRUE) 
Result.rda$model<- step(Result.rda$null, scope = formula(Result.rda$full), test = "perm")
anova.cca(Result.rda$model,permutations = how(nperm=9999))
plot(Result.rda$model)

qplot(data = data, x = orgC, y = TlrWg_mean, col = group)
qplot(data = data, x = Cu, y = TlrWg_mean, col = group)


#Nov-ref
data <- datareadln()%>% 
  filter(SplMonth == "Nov") %>%
  filter(group != "EA")

Env <-  data %>%
  select(N,S,P,orgC,Al,Cr,Mn,Fe,Ni,Cu,Zn,As,Cd,Pb,PhyStress)
Growth <-data %>%
  select(density,LvCtG_mean,LvCtY_mean, Seed_rate,
         Height_mean,BsDmr_mean,LvThk_mean,
         TlrWg_mean,StWg_mean,LvWgG_mean)

Result.rda$null <- rda(Growth~1,Env,scale=TRUE) 
Result.rda$full <- rda(Growth~.,Env,scale=TRUE) 
Result.rda$model<- step(Result.rda$null, scope = formula(Result.rda$full), test = "perm")
anova.cca(Result.rda$model,permutations = how(nperm=9999))
plot(Result.rda$model)


