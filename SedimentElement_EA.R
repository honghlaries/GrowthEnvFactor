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
  read.csv("./Data/Result_Sediment.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    dplyr::filter(group != "NV") %>% 
    dplyr::select(N:Pb,SiteID,SplMonth,group) 
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
  ggsave(filename = paste("sediment/plot/bar_",tag,"_EA.png",sep = ""), width = 6, height = 4)
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
  dplyr::select(SiteID, SplMonth, N:Pb)  %>% 
  dplyr::mutate(C_N = orgC / N) %>%
  tidyr::gather(trait, value, N:C_N) 

ggplot(data = dat%>%
         dplyr::group_by(SiteID,SplMonth,trait) %>%
         dplyr::summarise(mean = mean(value,na.rm = T),se = (sd(value,na.rm = T)/sum(!is.na(value))))) + 
  geom_bar(aes(x = SplMonth, y = mean, fill = SiteID),stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(x = SplMonth, ymin = mean - se, ymax = mean + se, fill = SiteID),
                position = position_dodge(0.9), width = 0.3) +
  xlab("") + ylab("") +
  scale_fill_discrete("Location",breaks = c("EA1","EA2"),labels = c("Seaward(EA1)","Landward(EA2)")) + 
  facet_wrap(facets = ~ trait, nrow = 4, scales = "free") 
ggsave(filename = "sediment/plot/bar_SedimentElement_EA.png")

t.summary <- NULL
tagList <- unique(dat$trait)
for(i in 1: length(tagList)){
  t.summary <- rbind(t.summary,doTandPlot(dat,tagList[i]))
} %>% write.csv("sediment/log/t.test_EA.csv",row.names = F)

read.csv("sediment/log/t.test_EA.csv") %>% print(row.names = F)

## ---------------------------------------------------
## scaling and re-analysis
##
## ---------------------------------------------------

# Al
tag <- "Al"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# As
tag <- "As"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# C
tag <- "C"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
dat1$value =(atan( ((dat1$value - 16400)/2000) ))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# C:N
tag <- "C_N"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = atan(sqrt(value-7)-1.4) )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Cd
tag <- "Cd"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Cr
tag <- "Cr"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Cu
tag <- "Cu"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Fe
tag <- "Fe"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Mn
tag <- "Mn"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = (value))
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# N
tag <- "N"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = atan( (value-650)/100) )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Ni
tag <- "Ni"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = value )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# orgC
tag <- "orgC"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value = atan( (value-6000)/200 ) )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# P
tag <- "P"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value =  atan( (value-550)/200 ) )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Pb
tag <- "Pb"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value =  value )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# S
tag <- "S"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value =  value )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)

# Zn
tag <- "Zn"

dat1 <- dat %>% filter(trait == tag) %>% 
  filter(!is.na(value)) %>% 
  mutate(SplMonth = as.character(SplMonth),
         value =  value )
shapiroTest(dat1)
leveneTest(dat1)

doadjustedTest(dat1,tag)
