## ---------------------------------------------------
## init
##
## ---------------------------------------------------
rm(list = ls())
library(dplyr);library(tidyr);library(ggplot2);library(car)

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
    dplyr::select(group,Pb,Cd,As,Cr,Cu,Zn,Ni) 
}

bkreadln <- function() {
  read.csv("./Data/meta_background.csv")
}

## ---------------------------------------------------
## Calculation 
##
## ---------------------------------------------------
dat <- datareadln()
bk <- bkreadln() 

dat.sum <- dat %>% 
  tidyr::gather(element,value,Pb:Ni) %>%
  dplyr::group_by(group,element) %>%
  dplyr::summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                                 format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                                 sep = " ± ")) %>%
  tidyr::spread(group, value) %>%
  dplyr::rename(CL.aft = CL, EA.aft  = EA)  %>% dplyr::select(-WE)

bk.sum <- bk %>% 
  tidyr::gather(element,value,Pb:Ni) %>%
  dplyr::group_by(group,element) %>%
  dplyr::summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                                 format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                                 sep = " ± ")) %>%
  tidyr::spread(group, value) %>%
  dplyr::rename(CL.bef = CL, EA.bef  = EA)

dat.adjust <- dat %>% 
  tidyr::gather(element,value,Pb:Ni) %>%
  dplyr::inner_join(bk %>% tidyr::gather(element,value,Pb:Ni) %>%
                      dplyr::group_by(group,element) %>%
                      dplyr::summarise(bk = mean(value,na.rm = T)),
                    by = c("element" = "element", "group" = "group")) %>%
  dplyr::mutate(value = value/bk) %>% dplyr::select(-bk) 

dat.adj.sum <- dat.adjust%>%
  dplyr::group_by(group,element) %>%
  dplyr::summarise(value = paste(format(mean(value,na.rm = T),digit = 3), 
                                 format(sd(value,na.rm = T)/sqrt(n()-sum(is.na(value))),digit = 3), 
                                 sep = " ± ")) %>%
  tidyr::spread(group, value) %>%
  dplyr::rename(CL.ratio = CL, EA.ratio  = EA) 


tags <- unique(dat.adjust$element)

rst <- NULL
for(i in 1: length(tags)) {
  dat1 <- dat.adjust %>% dplyr::filter(element == tags[i])
  t <- t.test(value~group,data = dat1)
  if (t$p.value >= 0.05) p <- format(t$p.value,digit = 3) else 
    if (t$p.value >= 0.01) p <- "*" else 
      if (t$p.value >= 0.001) p <- "**" else 
        p <- "***"
  rst <- rbind(rst, data.frame(element = tags[i],t = t$statistic,p = p))
}

output <- bk.sum %>%
  dplyr::inner_join(dat.sum) %>%
  dplyr::inner_join(dat.adj.sum) %>%
  dplyr::inner_join(rst) 

write.csv(output,"sediment/log/comparation.csv")