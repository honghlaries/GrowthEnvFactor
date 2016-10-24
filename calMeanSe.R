calMeanSe <- function  (dat, group, trait.input, trait.output = trait.input) {
  library(dplyr)
  library(tidyr)
  library(lazyeval)
  #dots = list(interp(~ mean(var, na.rm = T), var = as.name(trait.input)),
  #            interp(~ mean(var, na.rm = T), var = as.name(trait.input)))
  dots <- list(~mean(Number), ~n)
  mean <- dat %>% 
    group_by_(group) %>%
    summarise_(.dots = dots)
  #colnames(mean) <- paste(trait.output,c("_mean","_se"), sep = "")
  
  return(mean)
}

