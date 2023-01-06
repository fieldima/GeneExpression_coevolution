#Temporary BMS script cause of library issue

#Get fitContinuous and plot AIC of each protein

library(ape)
library(geiger)
library(phytools)
library(tidyverse)
library(arbutus)
library(OUwie)
library(motmot)

print("Coevolution")

class_fix <- function(fitobj){
  if(length(fitobj) == 5) class(fitobj) <- "gfit"
  if(length(fitobj) == 27) class(fitobj) <- "OUwie"
  fitobj
}

run_arb <- function (fits){
  arby <- lapply(fits, function(i) tryCatch(arbutus(i), error = function(x) NA))
  arby <- arby[!is.na(arby)]
  arby
}

fitResults <- readRDS("BMS_fits")

arb_result <- map(fitResults, class_fix) %>% run_arb()

pvals <- map_df(arb_result, pvalue_arbutus)
saveRDS(pvals, file = "BMS_pvals")

p_piv <- pvals %>% pivot_longer(cols = everything(), names_to = "tstat")

p_piv %>% ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
ggsave("BMS_results.png")
