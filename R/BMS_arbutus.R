#Get fitContinuous and plot AIC of each protein

library(ape)
library(geiger)
library(phytools)
library(tidyverse)
library(arbutus)
library(OUwie)
library(motmot)

CleanData <- function(phy, data) {
  cleaned <- treedata(phy,data,warnings=T)
  return(cleaned)
}

tree.file <- "Data/fungi_tree.tre"
gene.exp <- "Data/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv"

exp.df<- t(read.table(gene.exp,sep="\t",header=T,row.names = 1,stringsAsFactors = F))
tree <- read.tree(tree.file)
cleaned <- CleanData(phy = tree,exp.df)

motmot <- readRDS("motmot1.rds")

#Make OUwie df
OUwie_reg <- c(rep(2,17), 1)
daf <- data.frame(Genus_species = tree$tip.label, Reg = (OUwie_reg))

#For each gene, trim tips off the phylogeny when that species is missing a gene
phy <- cleaned$phy
data <- cleaned$data

fitResults <- vector(mode = "list", length = ncol(data))


for(j in 1:ncol(data)){
    todrop <- c()
for(i in 1:nrow(data)){
      if(is.na(data[[i,j]])){
        todrop <- append(todrop, rownames(data)[i])
      }
    }
    tmp <- drop.tip(phy, todrop)
    trimmed <- treedata(tmp, data[,j], sort = TRUE)
    fitBM <- fitContinuous(trimmed$phy, trimmed$data, model = "BM")
    fitOU <- fitContinuous(trimmed$phy, trimmed$data, model = "OU")
    fitEB <- fitContinuous(trimmed$phy, trimmed$data, model = "EB")
    trimmed$phy$node.label <- c(1, rep(2,16))
    OUwie_df <- daf %>% filter(Genus_species %in% rownames(trimmed$data)) %>% mutate(X = trimmed$data)
    fitBMS <- tryCatch(OUwie(trimmed$phy, OUwie_df, model = "BMS"), error = function(x)list(AIC = Inf))
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]], fitBMS$AIC)
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         ifelse(min(aic) == aic[3], list(c(fitEB, model = "EB")),
                                list(c(fitBMS, model = "BMS")))))
    fitResults[j] <- fit
  }

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  bms = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, ifelse(vec$model == "EB", eb <- eb + 1, bms <- bms + 1)))
  }
  df <- data.frame(OU = ou, BM = bm, EB = eb, BMS = bms)
  b <- df %>% pivot_longer(c(OU, BM, EB, BMS), names_to = "model")
  b
}

df <- model_count(fitResults)
df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
ggsave("BMS_AIC.png")

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

saveRDS(fitResults, file = "BMS_fits")

arb_result <- map(fitResults, class_fix) %>% run_arb()

pvals <- map_df(arb_result, pvalue_arbutus)
saveRDS(pvals, file = "BMS_pvals")

p_piv <- pvals %>% pivot_longer(cols = everything(), names_to = "tstat")

p_piv %>% ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
ggsave("BMS_results.png")
