#fit and arbutus for just BM

library(ape)
library(geiger)
library(phytools)
library(tidyverse)
library(arbutus)

CleanData <- function(phy, data) {
  cleaned <- treedata(phy,data,warnings=T)
  return(cleaned)
}

tree.file <- "Data/fungi_tree.tre"
gene.exp <- "Data/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv"

exp.df<- t(read.table(gene.exp,sep="\t",header=T,row.names = 1,stringsAsFactors = F))
tree <- read.tree(tree.file)
cleaned <- CleanData(phy = tree,exp.df)

dat <- cleaned$data

arbutusBM <- vector(mode = "list", length = ncol(dat))
for(j in 1:ncol(dat)){
  tmp <- cleaned$phy
  todrop <- c()
  for(i in 1:nrow(dat)){
    if(is.na(dat[[i,j]])){
      todrop <- append(todrop, rownames(dat)[i])
    }
  }
  tmp <- drop.tip(tmp, todrop)
  trimmed <- treedata(tmp, dat[,j], sort = TRUE)
  fitBM <- fitContinuous(trimmed$phy, trimmed$data, model = "BM")
  arbutusBM[[j]] <- arbutus(fitBM)
}

BM_pvals <- map_df(arbutusBM, pvalue_arbutus)

saveRDS(BM_pvals, "BM_pvals")
BM_piv <- BM_pvals %>% pivot_longer(cols = everything(), names_to = "tstat")
BM_piv %>% ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
ggsave("arbutus_BM.png")