#Load libraries
library(motmot)
library(tidyverse)
library(geiger)

print("Coevolution")

CleanData <- function(phy, data) {
  cleaned <- treedata(phy,data,warnings=T)
  return(cleaned)
}

#Load data and tree
path <- "Data/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv"
path2 <- "Data/fungi_tree.tre"

data <- t(read.table(path,sep="\t",header=T,row.names = 1,stringsAsFactors = F))
tree <- read.tree(path2)
cleaned <- CleanData(phy = tree,data)


#Motmot data
motmot <- cleaned$data %>% t() %>% as.data.frame() %>% drop_na() %>% t() %>% as.matrix() 

#Find rate shifts
tm2 <- transformPhylo.ML(motmot[,sample(1:ncol(motmot), 15)], tree, model = "tm2", minCladeSize = 3, nSplits = 2)
sum_tm2 <- summary(tm2)
plot(sum_tm2)
cutOff <- calcCutOff(tree, n = 1000, mc.cores = 16, model = "tm2")
res1 <- summary(tm2, cutoff = cutOff)
saveRDS(res1, "motmot1.rds")