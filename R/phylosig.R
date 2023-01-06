library(ape)
library(geiger)
library(phytools)
library(corrplot)
library(tidyverse)

gene.exp <- "Data/gene_expression_tpm_matrix_updated_Standard.LogNorm.tsv"
exp.df<- t(read.table(gene.exp,sep="\t",header=T,row.names = 1,stringsAsFactors = F))
tree <- read.tree(tree.file)
cleaned <- CleanData(phy = tree,exp.df)

dat <- cleaned$data
pvals <- readRDS("pvalues_df")

#Now checking the phylogenetic signal of one set of genes and comparing to other
num <- 1:dim(pvals)[1]
only_NA <- pvals %>% mutate(num = num) %>% filter(is.na(s.hgt)) %>% pull(num)
not_NA <- pvals %>% mutate(num = num) %>% filter(!is.na(s.hgt)) %>% pull(num)

only_NA_df <- dat[,only_NA]
not_NA_df <- dat[,not_NA]

#Null hypothesis of randomization test is that there is no phylogenetic signal,
#so low p-value indicates that there IS phylogenetic signal

run_phylosig <- function (phy, df) {
  sig_results <- vector("list", length = ncol(df))
  for(i in 1:length(sig_results)){
    res <- phylosig(phy, df[,i], method = "K", test = TRUE)
    sig_results[[i]] <- res
  }
  sig_results
}

only_NA_sig <- run_phylosig(cleaned$phy, only_NA_df)
not_NA_sig <- run_phylosig(cleaned$phy, not_NA_df)

get_K_and_P <- function (sig) {
  K <- sig$K
  P <- sig$P
  list("K" = K, "P-value" = P)
}

only_NA_vals <- map_df(only_NA_sig, get_K_and_P) %>% mutate(subset = "only_NA")
not_NA_vals <- map_df(not_NA_sig, get_K_and_P) %>% mutate(subset = "not_NA")

K_P_vals <- bind_rows(only_NA_vals, not_NA_vals) %>% pivot_longer(cols = c(K, `P-value`), names_to = "Statistic")

K_P_vals %>% ggplot(aes(x = subset, y = value, fill = subset)) + geom_violin() + 
  geom_boxplot() + theme_bw() + theme(legend.position = "none") + facet_wrap(~Statistic) +
  xlab("Subset") + ylab("Blomberg's K Statistic")
ggsave("phylosig.png")

K_P_vals %>% filter(Statistic == "K") %>% ggplot(aes(x = subset, y = value, fill = subset)) + geom_violin() + 
  geom_boxplot() + theme_bw() + theme(legend.position = "none") +
  xlab("Subset") + ylab("Blomberg's K Statistic")
ggsave("Supp_Fig_2.png")

#seeing differences between genes without and with phylogenetic signal
pvals %>% filter(is.na(s.hgt)) %>% pivot_longer(cols = c(-s.hgt)) %>% ggplot(aes(y = value, fill = name)) + geom_boxplot()
pvals %>% filter(!is.na(s.hgt)) %>% pivot_longer(cols = c(-s.hgt)) %>% ggplot(aes(y = value, fill = name)) + geom_boxplot()
