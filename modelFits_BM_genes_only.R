#Model fits just for genes that don't violate BM

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

#Vector of strings with names of genes consistent with BM
bm_prot <- read_tsv(file = "Data/control_set_consistent_w_bm.txt") %>%
  pivot_longer(cols = !c(Type, Score), values_to = "Prot") %>%
  select(!name) %>% unique() %>% pull(Prot)

bm_only <- dat[, bm_prot]

std.error <- function ( x ) {
  sd(x) / sqrt(length(x))
}

#First do arbutus with BM model
arbutusBM <- vector(mode = "list", length = ncol(bm_only))
for(j in 1:ncol(bm_only)){
  tmp <- cleaned$phy
  todrop <- c()
  for(i in 1:nrow(bm_only)){
    if(is.na(bm_only[[i,j]])){
      todrop <- append(todrop, rownames(bm_only)[i])
    }
  }
  tmp <- drop.tip(tmp, todrop)
  trimmed <- treedata(tmp, bm_only[,j], sort = TRUE)
  SE = std.error(trimmed$data)
  fitBM <- fitContinuous(trimmed$phy, trimmed$data, SE, model = "BM")
  arbutusBM[[j]] <- arbutus(fitBM)
}

BM_pvals <- map_df(arbutusBM, pvalue_arbutus)

saveRDS(BM_pvals, "only_BM_genes_pvals")
BM_piv <- BM_pvals %>% pivot_longer(cols = everything(), names_to = "tstat")
BM_piv %>% ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
ggsave("arbutus_only_BM_genes.png")

#Now pick best-fit model
fitResults <- vector(mode = "list", length = ncol(bm_only))
for(j in 1:ncol(bm_only)){
  tmp <- cleaned$phy
  todrop <- c()
  for(i in 1:nrow(bm_only)){
    if(is.na(bm_only[[i,j]])){
      todrop <- append(todrop, rownames(bm_only)[i])
    }
  }
  tmp <- drop.tip(tmp, todrop)
  trimmed <- treedata(tmp, bm_only[,j], sort = TRUE)
  SE = std.error(trimmed$data)
  fitBM <- fitContinuous(trimmed$phy, trimmed$data, SE, model = "BM")
  fitOU <- fitContinuous(trimmed$phy, trimmed$data, SE, model = "OU")
  fitEB <- fitContinuous(trimmed$phy, trimmed$data, SE, model = "EB")
  aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
  fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                       list(c(fitEB, model = "EB"))))
  fitResults[j] <- fit
}

saveRDS(fitResults, file = "only_bm_genes_allfit")

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, eb <- eb + 1))
  }
  data.frame(OU = ou, BM = bm, EB = eb)
}

df <- model_count(fitResults)

b <- df %>% pivot_longer(c(OU, BM, EB), names_to = "model")

b %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
ggsave("only_BM_allfit_AIC.png")

run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(f in fits){
    class(f) <- "gfit"
    arby[[count]] <- arbutus(f)
    count = count + 1
  }
  arby
}

arb_result <- run_arb(fitResults)
pvals <- map_df(arb_result, pvalue_arbutus)
saveRDS(pvals, file = "only_BM_allfit_pvalues_df")

p_piv <- pvals %>% pivot_longer(cols = everything(), names_to = "tstat")
saveRDS(p_piv, file = "only_BM_allfit_pvalues_table")

p_piv %>% ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
ggsave("only_BM_allfit_arbutus_results.png")

#Now checking the phylogenetic signal of one set of genes and comparing to other
num <- 1:dim(pvals)[1]
only_NA <- pvals %>% mutate(num = num) %>% filter(is.na(s.hgt)) %>% pull(num)
not_NA <- pvals %>% mutate(num = num) %>% filter(!is.na(s.hgt)) %>% pull(num)

only_NA_df <- bm_only[,only_NA]
not_NA_df <- bm_only[,not_NA]

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
  geom_boxplot() + facet_wrap(~Statistic) + theme_bw() + ggtitle("Phylogenetic signal of BM-only genes")+
  xlab("Subset of data")
ggsave("phylosig.png")
