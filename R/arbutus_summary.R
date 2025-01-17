#Summary stats for Arbutus results
library(tidyverse)

p_df <- readRDS("pvalues_df")

#remove msig 
p_analyze <- p_df %>% select(!m.sig) 

#Find which are less than 0.05
p_analyze2 <- p_analyze %>% transmute(c.less = c.var <= 0.05, sv.less = s.var <= 0.05, sa.less = s.asr <= 0.05, sh.less = s.hgt <= 0.05, d.less = d.cdf <= 0.05)

#Sum to get those that are less than 0.05
p_analyze3 <- p_analyze2 %>% transmute(inade = c.less + sv.less + sa.less + sh.less + d.less)

#Count
arb_summary <- p_analyze3 %>% count(inade) %>% mutate(prop = n/sum(n)) %>% mutate(inade = as.character(inade))
print(arb_summary)

arb_summary %>% ggplot(aes(x = inade, y = n, fill = inade)) + geom_bar(stat = "identity") + 
  xlab("Number of inadequacies") + ylab("Number of genes") + ggtitle("Amount of genes by number of inadequacies") + theme_bw()

#Same for author-chosen genes
chosen_df <- readRDS("only_BM_allfit_pvalues_df")

arb_summary_chosenonly <- chosen_df %>% select(!m.sig) %>% 
  transmute(c.less = c.var <= 0.05, sv.less = s.var <= 0.05, sa.less = s.asr <= 0.05, sh.less = s.hgt <= 0.05, d.less = d.cdf <= 0.05) %>% 
  transmute(inade = c.less + sv.less + sa.less + sh.less + d.less) %>% count(inade) %>% mutate(prop = n/sum(n)) %>% mutate(inade = as.character(inade))

arb_summary_chosenonly %>% ggplot(aes(x = inade, y = n, fill = inade)) + geom_bar(stat = "identity") + 
  xlab("Number of inadequacies") + ylab("Number of genes") + ggtitle("Amount of genes by number of inadequacies") + theme_bw()

#Now for if we don't do best-fit
BM_df <- readRDS("BM_pvals")

arb_summary_BMonly <- BM_df %>% select(!m.sig) %>% 
  transmute(c.less = c.var <= 0.05, sv.less = s.var <= 0.05, sa.less = s.asr <= 0.05, sh.less = s.hgt <= 0.05, d.less = d.cdf <= 0.05) %>% 
  transmute(inade = c.less + sv.less + sa.less + sh.less + d.less) %>% count(inade) %>% mutate(prop = n/sum(n)) %>% mutate(inade = as.character(inade))

arb_summary_BMonly %>% ggplot(aes(x = inade, y = n, fill = inade)) + geom_bar(stat = "identity") + 
  xlab("Number of inadequacies") + ylab("Number of genes") + ggtitle("Amount of genes by number of inadequacies") + theme_bw()

#Now not best-fit AND only genes authors looked at
BM_chosen_df <- readRDS("only_BM_genes_pvals")

arb_summary_BMchosenonly <- BM_chosen_df %>% select(!m.sig) %>% 
  transmute(c.less = c.var <= 0.05, sv.less = s.var <= 0.05, sa.less = s.asr <= 0.05, sh.less = s.hgt <= 0.05, d.less = d.cdf <= 0.05) %>% 
  transmute(inade = c.less + sv.less + sa.less + sh.less + d.less) %>% count(inade) %>% mutate(prop = n/sum(n)) %>% mutate(inade = as.character(inade))

arb_summary_BMchosenonly %>% ggplot(aes(x = inade, y = n, fill = inade)) + geom_bar(stat = "identity") + 
  xlab("Number of inadequacies") + ylab("Number of genes") + ggtitle("Amount of genes by number of inadequacies when fitting a BM model") + theme_bw()


#Compare Pvals from just BM to all fit
BM_only <- readRDS("BM_pvals") %>% mutate(models = "BM") %>% select(!m.sig) %>% pivot_longer(cols = -c(models), names_to = "TestStatistic", values_to = "P-values")
all_models <- readRDS("pvalues_df") %>% mutate(models = "Best Fit") %>% select(!m.sig) %>% pivot_longer(cols = -c(models), names_to = "TestStatistic", values_to = "P-values")

all_pvals_data <- full_join(BM_only, all_models)

all_pvals_data %>% ggplot(aes(x = `P-values`, fill = TestStatistic)) + geom_histogram() + 
  facet_grid(rows = vars(models), cols = vars(TestStatistic)) + theme_bw() + geom_vline(xintercept = 0.05) +
  theme(legend.position = "none", panel.grid = element_blank()) + scale_fill_brewer(palette ="Set1") 
