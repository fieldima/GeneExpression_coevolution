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
arb_summary <- p_analyze3 %>% count(inade) %>% mutate(prop = n/sum(n))
print(arb_summary)