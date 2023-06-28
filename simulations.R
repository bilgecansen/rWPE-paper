
library(foreach)
library(doSNOW)

source("functions_pe.R")
source("functions_pop.R")


# Permutations with white noise -------------------------------------------

ss <- c(10, 20, 30, 40)

# Chaos and point cycles
r_ch <- c(2.4, 2.6, 2.85, 3.0, 3.14, 3.3) 

wpe_ch <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.2, 0.4), r = r_ch, perm = "white")

# Equilibrium dynamics
r_eq <- seq(0.1, 1.9, 0.3)

wpe_eq <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.2, 0.4), r = r_eq, perm = "white")

# Equilibrium dynamics with trend
wpe_eq_cc <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick_cc(ss[k], sd = c(0.001, 0.01, 0.1, 0.2, 0.4), r = r_eq, perm = "white")

# Trends
r_tr <- c(0, 0.01, 0.1, 0.2, 0.4) 

wpe_tr <- foreach(k = 1:length(ss)) %do% 
  get_pe_exp(ss[k], sd = c(0.001, 0.01, 0.1, 0.2, 0.4), r = r_tr, perm = "white")

# White noise
wpe_wh <- foreach(k = 1:length(ss)) %do%  
  get_pe_white(ss[k], sd = c(0.001, 0.01, 0.1, 0.2, 0.4), perm = "white")

res <- list(wpe_ch = wpe_ch,
            wpe_eq = wpe_eq,
            wpe_eq_cc = wpe_eq_cc,
            wpe_tr = wpe_tr,
            wpe_wh = wpe_wh)

saveRDS(res, "results_sim.rds")


acf_tr <- foreach(k = 1:length(ss)) %do% 
  get_acf_exp(ss[1], sd = c(0.001, 0.01, 0.1, 0.2, 0.4), r = r_tr[1])
