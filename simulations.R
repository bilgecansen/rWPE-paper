
library(foreach)
library(doSNOW)

source("functions_pe.R")
source("functions_pop.R")


# Permutations with white noise -------------------------------------------

ss <- c(10,20,30,40)

# Chaos
r_ch <- c(2.8, 2.9, 3.0, 3.1, 3.3, 3.4) 

wpe_ch <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, perm = "white")

pe_ch <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, perm = "white", type = "PE")

# Point cycles
r_pc <- c(2.4, 2.5, 2.6) 

wpe_pc <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, perm = "white")

pe_pc <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, perm = "white", type = "PE")

# Equilibrium dynamics
r_eq <- seq(0.1, 1.9, 0.3)

wpe_eq <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, perm = "white")

pe_eq <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, perm = "white", type = "PE")

# Trends
r_tr <- seq(0, 1, 0.2) 

wpe_tr <- foreach(k = 1:length(ss)) %do% 
  get_pe_exp(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, perm = "white")

pe_tr <- foreach(k = 1:length(ss)) %do% 
  get_pe_exp(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, perm = "white", type = "PE")

# White noise
wpe_wh <- foreach(k = 1:length(ss)) %do%  
  get_pe_white(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), perm = "white")

pe_wh <- foreach(k = 1:length(ss)) %do%  
  get_pe_white(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), perm = "white", type = "PE")

res <- list(wpe_ch = wpe_ch,
            pe_ch = pe_ch,
            wpe_pc = wpe_pc,
            pe_pc = pe_pc,
            wpe_eq = wpe_eq,
            pe_eq = pe_eq,
            wpe_tr = wpe_tr,
            pe_tr = pe_tr,
            wpe_wh = wpe_wh,
            pe_wh = pe_wh)

saveRDS(res, "results_sim.rds")
