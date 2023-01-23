
library(foreach)
library(doSNOW)

source("functions_pe.R")
source("functions_pop.R")


# Permutations with white noise -------------------------------------------

ss <- c(10,20,30,40,50)

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


# Permutations with random walk -------------------------------------------

# Chaos
wpe_ch_r <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, perm = "random")

pe_ch_r  <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, perm = "random", type = "PE")

# Point cycles
wpe_pc_r  <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, perm = "random")

pe_pc_r  <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, perm = "random", type = "PE")

# Equilibrium dynamics
wpe_eq_r  <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, perm = "random")

pe_eq_r  <- foreach(k = 1:length(ss)) %do% 
  get_pe_rick(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, perm = "random", type = "PE")

# Trends
wpe_tr_r  <- foreach(k = 1:length(ss)) %do% 
  get_pe_exp(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, perm = "random")

pe_tr_r  <- foreach(k = 1:length(ss)) %do% 
  get_pe_exp(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, perm = "random", type = "PE")

# Random walk
wpe_wh_r  <- foreach(k = 1:length(ss)) %do%  
  get_pe_walk(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), perm = "random")

pe_wh_r  <- foreach(k = 1:length(ss)) %do%  
  get_pe_walk(ss[k], sd = c(0.001, 0.01, 0.1, 0.5), perm = "random", type = "PE")


# Compare two WPEs --------------------------------------------------------

compare_wpe <- function(n, sd = NA, r1 = 0, r2 = 0, x1, x2) {
  
  pb <- txtProgressBar(max = length(sd), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  p_wpe <- foreach (h = 1:length(sd), .combine = "c", .options.snow = opts, .export = c("rick_grow", "exp_grow")) %dopar% {
      
      p <- c()
      for (j in 1:100) {
        
        source("fit_pe.R")
        
        if (x1 == "exp") {
          x1 <- exp_grow(mu = r1, sd = sd[h], xi =  10, ts = n)
        }
        
        if (x1 == "ricker") {
          x1 <- rick_grow(mu_r = r1, sd = sd[h], K = 1000, xi =  100, ts = 1000, burnin = 1000 - n)
        }
        
        if (x1 == "white") {
          x1 <- rnorm(n, log(100), sd[h])
        }
        
        if (x2 == "exp") {
          x2 <- exp_grow(mu = r2, sd = sd[h], xi =  10, ts = n)
        }
        
        if (x2 == "ricker") {
          x2 <- rick_grow(mu_r = r2, sd = sd[h], K = 1000, xi =  100, ts = 1000, burnin = 1000 - n)
        }
        
        if (x2 == "white") {
          x2 <- rnorm(n, log(100), sd[h])
        }
        
       
        x1_ran <-  replicate(100, sample(x1)) 
        x2_ran <-  replicate(100, sample(x2)) 
        
        x1_wpe <- PE(x = exp(x1), weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1)
        x2_wpe <- PE(x = exp(x2), weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1)
        
        diff <- x1_wpe - x2_wpe
        
        
        null_wpe1 <- apply(exp(x1_ran), 2, 
                           function(x) PE(x = x, weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1))
        null_wpe2 <- apply(exp(x2_ran), 2, 
                           function(x) PE(x = x, weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1))
        
        null_diff <- null_wpe1 - null_wpe2 
        
        p[j] <- (length(which(null_diff > abs(diff))) + length(which(null_diff < -abs(diff))))/ length(null_diff) 
        
      }
      
      length(which(p < 0.05))
      
    }
  
  stopCluster(cl)
  
  names(p_wpe) <- sd
  
  p_wpe
  
}

# Chaos vs Chaos
diff_wpe10 <- compare_wpe(10, r1 = 3.5, r2 = 3.5, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "ricker", x2 = "ricker")
diff_wpe50 <- compare_wpe(50, r1 = 3.5, r2 = 3.5, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "ricker", x2 = "ricker")

# Trend vs Trend
diff_wpe10 <- compare_wpe(10, r1 = 0.5, r2 = 0.5, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "exp", x2 = "exp")
diff_wpe50 <- compare_wpe(50, r1 = 0.5, r2 = 0.5, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "exp", x2 = "exp")

# Chaos vs White Noise
diff_wpe10 <- compare_wpe(10, r1 = 3.5, r2 = 0, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "ricker", x2 = "white")
diff_wpe50 <- compare_wpe(50, r1 = 3.5, r2 = 0, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "ricker", x2 = "white")

# Chaos vs trend
diff_wpe10 <- compare_wpe(10, r1 = 3.5, r2 = 0.5, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "ricker", x2 = "exp")
diff_wpe50 <- compare_wpe(50, r1 = 3.5, r2 = 0.5, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "ricker", x2 = "exp")

# Trend vs White Noise
diff_wpe10 <- compare_wpe(10, r1 = 0.5, r2 = 0, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "exp", x2 = "white")
diff_wpe50 <- compare_wpe(50, r1 = 0.5, r2 = 0, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "exp", x2 = "white")

# White noise vs white noise
diff_wpe10 <- compare_wpe(10, r1 = 0, r2 = 0, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "white", x2 = "white")
diff_wpe50 <- compare_wpe(50, r1 = 0, r2 = 0, sd = c(0.001, 0.01, 0.1, 0.5), x1 = "white", x2 = "white")

