
#library(purrr)
library(foreach)
library(doSNOW)

source("fit_pe.R")


# Ricker model ------------------------------------------------------------

rick_grow <- function(mu_r, sd, K, xi, ts, burnin) {
  
  R <- c()
  x <- c()
  x[1] <- xi
  
  eps <- rnorm(ts, 0, sd)
  
  for (t in 1:ts) {
    
    R[t] <- exp(mu_r*(1-(x[t]/K)) + eps[t]) 
    
    x[t+1] <- x[t]*R[t]
    
  }
  
  return(x[-c(1:burnin)])
  
}

r <- c(0.1, 0.5, 1, 1.5, 2, 2.4, 3, 3.4)
sd <- c(0.001, 0.1, 0.25, 0.5)
ts <- c(10, 15, 20, 25, 30, 40, 50)

cl <- makeCluster(28, types = "SOCK")
registerDoSNOW(cl)

p_wpe <- foreach (i = 1:length(r)) %:%
  foreach (h = 1:length(sd)) %:%
  foreach (k = 1:length(ts)) %:%
  foreach (j = 1:100, .combine = "c", .packages = "foreach") %dopar% {
    
    x1 <- rick_grow(mu_r = r[i], sd = sd[h], K = 1000, xi =  100, ts = 1000, burnin = 1000 - ts[k])
    while (sd(x1) < 0.001) {
      x1 <- rick_grow(mu_r = r[i], sd = sd[h], K = 1000, xi =  100, ts = 1000, burnin = 1000 - ts[k])
    } 
    #plot(x1, type = "l")
    
    x2 <- list()
    for (l in 1:1000) {
      x2[[l]] <- sample(x1, length(x1))
    }
    
    x1_wpe <- PE(x = x1, weighted = T,  word_length = 3, tau = 1, tie_method = "average")
    null_wpe <- foreach(l = 1:length(x2), .combine = "c") %do% 
      PE(x = x2[[l]], weighted = T,  word_length = 3, tau = 1, tie_method = "average")
    null_diff <- x1_wpe - null_wpe 
    p <- length(which(null_diff < 0)) / length(null_diff)
    
    return(p)
    
  }
    
stopCluster(cl)

saveRDS(p_wpe, "p_wpe.RDS")
