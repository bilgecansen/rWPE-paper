
#library(purrr)
library(foreach)
library(doSNOW)

source("fit_pe.R")



# Population dynamics -----------------------------------------------------

exp_grow <- function(mu, sd, xi, ts) {
  
  x <- c()
  x[1] <- xi
  
  r <- rnorm(ts, mu, sd)
  
  for (t in 1:(ts-1)) {
    
    x[t+1] <- x[t] + r[t] 
    
  }
  
  return(x)
}

rick_grow <- function(mu_r, sd, K, xi, ts, burnin) {
  
  r <- c()
  x <- c()
  x[1] <- log(xi)
  
  eps <- rnorm(ts, 0, sd)
  
  for (t in 1:(ts-1)) {
    
    r[t] <- mu_r*(1-(exp(x[t])/K))
    
    x[t+1] <- x[t] + r[t] + eps[t] 
    
  }
  
  return(x[-c(1:burnin)])
  
}

get_wpe <- function(n, sd = NA, r = 0, grow, sample = "white") {
  
  pb <- txtProgressBar(max = length(sd)*length(r), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
    
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  p_wpe <- foreach (h = 1:length(sd), .combine = "cbind") %:%
    foreach (i = 1:length(r), .combine = "c", .options.snow = opts, .export = c("rick_grow", "exp_grow")) %dopar% {
      
      p <- c()
      for (j in 1:100) {
      
        source("fit_pe.R")
        
        if (grow == "exp") {
          x1 <- exp_grow(mu = r[i], sd = sd[h], xi =  10, ts = n)
        }
        
        if (grow == "ricker") {
          x1 <- rick_grow(mu_r = r[i], sd = sd[h], K = 1000, xi =  100, ts = 1000, burnin = 1000 - n)
        }
        
        if (grow == "white") {
          x1 <- rnorm(n, log(100), sd[h])
        }
        
        if (sample == "white") {
          x2 <-  replicate(100, sample(x1)) 
        }
        
        if (sample == "random") {
          r <- x1[2:n] - x1[1:(n-1)]
          x2 <- exp_grow(mu = 0, sd = sd[h], xi =  x1[1], ts = n)
        }
        
        x1_wpe <- PE(x = exp(x1), weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1)
        null_wpe <- apply(exp(x2), 2, function(x) PE(x = x, weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1))
        null_diff <- x1_wpe - null_wpe 
        p[j] <- length(which(null_diff > 0)) / length(null_diff)
        
      }
      
      length(which(p < 0.05))
      
    }
  
  stopCluster(cl)
  
  rownames(p_wpe) <- r
  colnames(p_wpe) <- sd
  
  p_wpe
  
}

# Chaos
r_ch <- c(2.7, 2.8, 2.9, 3.0, 3.2, 3.3, 3.4, 3.5) 

wpe10_ch <- get_wpe(10, sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, grow = "ricker")
wpe20_ch <- get_wpe(20, sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, grow = "ricker")
wpe30_ch <- get_wpe(30, sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, grow = "ricker")
wpe40_ch <- get_wpe(40, sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, grow = "ricker")
wpe50_ch <- get_wpe(50, sd = c(0.001, 0.01, 0.1, 0.5), r = r_ch, grow = "ricker")

# Point cycles
r_pc <- c(2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6) 

wpe10_pc <- get_wpe(10, sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, grow = "ricker")
wpe20_pc <- get_wpe(20, sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, grow = "ricker")
wpe30_pc <- get_wpe(30, sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, grow = "ricker")
wpe40_pc <- get_wpe(40, sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, grow = "ricker")
wpe50_pc <- get_wpe(50, sd = c(0.001, 0.01, 0.1, 0.5), r = r_pc, grow = "ricker")

# Trends
r_tr <- seq(0, 1, 0.2) 

wpe10_tr <- get_wpe(10, sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, grow = "exp")
wpe20_tr <- get_wpe(20, sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, grow = "exp")
wpe30_tr <- get_wpe(30, sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, grow = "exp")
wpe40_tr <- get_wpe(40, sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, grow = "exp")
wpe50_tr <- get_wpe(50, sd = c(0.001, 0.01, 0.1, 0.5), r = r_tr, grow = "exp")

# Equilibrium dynamics
r_eq <- seq(0, 1.9, 0.3)

wpe10_eq <- get_wpe(10, sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, grow = "ricker")
wpe20_eq <- get_wpe(20, sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, grow = "ricker")
wpe30_eq <- get_wpe(30, sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, grow = "ricker")
wpe40_eq <- get_wpe(40, sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, grow = "ricker")
wpe50_eq <- get_wpe(50, sd = c(0.001, 0.01, 0.1, 0.5), r = r_eq, grow = "ricker")

# White noise
wpe10_wh <- get_wpe(10, sd = c(0.001, 0.01, 0.1, 0.5), grow = "white")
wpe20_wh <- get_wpe(20, sd = c(0.001, 0.01, 0.1, 0.5), grow = "white")
wpe30_wh <- get_wpe(30, sd = c(0.001, 0.01, 0.1, 0.5), grow = "white")
wpe40_wh <- get_wpe(40, sd = c(0.001, 0.01, 0.1, 0.5), grow = "white")
wpe50_wh <- get_wpe(50, sd = c(0.001, 0.01, 0.1, 0.5), grow = "white")


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


# Autoregressive dynamics -------------------------------------------------

ar_1 <- function(xi, beta, alpha, sd, ts, burnin) {
  
  x <- c()
  x[1] <- xi
  
  eps <- rnorm(ts, 0, sd)
  
  for (t in 1:(ts-1)) {
    x[t+1] = alpha + beta*x[t] + eps[t]
  }
  
  return(x[-c(1:burnin)])
  
}

get_wpe_ar <- function(n, beta, cv) {
  
  pb <- txtProgressBar(max = length(sd)*length(beta), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  p_wpe <- foreach (h = 1:length(cv), .combine = "cbind") %:%
    foreach (i = 1:length(beta), .combine = "c", .options.snow = opts, .export = "ar_1") %dopar% {
      
      source("fit_pe.R")
      
      p <- c()
      for (j in 1:100) {
        x1 <- ar_1(xi = 100, alpha = 50, beta = beta[i], sd = (10/(1-beta[i]))*cv[h], ts = 1000, burnin = 1000-n)
        
        x2 <-  replicate(100, sample(x1))
        
        x1_wpe <- PE(x = x1, weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1)
        null_wpe <- apply(x2, 2, function(x) PE(x = x, weighted = T,  word_length = 3, tau = 1, tie_method = "noise", noise_amount = 1))
        null_diff <- x1_wpe - null_wpe 
        p[j] <- length(which(null_diff > 0)) / length(null_diff)
        
      }
      
      length(which(p < 0.05))
      
    }
  
  stopCluster(cl)
  
  rownames(p_wpe) <- beta
  colnames(p_wpe) <- cv
  
  p_wpe
  
}

wpe_ar50 <- get_wpe_ar(n = 50, sd = c(0.1, 0.5, 1, 1.5, 2), beta = seq(0, 0.8, 0.2))

