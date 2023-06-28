
exp_grow <- function(ri, sd, xi, ts) {
  
  x <- c()
  x[1] <- rnorm(1, log(xi), sd)
  
  r <- rnorm(ts, ri, sd)
  
  for (t in 1:(ts-1)) {
    
    x[t+1] <- x[t] + r[t] 
    
  }
  
  return(exp(x))
}

rick_grow <- function(ri, sd, K, xi, ts, burnin) {
  
  r <- c()
  x <- c()
  x[1] <- rnorm(1, log(xi), sd)
  
  eps <- rnorm(ts, 0, sd)
  
  for (t in 1:(ts-1)) {
    
    r[t] <- ri*(1-(exp(x[t])/K))
    
    x[t+1] <- x[t] + r[t] + eps[t]  
    
  }
  
  return(exp(x[-c(1:burnin)]))
  
}

rick_grow_cc <- function(ri, sd, K, xi, ts, burnin) {
  
  r <- c()
  x <- c()
  x[1] <- rnorm(1, log(xi), sd)
  
  eps <- rnorm(ts, 0, sd)
  
  for (t in 1:(ts-1)) {
    
    r[t] <- ri*(1-(exp(x[t])/K))
    
    x[t+1] <- x[t] + r[t] + eps[t]
    
    K <- K + 100
    
  }
  
  return(exp(x[-c(1:burnin)]))
  
}

random_walk <- function(sd, xi, ts) {
  x <- c()
  x[1] <- rnorm(1, log(xi), sd)
  for (t in 1:(ts-1)) {
    x[t+1] <- x[t] + rnorm(1, 0, sd)
  }
  
  return(exp(x))
}

get_pe_white <- function(n, sd = NA, r = 0, perm = "white", type = "WPE", alpha = 0.05) {
  
  pb <- txtProgressBar(max = length(sd)*length(r), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  if (type == "WPE") func <- test_wpe
  if (type == "PE") func <- test_pe
  
  p_pe <- 
    foreach (h = 1:length(sd), .combine = "cbind") %:%
    foreach (i = 1:length(r), .combine = "c", .options.snow = opts) %dopar% {
      
      source("functions_pe.R")
      
      p <- c()
      
      for (j in 1:1000) {
        x1 <- exp(rnorm(n, log(100), sd[h]))
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < alpha))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

get_pe_walk <- function(n, sd = NA, r = 0, perm = "white", type = "WPE", alpha = 0.05) {
  
  pb <- txtProgressBar(max = length(sd)*length(r), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  if (type == "WPE") func <- test_wpe
  if (type == "PE") func <- test_pe
  
  p_pe <- 
    foreach (h = 1:length(sd), .combine = "cbind") %:%
    foreach (i = 1:length(r), .combine = "c", .options.snow = opts, .export = c("random_walk")) %dopar% {
      
      source("functions_pe.R")
      
      p <- c()
      
      for (j in 1:1000) {
        x1 <- exp(random_walk(sd_r = sd, xi = 100, ts = n))
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < alpha))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

get_pe_rick <- function(n, sd = NA, r = 0, perm = "white", type = "WPE", alpha = 0.05) {
  
  pb <- txtProgressBar(max = length(sd)*length(r), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  if (type == "WPE") func <- test_wpe
  if (type == "PE") func <- test_pe
  
  p_pe <- 
    foreach (h = 1:length(sd), .combine = "cbind") %:%
    foreach (i = 1:length(r), .combine = "c", .options.snow = opts, .export = c("rick_grow")) %dopar% {
      
      source("functions_pe.R")
      
      p <- c()
      
      for (j in 1:1000) {
        x1 <- rick_grow(ri = r[i], sd = sd[h], K = 1000, xi =  100, ts = 1000, burnin = 1000 - n)
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < alpha))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

get_pe_rick_cc <- function(n, sd = NA, r = 0, perm = "white", type = "WPE", alpha = 0.05) {
  
  pb <- txtProgressBar(max = length(sd)*length(r), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  if (type == "WPE") func <- test_wpe
  if (type == "PE") func <- test_pe
  
  p_pe <- 
    foreach (h = 1:length(sd), .combine = "cbind") %:%
    foreach (i = 1:length(r), .combine = "c", .options.snow = opts, .export = c("rick_grow", "rick_grow_cc")) %dopar% {
      
      source("functions_pe.R")
      
      p <- c()
      
      for (j in 1:1000) {
        x1 <- rick_grow(ri = r[i], sd = sd[h], K = 1000, xi =  100, ts = 950, burnin = 949)
        x2 <- rick_grow_cc(ri = r[i], sd = sd[h], K = 1000, xi =  x1, ts = n + 1, burnin =  1)
        p[j] <- func(x = x2, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < alpha))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

get_pe_exp <- function(n, sd = NA, r = 0, perm = "white", type = "WPE", alpha = 0.05) {
  
  pb <- txtProgressBar(max = length(sd)*length(r), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(8, types = "SOCK")
  registerDoSNOW(cl)
  
  if (type == "WPE") func <- test_wpe
  if (type == "PE") func <- test_pe
  
  p_pe <- 
    foreach (h = 1:length(sd), .combine = "cbind") %:%
    foreach (i = 1:length(r), .combine = "c", .options.snow = opts, .export = c("exp_grow")) %dopar% {
      
      source("functions_pe.R")
      
      p <- c()
      
      for (j in 1:1000) {
        x1 <- exp_grow(ri = r[i], sd = sd[h], xi =  10, ts = n)
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < alpha))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

