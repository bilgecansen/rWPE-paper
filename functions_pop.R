
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

random_walk <- function(sd_r, xi, ts) {
  x <- c()
  x[1] <- log(xi)
  for (t in 1:(ts-1)) {
    x[t+1] <- x[t] + rnorm(1, 0, sd_r)
  }
  
  return(x)
}

get_pe_white <- function(n, sd = NA, r = 0, perm = "white", type = "WPE") {
  
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
      
      for (j in 1:100) {
        x1 <- rnorm(n, log(100), sd[h])
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < 0.05))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

get_pe_walk <- function(n, sd = NA, r = 0, perm = "white", type = "WPE") {
  
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
      
      for (j in 1:100) {
        x1 <- random_walk(sd_r = sd, xi = 100, ts = n)
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < 0.05))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

get_pe_rick <- function(n, sd = NA, r = 0, perm = "white", type = "WPE") {
  
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
      
      for (j in 1:100) {
        x1 <- rick_grow(mu_r = r[i], sd = sd[h], K = 1000, xi =  100, ts = 1000, burnin = 1000 - n)
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < 0.05))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}

get_pe_exp <- function(n, sd = NA, r = 0, perm = "white", type = "WPE") {
  
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
      
      for (j in 1:100) {
        x1 <- exp_grow(mu = r[i], sd = sd[h], xi =  10, ts = n)
        p[j] <- func(x = x1, n_random = 100, perm = perm)[3]
      }
      
      length(which(p < 0.05))
      
    }
  
  stopCluster(cl)
  
  rownames(p_pe) <- r
  colnames(p_pe) <- sd
  
  p_pe
}