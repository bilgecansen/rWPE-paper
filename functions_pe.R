
calculate_pe <- function(x, m = 3, tau = 1) {
  
  if (m == 1 | m == 0) stop("m cannot be 1 or 0")
  if (tau == 0) stop("tau cannot be 0")
  
  eff_m <- (m-1)*(tau-1) + m
  
  if (eff_m >= length(x)) stop("Combination of m and tau is too large for the time series length")
  
  np <- length(x) - eff_m + 1
  
  idx_mat <- matrix(seq(1, eff_m, by = tau), ncol = m, nrow = np, byrow = T)
  idx_mat <- apply(idx_mat, 2, function(z) z + (0:(np-1)))
  
  ts_mat <- t(apply(idx_mat, 1, function(y) x[y]))
  idx_NA <- apply(ts_mat, 1, function(z) any(is.na(z)))
  ts_mat <- ts_mat[!idx_NA,]
  
  if (is.vector(ts_mat)) stop("Too many NAs in the time series realive to its length")
  if (nrow(ts_mat) == 0) stop("Too many NAs in the time series realive to its length")
  
  create_word <- function(y) {
    paste(rank(y, ties.method = "first"), collapse = "-")
  }
  
  wd <- apply(ts_mat, 1, create_word)
  
  wd_dis <- table(wd)
  
  p <- wd_dis/sum(wd_dis)
  
  pe <- -sum(p*log2(p))/log2(factorial(m))
  
  word_var <- apply(ts_mat, 1, function(x) var(x))
  word_df <- data.frame(words = wd, var = word_var)
  pw <- tapply(word_var, wd, sum)/sum(word_var)
  
  wpe <- -sum(pw*log2(pw))/log2(factorial(m))
  
  res <- c(wpe, pe, nrow(ts_mat))
  names(res) <- c("WPE", "PE", "NP")

  return(res)
}

random_grow <- function(x) {
  x_r <- x[2:length(x)] - x[1:(length(x)-1)]
  x_r <- sample(x_r)
  x2 <- c()
  x2[1] <- sample(x, 1)
  for (t in 1:(length(x)-1)) {
    x2[t+1] <- x2[t] + x_r[t]
  }
  
  return(x2)
}

test_wpe <- function(x, m = 3, tau = 1, n_random, perm = "white") {
  
  if (perm == "white") {
    x2 <-  replicate(n_random, sample(x))
  }
  
  if (perm == "random") {
    x2 <- replicate(n_random, random_grow(x))
  }
  
  y <- calculate_pe(x, m, tau)
  
  x_wpe <- y[1]
  x_np <- y[3]
  
  null_wpe <- apply(x2, 2, function(z) calculate_pe(z, m = m, tau = tau)[1])
  null_diff_wpe <- x_wpe - null_wpe 
  p_wpe <- length(which(null_diff_wpe > 0)) / length(null_diff_wpe)
  res <- c(x_wpe, x_np, p_wpe)
  names(res) <- c("WPE", "NP", "p-value(WPE)")
  
  return(res)
}

test_pe <- function(x, m = 3, tau = 1, n_random, perm = "white") {
  
  if (perm == "white") {
    x2 <-  replicate(n_random, sample(x))
  }
  
  if (perm == "random") {
    x2 <- replicate(n_random, random_grow(x))
  }
  
  y <- calculate_pe(x, m, tau)
  
  x_pe <- y[2]
  x_np <- y[3]
  
  null_pe <- apply(x2, 2, function(z) calculate_pe(z, m = m, tau = tau)[2])   
  null_diff_pe <- x_pe - null_pe 
  p_pe <- length(which(null_diff_pe > 0)) / length(null_diff_pe)
  res <- c(x_pe, x_np, p_pe)
  names(res) <- c("PE", "NP", "p-value(PE)")
  
  return(res)
}
