
library(foreach)
library(doSNOW)

source("functions_pe.R")
source("fit_pe.R")

dat_pr_50 <- filter(dat_rogers, Classification == "periodic" & TimeSeriesLength == 50)
id_pr_50 <- dat_pr_50$ID %>% unique()

dat_ch_50 <- filter(dat_rogers, Classification == "chaotic" & TimeSeriesLength == 50)
id_ch_50 <- dat_ch_50$ID %>% unique()

compare_wpe <- function(data, id) {
  
  pb <- txtProgressBar(max = length(id), style = 3) 
  foreach(i = 1:length(id), .combine = "c") %do% {
    z <- filter(data, ID == id[i]) %>%
      select(contains("Sim.")) %>%
      as.matrix()
    
    z_res1 <- apply(z, 2, function(x) calculate_pe(x)[1])
    z_res2 <- apply(z, 2, function(x) PE(x = x, weighted = T,  word_length = 3, tau = 1, tie_method = "first"))
    
    
    setTxtProgressBar(pb, i)
    
    all.equal(z_res1, z_res2)
  }
}

comp_ch_50 <- compare_wpe(dat_ch_50, id_ch_50)
comp_pr_50 <- compare_wpe(dat_pr_50, id_pr_50)

