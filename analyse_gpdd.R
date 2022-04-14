
library(tidyverse)
library(foreach)

gpdd_data <- read.csv("data.csv")

source("fit_pe.R")

# Filter GPDD data --------------------------------------------------------

# Time series length larger than 10
n <- table(gpdd_data$MainID)
idx_n <- which(n >= 10)

gpdd_data1 <- filter(gpdd_data, MainID %in% names(idx_n))

# No time skips
step_length <- function (x) {
  x[2:length(x)] - x[1:(length(x)-1)]
}

idx_step <-group_by(gpdd_data1, MainID) %>%
  summarise(steps = step_length(SampleYear)) %>%
  summarise(step_skip = any(steps > 1)) %>%
  ungroup() %>%
  filter(step_skip == FALSE) %>%
  select(MainID)

# No count replicated more than 15% 
gpdd_data2 <- filter(gpdd_data1, MainID %in% idx_step$MainID)

idx_rep <- group_by(gpdd_data2, MainID) %>%
  summarise(n = table(Population)/sum(table(Population))) %>%
  summarise(rep = any(n > 0.15)) %>%
  filter(rep == FALSE) %>%
  select(MainID)

gpdd_data3 <- filter(gpdd_data2, MainID %in% idx_rep$MainID)


# Calculate WPE for each time series --------------------------------------

ts <- unique(gpdd_data3$MainID)

wpe <- foreach(i = 1:length(ts), .combine = "c") %do% {
  
  x <- filter(gpdd_data3, MainID == ts[i])$Population
  
  PE(x = x, weighted = T,  word_length = 3, tau = 1, tie_method = "average")
  
}

ts_length <- group_by(gpdd_data3, MainID) %>%
  summarise(n = n())

plot(ts_length$n, wpe)


# Apply WPE randomization test to time series  ----------------------------

pb <- txtProgressBar(0, length(ts), style = 3)

p <- foreach(i = 1:length(ts), .combine = "c") %do% {
  
  x <- filter(gpdd_data3, MainID == ts[i])$Population
  
  x2 <- list()
  for (h in 1:1000) {
    x2[[h]] <- sample(x, length(x))
  }
  
  null_wpe <- map_dbl(x2, function(x) PE(x = x, weighted = T,  word_length = 3, tau = 1, tie_method = "average"))
  null_diff <- wpe[i] - null_wpe
  
  setTxtProgressBar(pb, i)
  
  length(which(null_diff < 0)) / length(null_diff)
  
}

plot(ts_length$n, p)
