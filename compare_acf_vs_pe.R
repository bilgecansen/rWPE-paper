
library(foreach)
library(doSNOW)
library(tidyverse)
library(patchwork)

source("functions_pe.R")
source("functions_pop.R")

r_ch <- c(2.4, 2.6, 2.85, 3.0, 3.14, 3.3)
r_tr <- c(0, 0.01, 0.1, 0.2, 0.4) 

x1a <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[1], sd = 0.1, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 990)
x1b <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[1], sd = 0.4, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 990)
x1c <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[1], sd = 0.1, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 960)
x1d <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[1], sd = 0.4, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 960)

x2a <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[6], sd = 0.1, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 990)
x2b <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[6], sd = 0.4, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 990)
x2c <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[6], sd = 0.1, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 960)
x2d <- foreach(i = 1:100) %do% 
  rick_grow(ri = r_ch[6], sd = 0.4, K = 1000, xi =  100, ts = 1000, burnin = 1000 - 960)

x3a <- foreach(i = 1:100) %do% 
  exp_grow(ri = r_tr[3], sd = 0.1, xi =  10, ts = 10)
x3b <- foreach(i = 1:100) %do% 
  exp_grow(ri = r_tr[3], sd = 0.4, xi =  10, ts = 10)
x3c <- foreach(i = 1:100) %do% 
  exp_grow(ri = r_tr[3], sd = 0.1, xi =  10, ts = 40)
x3d <- foreach(i = 1:100) %do% 
  exp_grow(ri = r_tr[3], sd = 0.4, xi =  10, ts = 40)

x4a <- foreach(i = 1:100) %do%
  rnorm(10,0,0.1)
x4b <- foreach(i = 1:100) %do%
  rnorm(10,0,0.4)
x4c <- foreach(i = 1:100) %do%
  rnorm(40,0,0.1)
x4d <- foreach(i = 1:100) %do%
  rnorm(40,0,0.4)

x_wpe1a <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x1a[[i]])[1]
x_acf1a <- foreach(i = 1:100, .combine = "c") %do% acf(x1a[[i]], plot = F)$acf[2]
x_wpe1b <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x1b[[i]])[1]
x_acf1b <- foreach(i = 1:100, .combine = "c") %do% acf(x1b[[i]], plot = F)$acf[2]
x_wpe1c <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x1c[[i]])[1]
x_acf1c <- foreach(i = 1:100, .combine = "c") %do% acf(x1c[[i]], plot = F)$acf[2]
x_wpe1d <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x1d[[i]])[1]
x_acf1d <- foreach(i = 1:100, .combine = "c") %do% acf(x1d[[i]], plot = F)$acf[2]

x_wpe2a <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x2a[[i]])[1]
x_acf2a <- foreach(i = 1:100, .combine = "c") %do% acf(x2a[[i]], plot = F)$acf[2]
x_wpe2b <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x2b[[i]])[1]
x_acf2b <- foreach(i = 1:100, .combine = "c") %do% acf(x2b[[i]], plot = F)$acf[2]
x_wpe2c <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x2c[[i]])[1]
x_acf2c <- foreach(i = 1:100, .combine = "c") %do% acf(x2c[[i]], plot = F)$acf[2]
x_wpe2d <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x2d[[i]])[1]
x_acf2d <- foreach(i = 1:100, .combine = "c") %do% acf(x2d[[i]], plot = F)$acf[2]

x_wpe3a <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x3a[[i]])[1]
x_acf3a <- foreach(i = 1:100, .combine = "c") %do% acf(x3a[[i]], plot = F)$acf[2]
x_wpe3b <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x3b[[i]])[1]
x_acf3b <- foreach(i = 1:100, .combine = "c") %do% acf(x3b[[i]], plot = F)$acf[2]
x_wpe3c <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x3c[[i]])[1]
x_acf3c <- foreach(i = 1:100, .combine = "c") %do% acf(x3c[[i]], plot = F)$acf[2]
x_wpe3d <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x3d[[i]])[1]
x_acf3d <- foreach(i = 1:100, .combine = "c") %do% acf(x3d[[i]], plot = F)$acf[2]

x_wpe4a <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x4a[[i]])[1]
x_acf4a <- foreach(i = 1:100, .combine = "c") %do% acf(x4a[[i]], plot = F)$acf[2]
x_wpe4b <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x4b[[i]])[1]
x_acf4b <- foreach(i = 1:100, .combine = "c") %do% acf(x4b[[i]], plot = F)$acf[2]
x_wpe4c <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x4c[[i]])[1]
x_acf4c <- foreach(i = 1:100, .combine = "c") %do% acf(x4c[[i]], plot = F)$acf[2]
x_wpe4d <- foreach(i = 1:100, .combine = "c") %do% calculate_pe(x4d[[i]])[1]
x_acf4d <- foreach(i = 1:100, .combine = "c") %do% acf(x4d[[i]], plot = F)$acf[2]

plot_comp <- function(x, y, title) {
  theme_set(theme_bw())
  ggplot() +
    geom_point(aes(x = x, y = y)) +
    labs(x = "Autocorrelation", y = "WPE", title = title)
}

g1a <- plot_comp(x_acf1a, x_wpe1a, title = "Periodic, sd = 0.1, t = 10")
g1b <- plot_comp(x_acf1b, x_wpe1b, title = "Periodic, sd = 0.4, t = 10")
g1c <- plot_comp(x_acf1c, x_wpe1c, title = "Periodic, sd = 0.1, t = 40")
g1d <- plot_comp(x_acf1d, x_wpe1d, title = "Periodic, sd = 0.4, t = 40")

(g1a | g1b) / (g1c | g1d)
ggsave("plot_acf_vs_wpe1.jpeg", width = 8, height = 8, units = "in")

g2a <- plot_comp(x_acf2a, x_wpe2a, title = "Chaotic, sd = 0.1, t = 10")
g2b <- plot_comp(x_acf2b, x_wpe2b, title = "Chaotic, sd = 0.4, t = 10")
g2c <- plot_comp(x_acf2c, x_wpe2c, title = "Chaotic, sd = 0.1, t = 40")
g2d <- plot_comp(x_acf2d, x_wpe2d, title = "Chaotic, sd = 0.4, t = 40")

(g2a | g2b) / (g2c | g2d)
ggsave("plot_acf_vs_wpe2.jpeg", width = 8, height = 8, units = "in")

g3a <- plot_comp(x_acf3a, x_wpe3a, title = "Linear, sd = 0.1, t = 10")
g3b <- plot_comp(x_acf3b, x_wpe3b, title = "Linear, sd = 0.4, t = 10")
g3c <- plot_comp(x_acf3c, x_wpe3c, title = "Linear, sd = 0.1, t = 40")
g3d <- plot_comp(x_acf3d, x_wpe3d, title = "Linear, sd = 0.4, t = 40")

(g3a | g3b) / (g3c | g3d)
ggsave("plot_acf_vs_wpe3.jpeg", width = 8, height = 8, units = "in")

g4a <- plot_comp(x_acf4a, x_wpe4a, title = "White Noise, sd = 0.1, t = 10")
g4b <- plot_comp(x_acf4b, x_wpe4b, title = "White Noise, sd = 0.4, t = 10")
g4c <- plot_comp(x_acf4c, x_wpe4c, title = "White Noise, sd = 0.1, t = 40")
g4d <- plot_comp(x_acf4d, x_wpe4d, title = "White Noise, sd = 0.4, t = 40")

(g4a | g4b) / (g4c | g4d)
ggsave("plot_acf_vs_wpe4.jpeg", width = 8, height = 8, units = "in")
