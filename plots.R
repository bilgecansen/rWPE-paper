
library(tidyverse)
library(patchwork)

source("functions_pop.R")
source("functions_pe.R")

theme_set(theme_bw())

sd_sim <-  c(0.001, 0.01, 0.1, 0.2, 0.4)
r_tr <- c(0, 0.01, 0.1, 0.2, 0.4)
r_eq <- seq(0.1, 1.9, 0.3)
r_ch <- c(2.4, 2.6, 2.85, 3.0, 3.14, 3.3)

res <- readRDS("results_sim.rds")

plot_pe_tiles <- function(df, title) {
  ggplot() +
    geom_tile(data = df, mapping = aes(x = x, y = y, fill = z), color = "black", alpha = 0.8) +
    geom_text(data = df, mapping = aes(x = x, y = y, label = z), color = "white")+ 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(limits = c(0,1),
                        low = "darkred", 
                        high = "blue4") +
    labs(title = paste ("t = ", title, sep = ""), y = "r", x = bquote(sigma)) +
    theme(panel.border = element_rect(size = 2),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
}


# Weighted Permutation Entropy --------------------------------------------

ts_length <- c(10, 20, 30, 40)

df_tr <- foreach(i = 1:4) %do% {
  data.frame(
    x = rep(as.factor(sd_sim), length(r_tr)),
    y = rep(as.factor(r_tr), each = length(sd_sim)),
    z = c(t(res$wpe_tr[[i]]/1000))
  )
}

df_ch <- foreach(i = 1:4) %do% {
  data.frame(
    x = rep(as.factor(sd_sim), length(r_ch)),
    y = rep(as.factor(r_ch), each = length(sd_sim)),
    z = round(c(t(res$wpe_ch[[i]]/1000)),3)
  )
}

df_eq <- foreach(i = 1:4) %do% {
  data.frame(
    x = rep(as.factor(sd_sim), length(r_eq)),
    y = rep(as.factor(r_eq), each = length(sd_sim)),
    z = c(t(res$wpe_eq[[i]]/1000))
  )
}

df_eq_cc <- foreach(i = 1:4) %do% {
  data.frame(
    x = rep(as.factor(sd_sim), length(r_eq)),
    y = rep(as.factor(r_eq), each = length(sd_sim)),
    z = round(c(t(res$wpe_eq_cc[[i]]/1000)),3)
  )
}

df_wh <- foreach(i = 1:4) %do% {
  data.frame(
    x = as.factor(sd_sim),
    y = rep(0, length(sd_sim)),
    z = res$wpe_wh[[i]][1,]/1000
  )
}
 
theme_set(theme_bw())

g_tr <- foreach(i = 1:4) %do% {
  plot_pe_tiles(filter(df_tr[[i]], y !=0), title = ts_length[i])
}

g_ch <- foreach(i = 1:4) %do% {
  plot_pe_tiles(df_ch[[i]], title = ts_length[i])
}

g_eq <- foreach(i = 1:4) %do% {
  plot_pe_tiles(df_eq[[i]], title = ts_length[i])
}

g_eq_cc <- foreach(i = 1:4) %do% {
  plot_pe_tiles(df_eq_cc[[i]], title = ts_length[i])
}

g_wh <- foreach(i = 1:4) %do% {
  plot_pe_tiles(df_wh[[i]], title = ts_length[i]) + 
    theme(axis.title.y = element_blank())
}

(g_tr[[1]] + g_tr[[2]]) / (g_tr[[3]] + g_tr[[4]])
ggsave("plot_wpe_tr.jpeg", width = 8, height = 8, units = "in")

(g_ch[[1]] + g_ch[[2]]) / (g_ch[[3]] + g_ch[[4]])
ggsave("plot_wpe_ch.jpeg", width = 8, height = 8, units = "in")

(g_eq[[1]] + g_eq[[2]]) / (g_eq[[3]] + g_eq[[4]])
ggsave("plot_wpe_eq.jpeg", width = 8, height = 8, units = "in")

(g_eq_cc[[1]] + g_eq_cc[[2]]) / (g_eq_cc[[3]] + g_eq_cc[[4]])
ggsave("plot_wpe_eq_cc.jpeg", width = 8, height = 8, units = "in")

(g_wh[[1]] + g_wh[[2]]) / (g_wh[[3]] + g_wh[[4]])
ggsave("plot_wpe_wh.jpeg", width = 8, height = 8, units = "in")
