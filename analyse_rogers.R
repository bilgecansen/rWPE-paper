
library(foreach)
library(doSNOW)
library(tidyverse)
library(patchwork)

source("functions_pe.R")
source("functions_pop.R")

dat_rogers <- read.csv("rogers_simulation_dataset_test.csv")

dat_pr_25 <- filter(dat_rogers, Classification == "periodic" & TimeSeriesLength == 25)
id_pr_25 <- dat_pr_25$ID %>% unique()

dat_pr_50 <- filter(dat_rogers, Classification == "periodic" & TimeSeriesLength == 50)
id_pr_50 <- dat_pr_50$ID %>% unique()

dat_ch_25 <- filter(dat_rogers, Classification == "chaotic" & TimeSeriesLength == 25)
id_ch_25 <- dat_ch_25$ID %>% unique()

dat_ch_50 <- filter(dat_rogers, Classification == "chaotic" & TimeSeriesLength == 50)
id_ch_50 <- dat_ch_50$ID %>% unique()

apply_wpe_test <- function(data, id) {
  
  pb <- txtProgressBar(max = length(id), style = 3) 
  foreach(i = 1:length(id)) %do% {
    z <- filter(data, ID == id[i]) %>%
      select(contains("Sim.")) %>%
      as.matrix()
    
    z_res <- apply(z, 2, function(x) test_wpe(x, n_random = 100)[3])
    
    setTxtProgressBar(pb, i)
    
    return(z_res)
  }
}

res_pr_25 <- apply_wpe_test(dat_pr_25, id_pr_25)
res_pr_50 <- apply_wpe_test(dat_pr_50, id_pr_50)
res_ch_25 <- apply_wpe_test(dat_ch_25, id_ch_25)
res_ch_50 <- apply_wpe_test(dat_ch_50, id_ch_50)

sig_pr_25 <- map_dbl(res_pr_25, function(x) which(x < 0.05) %>% length()/100)
sig_pr_50 <- map_dbl(res_pr_50, function(x) which(x < 0.05) %>% length()/100)
sig_ch_25 <- map_dbl(res_ch_25, function(x) which(x < 0.05) %>% length()/100)
sig_ch_50 <- map_dbl(res_ch_50, function(x) which(x < 0.05) %>% length()/100)

noise_level <- unique(dat_pr_25$NoiseLevel)
model_pr <- unique(dat_pr_25$Model)
model_ch <- unique(dat_ch_25$Model)

df_pr_25 <- data.frame(x = rep(as.factor(noise_level), times = length(model_pr)),
                       y = rep(model_pr, each = length(noise_level)),
                       z = sig_pr_25)

df_pr_50 <- data.frame(x = rep(as.factor(noise_level), times = length(model_pr)),
                       y = rep(model_pr, each = length(noise_level)),
                       z = sig_pr_50)

df_pr <- bind_rows(df_pr_25, df_pr_50) %>%
  add_column(ss = c(rep("25", 28), rep("50", 28)),
             model = c(rep("Periodic", 28), rep("Periodic", 28)))

df_ch_25 <- data.frame(x = rep(as.factor(noise_level), times = length(model_ch)),
                       y = rep(model_ch, each = length(noise_level)),
                       z = sig_ch_25)

df_ch_50 <- data.frame(x = rep(as.factor(noise_level), times = length(model_ch)),
                       y = rep(model_ch, each = length(noise_level)),
                       z = sig_ch_50)

df_ch <- bind_rows(df_ch_25, df_ch_50) %>%
  add_column(ss = c(rep("25", 28), rep("50", 28)),
             model = c(rep("Chaotic", 28), rep("Chaotic", 28)))

# Plots -------------------------------------------------------------------


plot_pe_tiles <- function(df, title) {
  
  theme_set(theme_bw())
  ggplot() +
    geom_tile(data = df, mapping = aes(x = x, y = y, fill = z), color = "black", alpha = 0.8) +
    geom_text(data = df, mapping = aes(x = x, y = y, label = z), color = "white") +
    facet_wrap(~ model + ss) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(limits = c(0,1),
                        low = "darkred", 
                        high = "blue4") +
    labs(x = "CV") +
    theme(panel.border = element_rect(size = 2),
          axis.title = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
}

g1 <- plot_pe_tiles(df_pr)
g2 <- plot_pe_tiles(df_ch)


(g1) / (g2)
ggsave("plot_rogers.jpeg", width = 8, height = 8, units = "in")


