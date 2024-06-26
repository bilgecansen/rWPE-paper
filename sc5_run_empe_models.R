
library(foreach)
library(tidyverse)
library(rstan)
library(MCMCvis)
library(loo)
library(grid)
library(patchwork)

# Install with devtools::install_github("hrbrmstr/ggalt")
library(ggalt)

source("functions_pe.R")
source("functions_pop.R")

empe_sites <- read.csv("data_penguin/empe_sites_new.csv")
data_pop <- readRDS("data_penguin/data_pop_empe.rds")


# Population model to estimate abundance ----------------------------------

data_stan_base <- readRDS("data_penguin/data_stan_empe.rds")

res_base <- stan(file = "model_base_site.stan",
                 data = data_stan_base,
                 chains = 4,
                 cores = 4,
                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                 init_r = 0.01)

saveRDS(res_base, "results/results_pop_empe.rds")

# Colony level abundances
params_N <- 
  foreach(i = 1:41) %:% 
  foreach(h = 1:10, .combine = "c") %do%
  paste("N", "[", i, ",", h, "]", sep = "")

N_all <- foreach(i = 1:41) %do% {
  MCMCsummary(res_base, params = params_N[[i]], ISB = F, exact = T)[,1]
}

wpe_empe_col <- foreach(i = 1:41) %do% {
  test_wpe(log(N_all[[i]]), n_random = 1000)
}

wpe_empe_col2 <- foreach(i = 1:41, .combine = "c") %do% {
  wpe_empe_col[[i]][1]
}

wpe_empe_col3 <- foreach(i = 1:41, .combine = "c") %do% {
  wpe_empe_col[[i]][3]
}

idx_empe <- order(wpe_empe_col2, decreasing = T)
wpe_sig <- as.factor(as.numeric(wpe_empe_col3[idx_empe] < 0.05))
empe_site_id <- data_pop$sat$site_id %>% unique()
empe_site_id <- empe_site_id[order(empe_site_id)]
empe_site_id <- factor(empe_site_id[idx_empe], levels = empe_site_id[idx_empe])


# Plots -------------------------------------------------------------------

theme_set(theme_bw())

g_col <- ggplot() +
  geom_col(aes(x = empe_site_id, y = wpe_empe_col2[idx_empe], fill = wpe_sig)) +
  scale_fill_manual(name = NULL,
                    labels = c("0" = "p > 0.05",
                               "1" = "p < 0.05"),
                    values = c("0" = "darkred", 
                               "1" = "blue4")) +
  labs(y = "WPE") +
  theme(legend.position = "bottom",
        title = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8, angle = 90),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

g_col
ggsave("plot_empe.jpeg", width = 10, height = 6, units = "in")  

