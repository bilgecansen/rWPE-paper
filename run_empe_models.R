
library(foreach)
library(tidyverse)
library(mapppdr)
library(rstan)
library(MCMCvis)
library(loo)
library(grid)
library(patchwork)

# Install with devtools::install_github("hrbrmstr/ggalt")
library(ggalt)

source("functions_pe.R")
source("functions_pop.R")

empe_sites <- read.csv("empe_sites_new.csv")
data_pop <- readRDS("data_pop_empe.rds")

## Fast ice regions
reg_data_ice <- select(data_pop$sat, site_id, site_number, ice_reg) %>% 
  distinct() %>%
  mutate(reg_no = as.numeric(as.factor(ice_reg)))

#reg_data_ice[which(reg_data_ice$site_id %in% c("DOLL", "SMTH", "JASN")),]$ice_reg <- "Weddell Sea"
#reg_data_ice[which(reg_data_ice$site_id %in% c("DOLL", "SMTH", "JASN")),]$reg_no <- 7

reg_n_ice <- foreach(i = 1:8) %do%
  filter(reg_data_ice, reg_no == i)$site_number

## pack ice regions
reg_data_p_ice <- select(data_pop$sat, site_id, site_number, p_ice_reg) %>% 
  distinct() %>%
  mutate(reg_no = as.numeric(as.factor(p_ice_reg)))

reg_n_p_ice <- foreach(i = 1:5) %do%
  filter(reg_data_p_ice, reg_no == i)$site_number


# Population model to estimate abundance ----------------------------------

data_stan_base <- readRDS("data_stan_empe.rds")

data_stan_base$reg_n_ice1 <- reg_n_ice[[1]]
data_stan_base$reg_n_ice2 <- reg_n_ice[[2]]
data_stan_base$reg_n_ice3 <- reg_n_ice[[3]]
data_stan_base$reg_n_ice4 <- reg_n_ice[[4]]
data_stan_base$reg_n_ice5 <- reg_n_ice[[5]]
data_stan_base$reg_n_ice6 <- reg_n_ice[[6]]
data_stan_base$reg_n_ice7 <- reg_n_ice[[7]]
data_stan_base$reg_n_ice8 <- reg_n_ice[[8]]

data_stan_base$reg_n_p_ice1 <- reg_n_p_ice[[1]]
data_stan_base$reg_n_p_ice2 <- reg_n_p_ice[[2]]
data_stan_base$reg_n_p_ice3 <- reg_n_p_ice[[3]]
data_stan_base$reg_n_p_ice4 <- reg_n_p_ice[[4]]
data_stan_base$reg_n_p_ice5 <- reg_n_p_ice[[5]]

res_base <- stan(file = "model_base_site.stan",
                 data = data_stan_base,
                 chains = 4,
                 cores = 4,
                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                 init_r = 0.01)

saveRDS(res_base, "results_pop_empe.rds")

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

