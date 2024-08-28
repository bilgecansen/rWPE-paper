
library(tidyverse)
library(MCMCvis)
library(ggrepel)
library(foreach)
library(patchwork)

source("functions_pe.R")

#####################################
# ADELIE PENGUINs
#####################################

data_pop <- readRDS("data_penguin/data_pop_adpe.rds")
data_stan_null <- readRDS("data_penguin/data_stan_null_adpe.rds")

res_null <- readRDS("results/results_pop_adpe.rds")

# Calculate long-term ice area average ------------------------------------

sites <- data_pop$site_list$site_id

ice <- readRDS("data_penguin/data_forced_finn_500km_adpe.rds") %>%
  filter(season > 1978) %>%
  filter(site_id %in% sites) %>%
  select(site_id, season, aice) %>%
  group_by(site_id) %>%
  summarise(aice_avg = mean(aice))

aice_std <- (ice$aice_avg - mean(ice$aice_avg))/sd(ice$aice_avg)

# Extract demographic parameters ------------------------------------------

n_sites <- data_stan_null$n_sites
n_seasons <- data_stan_null$n_seasons
sites <- data_pop$site_list$site_id

params_N <- 
  foreach(i = 1:n_sites) %:% 
  foreach(h = 1:n_seasons, .combine = "c") %do%
  paste("lz", "[", i, ",", h, "]", sep = "")

N_all <- foreach(i = 1:n_sites) %do% {
  MCMCsummary(res_null, params = params_N[[i]], ISB = F, exact = T)[,1]
}

season_relative <- foreach(i = 1:n_sites) %do% {
  min_year <- data_pop$abundance_initial$season_relative[i]
  max_year <- max(filter(data_pop$abundance_nests, site_id == sites[i])$season_relative)
  
  c(min_year, max_year)
}

names(season_relative) <- data_pop$site_list$site_id

seasons <-  foreach(i = 1:n_sites) %do% {
  min_year <- data_pop$abundance_initial$season[i]
  max_year <- max(filter(data_pop$abundance_nests, site_id == sites[i])$season)
  
  c(min_year, max_year)
}

names(seasons) <- data_pop$site_list$site_id

# subset model results with years with data
N_all <- foreach(i = 1:n_sites) %do% {
  N_all[[i]][season_relative[[i]][1]:(season_relative[[i]][2])]
}

# Replace abundance estimates with NA for years with no data
for(i in 1:n_sites) {
  x <- filter(data_pop$abundance_nests, site_id == sites[i])
  N_all[[i]][-c(1,x$season_relative)] <- NA
}

names(N_all) <- data_pop$site_list$site_id


# Calculate WPE of sites --------------------------------------------------

wpe_adpe <- foreach(i = 1:length(sites)) %do% {
  tryCatch({test_wpe(N_all[[i]], n_random = 1000)}, error=function(e){NA})
}
names(wpe_adpe) <- sites

idx_np <- which(map_dbl(wpe_adpe, function(x) x[2]) > 7)
wpe_adpe2 <- map_dbl(wpe_adpe, function(x) x[1])[idx_np]
wpe_sig <- map_dbl(wpe_adpe, function(x) x[3])[idx_np]
wpe_sig <- as.factor(as.numeric(wpe_sig < 0.05))

mu_r <- MCMCsummary(res_null, params = "mu_site")[,1]
beta <- MCMCchains(res_null, params = c("beta", "beta0"))
ice_sim <- seq(min(aice_std), max(aice_std), length.out = 100)
beta_idx <- sample(1:nrow(beta), 100)

z <- foreach(i = 1:length(beta_idx), .combine = "rbind") %do% {
  beta[beta_idx[i],3] + beta[beta_idx[i],1]*ice_sim + beta[beta_idx[i],2]*ice_sim^2
}

z2 <- pivot_longer(as.data.frame(z), cols = everything()) %>%
  add_column(iter = rep(1:100, each = 100),
             x = rep(ice_sim*sd(ice$aice_avg) + mean(ice$aice_avg), 100))

idx_wpe <- order(wpe_adpe2, decreasing = T)
sites_wpe <- factor(sites[idx_np][idx_wpe], levels = (sites[idx_np][idx_wpe]))

g_base <- ggplot() +
  geom_point(aes(x = factor(sites_wpe, levels = sites_wpe[15:1]), 
                 y = wpe_adpe2[idx_wpe], col = wpe_sig[idx_wpe]), size = 4) +
  geom_segment(aes(x = factor(sites_wpe, levels = sites_wpe[15:1]), 
                   xend = factor(sites_wpe, levels = sites_wpe[15:1]), 
                   y = 0, yend = wpe_adpe2[idx_wpe], col = wpe_sig[idx_wpe]),
               size = 1.1) +
  scale_color_manual(name = NULL,
                    labels = c("0" = "p > 0.05",
                               "1" = "p < 0.05"),
                    values = c("0" = "darkred", 
                               "1" = "darkblue")) +
  labs(y = "WPE", title = "Adélie Penguins") +
  theme(legend.position = "none",
        title = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  coord_flip()

g_quad <- ggplot() +
  geom_line(mapping = aes(x = z2$x,
                          y = z2$value,
                          group = z2$iter), alpha = 0.2) +
  geom_point(mapping = aes(x = ice$aice_avg[-idx_wpe], 
                           y = mu_r[-idx_wpe]), size = 3, alpha = 0.9, col = "grey") +
  geom_point(mapping = aes(x = ice$aice_avg[-idx_wpe], 
                           y = mu_r[-idx_wpe]), size = 3, shape =  1, col = "black") +
  geom_point(mapping = aes(x = ice$aice_avg[idx_wpe], 
                           y = mu_r[idx_wpe], col = wpe_sig), size = 4, alpha = 0.8) +
  geom_point(mapping = aes(x = ice$aice_avg[idx_wpe], 
                           y = mu_r[idx_wpe]), size = 4, shape = 1, col = "black") +
  #scale_y_continuous(limits = c(-0.07, 0.05)) +
  scale_color_manual(name = NULL,
                     labels = c("0" = "p > 0.05",
                                "1" = "p < 0.05"),
                     values = c("0" = "darkred", 
                                "1" = "darkblue")) +
  labs(x = "Ice Concentration", y = "Average Growth") +
  theme(legend.position = "none",
        title = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank())


#####################################
# EMPEROR PENGUINs
#####################################

res_base <- readRDS("results/results_pop_empe.rds")

empe_sites <- read.csv("data_penguin/empe_sites_new.csv")
data_pop <- readRDS("data_penguin/data_pop_empe.rds")

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
wpe_sig_empe <- as.factor(as.numeric(wpe_empe_col3[idx_empe] < 0.05))
empe_site_id <- data_pop$sat$site_id %>% unique()
empe_site_id <- empe_site_id[order(empe_site_id)]
empe_site_id <- factor(empe_site_id[idx_empe], levels = empe_site_id[idx_empe])


# Plots -------------------------------------------------------------------

theme_set(theme_bw())

g_col <- ggplot() +
  geom_point(aes(x = factor(empe_site_id, levels = empe_site_id[41:1]), 
                 y = wpe_empe_col2[idx_empe], col = wpe_sig_empe), size = 4) +
  geom_segment(aes(x = factor(empe_site_id, levels = empe_site_id[41:1]), 
                   xend = factor(empe_site_id, levels = empe_site_id[41:1]), 
                   y = 0, yend = wpe_empe_col2[idx_empe], col = wpe_sig_empe), 
               size = 1.1) +
  scale_color_manual(name = NULL,
                    labels = c("0" = "p > 0.05",
                               "1" = "p < 0.05"),
                    values = c("0" = "darkred", 
                               "1" = "blue4")) +
  labs(y = "WPE", title = "Emperor Penguins") +
  theme(legend.position = "none",
        title = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) + 
  coord_flip()


#####################################
# Final Plot
#####################################

((g_base / g_quad) | g_col) + 
  plot_annotation(tag_levels = "a")

ggsave("figures/plot_penguin.pdf", width = 180, height = 180, units = "mm") 

# to extract the legend
g_legend <- ggplot() +
  geom_line(mapping = aes(x = z2$x,
                          y = z2$value,
                          group = z2$iter), alpha = 0.2) +
  geom_point(mapping = aes(x = ice$aice_avg[-idx_wpe], 
                           y = mu_r[-idx_wpe]), size = 4, alpha = 0.9, col = "grey") +
  geom_point(mapping = aes(x = ice$aice_avg[-idx_wpe], 
                           y = mu_r[-idx_wpe]), size = 4, shape =  1, col = "black") +
  geom_point(mapping = aes(x = ice$aice_avg[idx_wpe], 
                           y = mu_r[idx_wpe], col = wpe_sig), size = 6, alpha = 0.8) +
  geom_point(mapping = aes(x = ice$aice_avg[idx_wpe], 
                           y = mu_r[idx_wpe]), size = 6, shape = 1, col = "black") +
  #scale_y_continuous(limits = c(-0.07, 0.05)) +
  scale_color_manual(name = NULL,
                     labels = c("0" = "p > 0.05",
                                "1" = "p < 0.05"),
                     values = c("0" = "darkred", 
                                "1" = "darkblue")) +
  labs(x = "Ice Concentration", y = "Average Growth") +
  theme(legend.position = "bottom",
        title = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid.minor = element_blank())

g_legend
ggsave("plot_penguin_legend.pdf", width = 8, height = 6, units = "in") 
