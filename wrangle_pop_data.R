
# Wrangle MAPPPDR and environmental data for population models  

library(tidyverse)
#devtools::install_github("CCheCastaldo/mapppdr")
library(mapppdr)
library(foreach)
library(doSNOW)
library(abind)

penguin_obs <- filter(mapppdr::penguin_obs, species_id == "ADPE")
idx <- which(penguin_obs$site_id == "HUKU" & penguin_obs$type == "nests")
penguin_obs <- penguin_obs[-idx,]

# Count Data --------------------------------------------------------------

prep_pop_data <- function(species, min_season, max_season, sites, add = F) {
  
  sites <- penguin_obs %>%
    filter(count > 0 & species_id == species & (type == "nests" | type == "adults") & site_id %in% sites) %>%
    group_by(site_id) %>%
    summarize(ts_length = n())
  
  penguin_obs <- penguin_obs %>%
    # keep all sites that have at least 1 count between min and max season
    dplyr::filter(count > 0 & species_id == species & site_id %in% sites$site_id)
  
  site_list <- penguin_obs %>%
    dplyr::filter(season >= min_season & season <= max_season) %>%
    # create relative season index 
    mutate(season_relative = season - min_season + 1) %>%
    # determine first season a count is observed for each site
    group_by(site_id) %>%
    summarise(initial_season = min(season_relative), last_season = max(season_relative)) %>%
    ungroup() %>%
    # join to get other site specific covariates for visualization purposes
    left_join(mapppdr::sites, by = "site_id") %>%
    # create site index for model and visualization
    #filter(ccamlr_id == "88.1") %>% 
    mutate(site = as.numeric(as.factor(site_id))) %>%
    dplyr::select(site_id, site_name, ccamlr_id, ccamlr_id, site, initial_season, last_season, latitude, longitude)
  
  n_sites <- nrow(site_list)
  n_seasons <- max_season - min_season + 1
  seasons <- min_season:max_season
  
  abundance <- penguin_obs %>%
    dplyr::filter(season >= min_season & season <= max_season) %>%
    # join to get site index and initial season
    right_join(site_list, by = "site_id") %>%
    # create relative season index 
    mutate(season_relative = season - min_season + 1) %>%
    # ASSUMPTION: increase accuracy category of all adult counts by + 3 with a max error of 5
    rowwise() %>%
    mutate(accuracy = replace(accuracy, type == "adults", base::min((accuracy[type == "adults"] + 3), 5))) %>%
    ungroup() %>%  
    mutate(type = replace(type, type == "adults", "nests")) %>%
    # ASSUMPTION: keep maximum nest and chick count reported each season for a site
    group_by(site_id, season, season_relative, type) %>%
    arrange(desc(count), accuracy) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    # ASSUMPTION: convert accuracy to the following errors/precisions
    mutate(sigma = case_when(
      accuracy == 1 ~ 0.025, 
      accuracy == 2 ~ 0.051,
      accuracy == 3 ~ 0.122, 
      accuracy == 4 ~ 0.226, 
      accuracy == 5 ~ 0.821)) %>%
    mutate(precision = case_when(
      accuracy == 1 ~ 1/0.025^2, 
      accuracy == 2 ~ 1/0.051^2,
      accuracy == 3 ~ 1/0.122^2, 
      accuracy == 4 ~ 1/0.226^2, 
      accuracy == 5 ~ 1/0.821^2)) %>%  
    dplyr::select(site_id, site, ccamlr_id, season, season_relative, initial_season, last_season, type, 
                  count, accuracy, sigma, precision) %>%
    arrange(site, season_relative, type, -count, accuracy, sigma, precision)
  
  if (add) {
    
    for (i in 1:5) {
      z <- filter(abundance, site_id == site_list$site_id[i])
      idx <- which(!seasons %in% z$season)
      
      if (length(idx) > 0) {
        
        for (h in 1:length(idx)) {
          
          abundance <- abundance %>%
            add_row(site_id = z$site_id[1],
                    site = z$site[1],
                    ccamlr_id = z$ccamlr_id[1],
                    season = seasons[idx[h]],
                    season_relative = idx[h],
                    initial_season = z$initial_season[1],
                    count = round(mean(c(z$count[idx[h]+1], z$count[idx[h]-1])),0),
                    accuracy = z$accuracy[1],
                    sigma = z$sigma[1],
                    precision = z$precision[1])
          
        }
      }
    }
    
    abundance <- arrange(abundance, site_id, season) %>%
      group_by(season, ccamlr_id, season_relative, initial_season, last_season, accuracy, sigma, precision) %>%
      summarise(count = sum(count)) %>%
      ungroup()
    
    abundance_initial <- abundance %>%
      # keep first observed count for each site's time series
      dplyr::filter(initial_season == season_relative) %>%
      # ASSUMPTION: if no nest count is available in the initial season and a chick count is then
      # assume chick count is 1:1 nest count
      group_by(season, season_relative) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      dplyr::select(season, ccamlr_id, season_relative, count, sigma, precision)
    
    abundance_nests <- abundance %>%
      # keep all nest counts after the initial season
      dplyr::filter(initial_season != season_relative) %>%
      dplyr::select(season, ccamlr_id, season_relative, count, sigma, precision)
  
  } else {
    
    abundance_initial <- abundance %>%
      # keep first observed count for each site's time series
      dplyr::filter(initial_season == season_relative) %>%
      # ASSUMPTION: if no nest count is available in the initial season and a chick count is then
      # assume chick count is 1:1 nest count
      group_by(site_id, season, site, season_relative) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      dplyr::select(site_id, season, site, ccamlr_id, season_relative, count, sigma, precision)
    
    abundance_nests <- abundance %>%
      # keep all nest counts after the initial season
      dplyr::filter(initial_season != season_relative & type == "nests") %>%
      dplyr::select(site_id, season, site, ccamlr_id, season_relative, count, sigma, precision)
  
  } 
  
  
  
  dat <- list(abundance_nests = abundance_nests,
              abundance_initial = abundance_initial,
              site_list = site_list,
              n_seasons = n_seasons,
              seasons = seasons)
  
  return(dat)
  
}

sites <- c("PCHA", "LLAN", "PTHO", "CORM", "ARDL", "BISC", "CHIS", "HUMB", "TORG", "PETE", "LITC", "DETA",
           "ROYD", "BRDN", "BRDM", "BRDS", "CRZE", "CRZW", "INEX", "NORF", "EDMO", "WHEA", "CHAL",
           "CMID", "CNOR", "CSOU", "FRAE", "BEAU", "PGEO","BENT", "ONGU", "MAME", "RUMP", "YTRE", 
           "HUKU", "MIZU", "NOKK", "TORI", "BECH")

data_pop <- prep_pop_data("ADPE", 1979, 2018, sites = sites)

saveRDS(data_pop, "data_pop_adpe.rds")

# Stan data
data_stan_null <- list(
  nests = nrow(data_pop$abundance_nests),
  y_n = log(data_pop$abundance_nests$count),
  site_n = data_pop$abundance_nests$site,
  sigma_n = data_pop$abundance_nests$sigma,
  season_n = data_pop$abundance_nests$season_relative,
  y_i = log(data_pop$abundance_initial$count),
  sigma_i = data_pop$abundance_initial$sigma,
  n_sites = nrow(data_pop$site_list),
  n_seasons = data_pop$n_seasons,
  s = as.vector(data_pop$site_list$initial_season),
  reg = as.numeric(as.factor(data_pop$site_list$ccamlr_id)))

saveRDS(data_stan_null, "data_stan_null_adpe.rds")
