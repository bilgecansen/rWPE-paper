
data {
  int<lower=1> n_sites;
  int<lower=1> n_seasons;
  int<lower=1> nests_sat;
  int<lower=1> nests_aer;
  int<lower=1> site_n_sat[nests_sat];
  int<lower=1> site_n_aer[nests_aer];
  int<lower=1> season_n_sat[nests_sat];
  int<lower=1> season_n_aer[nests_aer];
  real<lower=0> y_sat[nests_sat];
  int y_aer[nests_aer];
  real<lower=0> y_aer_test[nests_aer];
  int<lower=1> img_qual[nests_sat];
  real dos_aer[nests_aer];
  real dos_sat[nests_sat];
}

parameters {
  real alpha;
  real beta_dos1;
  real beta_dos2;
  real<lower=0> beta_sat[3];
  real<lower=0> mu_first;
  real mu_site_raw[n_sites];
  matrix[n_sites, n_seasons - 1] r_raw;
  real<lower=0> sigma_site;
  real<lower=0> sigma_aer;
  real<lower=0> sigma_r;
  real<lower=0> sigma_first;
  real<lower=0> sigma_sat[3];
  real lz_first[n_sites];
  real<lower=0> lambda[nests_aer];
}

transformed parameters {
  real<lower=0> sat_mean[nests_sat];
  real<lower=0> N[n_sites, n_seasons];
  real lz[n_sites, n_seasons];
  
  
  for (i in 1:n_sites) {
    lz[i,1] = lz_first[i];
    for (t in 1:(n_seasons-1)) {
      lz[i,t+1] = lz[i, t] + alpha + (mu_site_raw[i]*sigma_site) + (r_raw[i,t]*sigma_r);
    }
  }
  
  N = exp(lz);
  for (i in 1:nests_sat){
    sat_mean[i] = N[site_n_sat[i], season_n_sat[i]]*beta_sat[img_qual[i]]*exp(beta_dos1*dos_sat[i]);
  }
}

model {
  // Priors
  mu_first ~ uniform(0,50000);
  sigma_first ~ normal(0,2);
  lz_first ~ normal(log(mu_first), sigma_first);
  
  for (i in 1:n_sites) {
    r_raw[i,] ~ normal(0,1);
  }
  
  alpha ~ normal(0,0.2);
  beta_dos1 ~ normal(0,0.2);
  beta_dos2 ~ normal(0,0.2);
  beta_sat ~ normal(1,0.2);
  mu_site_raw ~ normal(0,1);
  sigma_site ~ normal(0,1);
  sigma_aer ~ normal(0,1);
  sigma_r ~ normal(0,1);
  sigma_sat ~ normal(0,1);
  
  //Nest counts after initial year
  for (i in 1:nests_aer) {
    lambda[i] ~ lognormal(lz[site_n_aer[i], season_n_aer[i]] + dos_aer[i]*beta_dos2 - 0.5*(sigma_aer^2), sigma_aer);
    y_aer[i] ~ poisson(lambda[i]);
  }
  
  for (i in 1:nests_sat){
    y_sat[i] ~ normal(sat_mean[i], sat_mean[i]*sigma_sat[img_qual[i]]);   
  }  
}

generated quantities {
  int<lower = 0, upper = 1> mean_gt_sat;
  int<lower = 0, upper = 1> sd_gt_sat;
  int<lower = 0, upper = 1> mean_gt_aer;
  int<lower = 0, upper = 1> sd_gt_aer;
  vector[nests_sat] y_rep_sat;
  vector[nests_aer] y_rep_aer;
  
  // Posterior predictive checks
  for (i in 1:nests_aer) {
    y_rep_aer[i] = poisson_rng(lambda[i]);
  }
  
  for (i in 1:nests_sat) {
    y_rep_sat[i] = normal_rng(sat_mean[i], sat_mean[i]*sigma_sat[img_qual[i]]);
  }
  
  mean_gt_sat = mean(y_rep_sat) > mean(y_sat);
  sd_gt_sat = sd(y_rep_sat) > sd(y_sat);
  mean_gt_aer = mean(y_rep_aer) > mean(y_aer_test);
  sd_gt_aer = sd(y_rep_aer) > sd(y_aer_test);
}

