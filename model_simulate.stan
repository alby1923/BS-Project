data {
  int<lower=1> N; // number of subjects
  int<lower=1> P; // number of covariates
  int<lower=1> T; // number of total observations
  int<lower=1> subj[T]; // subject id vector
 
  int<lower=0> n_obs; // number of observed y
  int<lower=0> n_miss; // number of missing y
  int<lower=1, upper=T> obs_indices[n_obs]; // indices of observed y
  int<lower=1, upper=T> miss_indices[n_miss]; // indices of missing y
  vector[n_obs] y_obs; // observed y values
  
  matrix[n_obs, P] X1; // predictors observed values y
  matrix[n_miss, P] X2; // predictors missing values y
}

parameters {
  vector[P] beta; // fixed intercept and slope
  vector[N] b; // subject intercepts
  real<lower=0> eta; // sd for subject intercepts
  real<lower=0> sigma_e; // error sd
  real mub; // mean for subject intercepts
}

model {
  // Priors
  mub ~ normal(0, 2); // subj random effects
  eta ~ inv_gamma(3, 2); // variance of intercepts
  b ~ normal(mub, eta); // subj random effects
  beta ~ normal(0, 5); // fixed effects
  sigma_e ~ inv_gamma(3, 2); // error variance

  // Linear predictor for observed y (just estimate mu on known target values)
  vector[n_obs] mu_obs;
  mu_obs = X1 * beta + b[subj[obs_indices]];

  // Likelihood for observed y
  y_obs ~ normal(mu_obs, sigma_e);
  
}

generated quantities {
  vector[n_obs] y_obs_hat;  // Predictions for observed y
  vector[n_miss] y_miss_hat; // Predictions for missing y
  vector[T] y_hat; // Combined predicted y for all observations
  
  // Fill in predictions for observed y
  y_obs_hat = X1 * beta + b[subj[obs_indices]];  // computes y_obs_hat from linear model
  y_hat[obs_indices] = y_obs_hat;

  // Compute mu_miss for missing y
  vector[n_miss] mu_miss;
  mu_miss = X2 * beta + b[subj[miss_indices]];

  // Since we want to simulate data, we sample y_miss_hat from a normal distribution with mean mu_miss and var sigma_e
  for (i in 1:n_miss) {
    y_miss_hat[i] = normal_rng(mu_miss[i], sigma_e);
  }

  // Fill in y_hat for missing indices
  y_hat[miss_indices] = y_miss_hat;
}



