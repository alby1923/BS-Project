// Model: Multivariate Hierarchical Regression with Shared Donor Effects and Sequential Dependencies

data {
  int<lower=1> N;          // Number of donors
  int<lower=1> P;          // Number of covariates
  int<lower=1> T;          // Total number of observations
  int<lower=1> K;          // Number of target variables
  
  int<lower=1, upper=N> subj[T];    // Donor ID for each observation
  matrix[T, P] X;                    // Predictor matrix
  matrix[T, K] y;                    // Outcome matrix (T observations, K targets)
  matrix[T, K] y_trasl; //autoregressive component
}

parameters {
  matrix[K, K] gamma_raw;   // Unconstrained gamma matrix for priors (includes upper triangular)
  matrix[P, K] beta;        // Coefficients for each covariate and target
  vector[N] alpha;          // Donor-specific intercepts
  vector[K] phi;            // Target-specific intercepts
  
  real mub;                 // Hyperparameter: mean of alpha
  real<lower=0> eta;        // Hyperparameter: SD of alpha
  vector[K] autoreg_coef_raw; // Y_trasl will be a vector of zeros for three targets
  real<lower=0> sigma_e[K]; // Error SDs for each target
}

transformed parameters {
  matrix[K, K] gamma;       // Lower triangular gamma matrix
  vector[K] autoreg_coef; // Set to zero coefficients for targets with no autoregressive component
  
  // Set lower triangular values to sampled values, upper triangular to 0
  for (k in 1:K) {
    for (j in 1:K) {
      if (j > k) {
        gamma[k, j] = gamma_raw[k, j]; // Retain sampled value for upper triangular
      } else {
        gamma[k, j] = 0;              // Set lower triangular to 0
      }
    }
  }
  
  for(k in 1:K){
    if (k==2){
      autoreg_coef[k] = autoreg_coef_raw[k];
    }
    else {
      autoreg_coef[k] = 0;
    }
  }
}

model {
  // Priors
  mub ~ normal(0, 2);
  eta ~ inv_gamma(3, 2);
  
  alpha ~ normal(mub, eta);            // Donor-specific effects
  phi ~ normal(0, 5);                  // Target-specific intercepts
  to_vector(beta) ~ normal(0, 5);      // Coefficients for covariates
  
  // Priors for gamma_raw (upper triangular)
  for (k in 1:K) {
    for (j in 1:K) {
      if (j > k) {
        gamma_raw[k, j] ~ normal(0, 5); // Prior for upper triangular values
      }
    }
  }
  
  // Priors for error standard deviations
  for (k in 1:K) {
    sigma_e[k] ~ inv_gamma(3, 2);       // Error SD priors
  }
  
  // Prior for autoregressive coefficient
  for (k in 1:K) {
   autoreg_coef_raw[k] ~ normal(0, 2); //we know from previous experiments this value is really near zero
  }
  
  // Likelihood
  for (k in 1:K) {
    // Initialize the linear predictor
    vector[T] mu_k = X * beta[, k] + alpha[subj] + phi[k] + autoreg_coef[k] * to_vector(y_trasl[, k]) + y * gamma[,k];
    
    // Define the likelihood
    y[, k] ~ normal(mu_k, sigma_e[k]);
  }
}

generated quantities {
  matrix[T, K] y_hat;        // Predicted outcomes
  
  for (k in 1:K) {
    // Initialize the linear predictor
    y_hat[, k] = X * beta[, k] + alpha[subj] + phi[k] + autoreg_coef[k] * to_vector(y_trasl[, k]);
  }
  
  // Add contributions from dependencies on previous targets
  for (k in 2:K) {
    y_hat[:, k] += y_hat * gamma[, k];
  }
  
}
