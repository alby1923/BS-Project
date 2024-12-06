// Model: Multivariate Hierarchical Regression with Shared Donor Effects and Sequential Dependencies

data {
  int<lower=1> N;          // Number of donors
  int<lower=1> P;          // Number of covariates
  int<lower=1> T;          // Total number of observations
  int<lower=1> K;          // Number of target variables
  
  int<lower=1, upper=N> subj[T];    // Donor ID for each observation
  matrix[T, P] X;                    // Predictor matrix
  matrix[T, K] y;                    // Outcome matrix (T observations, K targets)
}

parameters {
  matrix[K, K] gamma_raw;   // Unconstrained gamma matrix for priors (includes upper triangular)
  matrix[P, K] beta;        // Coefficients for each covariate and target
  vector[N] alpha;          // Donor-specific intercepts
  vector[K] phi;            // Target-specific intercepts
  
  real mub;                 // Hyperparameter: mean of alpha
  real<lower=0> eta;        // Hyperparameter: SD of alpha
  
  real<lower=0> sigma_e[K]; // Error SDs for each target
}

transformed parameters {
  matrix[K, K] gamma;       // Lower triangular gamma matrix
  
  // Set lower triangular values to sampled values, upper triangular to 0
  for (k in 1:K) {
    for (j in 1:K) {
      if (j < k) {
        gamma[k, j] = gamma_raw[k, j]; // Retain sampled value for lower triangular
      } else {
        gamma[k, j] = 0;              // Set upper triangular to 0
      }
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
  
  // Priors for gamma_raw (lower triangular)
  for (k in 1:K) {
    for (j in 1:K) {
      if (j < k) {
        gamma_raw[k, j] ~ normal(0, 5); // Prior for lower triangular values
      }
    }
  }
  
  // Priors for error standard deviations
  for (k in 1:K) {
    sigma_e[k] ~ inv_gamma(3, 2);       // Error SD priors
  }
  
  // Likelihood
  for (k in 1:K) {
    // Initialize the linear predictor
    vector[T] mu_k = X * beta[, k] + alpha[subj] + phi[k];
    
    // Add dependencies on previous Y's only if k > 1
    if (k > 1) {
      for (j in 1:(k-1)) {
        mu_k += gamma[k, j] * y[, j];
      }
    }
    
    // Define the likelihood
    y[, k] ~ normal(mu_k, sigma_e[k]);
  }
}

generated quantities {
  matrix[T, K] y_hat;        // Predicted outcomes
  
  for (k in 1:K) {
    // Initialize the linear predictor
    y_hat[, k] = X * beta[, k] + alpha[subj] + phi[k];
    
    // Add dependencies on previous Y's only if k > 1
    if (k > 1) {
      for (j in 1:(k-1)) {
        y_hat[, k] += gamma[k, j] * y_hat[, j];
      }
    }
  }
}
