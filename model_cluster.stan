data {
  int<lower=1> N;                   // number of observations
  int<lower=1> D;                   // dimension of the target variables (e.g., D=5)
  int<lower=1> G;                   // number of groups (e.g., 16)
  int<lower=1> K;                   // truncation level for DP
  vector[D] Y[N];                   // observed target variables
  int<lower=1,upper=G> group[N];    // group indicator for each observation
  real<lower=0> alpha;              // DP concentration parameter
}

parameters {
  // Stick-breaking parameters
  real<lower=0,upper=1> v[K];

  // Cluster-level parameters
  matrix[D, K] mu_raw;      // raw cluster means
  vector[D] mu0;            // global prior mean
  cholesky_factor_corr[D] L; // Cholesky of correlation matrix for Y
  vector<lower=0>[D] sigma;  // std dev for each dimension

  // Group-level offsets
  matrix[D, G] gamma_raw;  // raw group offsets
}

transformed parameters {
  vector[K] pi;
  {
    vector[K] stick_elems;
    stick_elems[1] = v[1];

    for (k in 2:K) {
      real prod_term = 1.0;
      for (j in 1:(k-1)) {
        prod_term *= (1 - v[j]);
      }
      stick_elems[k] = v[k] * prod_term;
    }

    pi = stick_elems;  // mixture weights
  }

  // Construct actual means and offsets
  matrix[D,K] mu = mu_raw; // If needed, scale or shift mu_raw
  matrix[D,G] gamma = gamma_raw; // If needed, scale/shift offsets
}

model {
  // Priors
  // DP concentration
  // alpha can be fixed or have a prior (here fixed as data)
  
  // Stick-breaking
  for (k in 1:K) {
    v[k] ~ beta(1, alpha);
  }

  // Priors on cluster means
  // Let's assume mu_raw ~ normal(mu0, 1) and gamma_raw ~ normal(0,1) for simplicity
  mu0 ~ normal(0,5);
  to_vector(mu_raw) ~ normal(0,5);
  
  // Priors on group offsets
  to_vector(gamma_raw) ~ normal(0,1);

  // Covariance structure
  sigma ~ cauchy(0,2);
  L ~ lkj_corr_cholesky(2);

  // Likelihood
  {
    // For each observation, mixture likelihood
    // Marginalize cluster assignments:
    vector[K] lp;
    for (n in 1:N) {
      for (k in 1:K) {
        vector[D] mean_nk = mu[,k] + gamma[, group[n]];
        lp[k] = log(pi[k]) + multi_normal_cholesky_lpdf(Y[n] | mean_nk, diag_pre_multiply(sigma, L));
      }
      target += log_sum_exp(lp);
    }
  }
}
