data {
  int<lower=1> N; // number of subjects
  int<lower=1> P; // number of covariates
  int<lower=1> T; // number of total observations
  int<lower=1> subj[T]; // subject id vector train
  vector[T] y; // train y values
  matrix[T, P] X; // predictors train values y
  
  matrix[N, P] Xtest; // predictors test values y
  int<lower=1> subj_test[N]; // subject id vector test
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
  vector[T] mu = X * beta + b[subj];

  // Likelihood for test y
  y ~ normal(mu, sigma_e);
  
}

generated quantities {
  vector[T] y_train;  // Predictions for train y
  vector[N] y_test; // Predictions for test y
  
  //y train predictions
  y_train = X * beta + b[subj];

  //y test predictions
  y_test = Xtest * beta + b[subj_test];
}

