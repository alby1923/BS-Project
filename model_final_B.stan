data {
  int<lower=1> N; // number of subjects
  int<lower=1> P; // number of covariates
  int<lower=1> T; // number of total observations
  int<lower=1> K; // previous targets
  int<lower=1> subj[T]; // subject id vector train
  vector[T] y; // train y values
  matrix[T, P] X; // predictors train values y
  int<lower=0> lag; // threshold for autocorrelation
  matrix[T,K] y_std; // previous targets standardized
  matrix[T,K] y_pred_train_std; // previous predicted targets standardized
  vector[T] y_trasl; // autoregression values
  
  matrix[N, P] Xtest; // predictors test values y
  int<lower=1> subj_test[N]; // subject id vector test
  matrix[N,K] y_pred_test_std; // previous predicted targets standardized
  vector[N] y_trasl_test; // autoregression values
}

parameters {
  vector[P] beta; // fixed intercept and slope
  vector[N] b; // subject intercepts
  vector[K] gamma; // targets parameters
  real<lower=0> eta; // sd for subject intercepts
  real<lower=0> sigma_e; // error sd
  real mub; // mean for subject intercepts
  real phi_raw; //autoregressive coefficient
}

transformed parameters {
  real phi;
  
  //if no autoregressive component set it to zero
  if(lag!=0){
      phi = phi_raw;
    }
  else {
    phi = 0;
    }
}

model {
  // Priors
  mub ~ normal(0, 2); // subj random effects
  eta ~ inv_gamma(3, 2); // variance of intercepts
  b ~ normal(mub, eta); // subj random effects
  beta ~ normal(0, 5); // fixed effects
  sigma_e ~ inv_gamma(3, 2); // error variance
  gamma ~ normal(0, 5); // previous targets effects
  phi_raw ~ normal(0, 1); // from previous experiments we know this value is really low
  
  // Linear predictor for observed y (just estimate mu on known target values)
  vector[T] mu = X * beta + b[subj] + (phi * y_trasl) + y_std * gamma; 
  
  //values associated to autoregressive component NOT standardized since output is log-transformed and here we are computing autoregressive value
  
  // Likelihood for test y
  y ~ normal(mu, sigma_e);
  
}

generated quantities {
  vector[T] y_train;  // Predictions for train y
  vector[N] y_test; // Predictions for test y
  
  //y train predictions
  y_train = X * beta + b[subj] + (phi * y_trasl) + y_pred_train_std * gamma;

  //y test predictions
  y_test = Xtest * beta + b[subj_test] + (phi * y_trasl_test) + y_pred_test_std * gamma;
}

