data {
  int<lower=1> N; //number of subjects
  int<lower=1> P; //number of covariates
  int<lower=1> T; //number of total observations = dates(first patient) + dates(second) + ...
  int<lower=1> subj[N]; //subject id vector
  int<lower=1> counter[N]; //exact number of observations for every patient
  real y[T]; //outcome
  matrix[T,P] X; //predictors
}

parameters {
  vector[P] beta; //fixed intercept and slope
  vector[N] b; //subject intercepts (note: fixed for each patient, not time-dependent)
  real<lower=0> eta; //sd for subject intercepts
  real<lower=0> sigma_e; //error sd
  real mub;
}

model {
  //priors
  mub ~ normal(0, 2); //subj random effects
  eta ~ inv_gamma(3,2); //variance of intercepts
  b ~ normal(mub, eta); //subj random effects
  //w ~ normal(0, sigma_w); //item random effects
  beta ~ normal(0,5);
  sigma_e ~ normal(0,5);
  
  
  vector[T] mu;
  int index = 1;
  for (t in 1:T) {
    mu[t] = dot_product(X[t], beta) + b[subj[index]]; // Media prevista
    if(t==counter[index]){
     index = index + 1;
    }
  }
  
  y ~ normal(mu, sigma_e);
}
generated quantities {
  vector[T] log_lik; //log-likelihood
  vector[T] y_hat; 
 // vector[T] mu;
  //vector[T] y;
  
  int index = 1;
  for (t in 1:T) {
  //mu[t] = dot_product(X[t], beta) + b[subj[index]]; // Media prevista
  //y[t] ~ normal(mu[t], sigma_e);
   y_hat[t] = dot_product(X[t], beta) + b[subj[index]];
    log_lik[t] = normal_lpdf(y[t] | y_hat[t], sigma_e);
     if(t==counter[index]){
      index = index + 1;
    }
  }
}
