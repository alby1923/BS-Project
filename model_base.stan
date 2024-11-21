data {
  int<lower=1> N; //number of subjects
  int<lower=1> P; //number of covariates
  int<lower=1> T; //number of total observations = dates(first patient) + dates(second) + ...
  int<lower=1> subj[T]; //subject id vector
  //int<lower=1> counter[N]; //exact number of observations for every patient
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
  mu = X * beta + b[subj];
  
  y ~ normal(mu, sigma_e);
}
generated quantities {
  //vector[T] log_lik; //log-likelihood
  //real log_lik;
  vector[T] y_hat; 
   vector[T] mu;

   y_hat = X * beta + b[subj];
   //log_lik = normal_lpdf(y | y_hat, sigma_e); //prima passava tutto individuale, ora la somma, è ok?

}
