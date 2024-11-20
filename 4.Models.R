library(progress)
library(ggplot2)
library(rstan)

#STAN MODEL (GLUCOSIO)-------------
df_glucosio <- readRDS('notNAresponses.rds')

#remove all NA covariates
#df_glucosio <- df_glucosio[,-c(34,35,36,37)]

df_glucosio <- df_glucosio[order(df_glucosio[, 1]), ]
#data needed for the model
patients_glucosio_mat = as.matrix(unique(df_glucosio[,1]))
patients_glucosio <- as.vector(patients_glucosio_mat)
N_glucosio = dim(patients_glucosio_mat)[1]
T_glucosio = dim(df_glucosio)[1]
data = df_glucosio$Glucosio #save response data
covariates = df_glucosio[,-c(1,2,5,16,20,25,27)] #remove response (and useless) covariates

tot_obs_patient <- table(df_glucosio$CAI)
tot_obs_patient <- as.vector(tot_obs_patient[match(patients_glucosio_mat, names(tot_obs_patient))])


data_glucosio = list(N = N_glucosio, P = dim(covariates)[2], T = T_glucosio, 
                     subj = patients_glucosio, 
                     y = data, counter = tot_obs_patient, X = covariates)

fit_glucosio = stan(file = 'model_base.stan', 
                    data = data_glucosio, 
                    chains = 1, 
                    iter = 500, 
                    warmup = 100, 
                    thin = 1, 
                    seed = 123)
