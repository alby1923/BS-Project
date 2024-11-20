library(progress)
library(ggplot2)
library(rstan)

df_stan <- readRDS('filled_dataset.rds')
#sort patients by CAI
#df_stan <- df_stan[order(df_stan[, 1]), ]
#STAN MODEL -------------

#need to keep only patients with at least 4 observations (TO DO)

#linear model (or maybe other ideas for variable selection) to cut non-significant variables (TO DO)

data_stan <- function(df,variables_to_keep,response){
#data needed for the model
patients_mat <- as.matrix(unique(df[,1]))
patients <- as.vector(patients_mat)
N <- dim(patients_mat)[1]
T <- dim(df)[1]
data <- response #save response data
covariates <- df[,variables_to_keep] #remove response (and useless) covariates

#tot_obs_patient <- table(df$CAI)
tot_obs_patient <- as.vector(table(df$CAI)[match(patients_mat, names(table(df$CAI)))])


data_model = list(N = N, P = dim(covariates)[2], T = T, 
                     subj = patients, 
                     y = data, counter = tot_obs_patient, X = covariates)

return(data_model)
}

#keep just significan covariates and cut responses
chosen_columns <- c(
  "Alanina_aminotransferasi_alt",
  #"Albumina",
  "Altezza",
  "Colesterolo_Hdl",
  "Colesterolo_totale",
  "Creatinina",
  "Distribuzione_di_volume",
  "Ematocrito_hct",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_hb",
  "Emoglobina_massa_media_mch",
  "Eosinofili_perc",
  "Eritrociti_rbc",
  "Ferritina",
  #"Ferro_totale",
  #"Glucosio",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Monociti_perc",
  "PMAX",
  "Peso",
  "Piastrine",
  "Polso",
  "Proteine_totali",
  #"S_alfa_1_globuline",
  #"S_alfa_2_globuline",
  #"S_beta_1_globuline",
  #"S_beta_2_globuline",
  #"S_gamma_globuline",
  "Trigliceridi",
  "Volume_medio",
  #  "Alcool",
  #  "Attivita_fisica",
  "Circonferenza_vita"
  #  "Fumo"
)

response <- df_stan$Glucosio #CHANGE THIS FOR OTHER RESPONSES
data_for_model <- data_stan(df_stan,chosen_columns,response)

#change values for simulation
fit = stan(file = 'model_base.stan', 
                    data = data_for_model, 
                    chains = 1, 
                    iter = 5, 
                    warmup = 1, 
                    thin = 1, 
                    seed = 19)

summary(fit, pars = c("beta", "sigma_e"))
traceplot(fit, pars = c(paste0("beta[", 1:20, "]"), "sigma_e"))
print(fit, pars = c("beta", "sigma_e"))
#parallel::detectCores()

#extract results
posterior_samples <- extract(fit)
#access to y_hat
y_hat <- posterior_samples$y_hat
y_hat_last_chain <- y_hat[4000,] #we forgot to ask about this...
#obtain residuals
res <- response - y_hat_last_chain

#group residuals and y_hat for each patient
patient_id <- df_stan$CAI
res_grouped <- split(res, patient_id)
#y_hat_grouped <- split(y_hat_last_chain, patient_id)

#study autocorrelation of residuals (timestamp needed) (TO DO)

#add something to see convergence of chains??

#stuff to assess how the model performed (need to add something)
rss <- sum(res^2)
rmse <- sqrt(rss/dim(df_stan)[1])
sd_y <- sd(response)
RMSE_to_sd_ratio <- rmse / sd_y
