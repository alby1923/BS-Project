library(progress)
library(ggplot2)
library(rstan)

df_stan <- readRDS('filled_dataset.rds')

#FILTER DATASET (OBSERVATIONS AND COVARIATES) -------------
#need to keep only patients with at least 4 observation
minimum_observations <- function(df, threshold) {
  
  # Group by patient ID (CAI) and count the number of observations for each patient
  patient_counts <- table(df$CAI)
  
  # Filter patients who have at least the threshold number of observations
  valid_patients <- names(patient_counts[patient_counts >= threshold])
  
  # Subset the dataset to include only those patients who meet the threshold
  df_filtered <- df[df$CAI %in% valid_patients, ]
  
  # Return the filtered dataset
  return(df_filtered)
}

df_stan_4obs <- minimum_observations(df_stan,4)

#linear model (or maybe other ideas for variable selection) to cut non-significant variables
model <- lm(Glucosio~ . - CAI - Date - delta_date - Glucosio - Colesterolo_Hdl - PMAX - Trigliceridi - Circonferenza_vita, data = df_stan)
summary(model)
#comment below discarded variable
#STAN MOEL -------
data_stan <- function(df,variables_to_keep,response){
#data needed for the model
patients_mat <- as.matrix(unique(df[,1]))
patients <- as.vector(patients_mat)
N <- dim(patients_mat)[1]
T <- dim(df)[1]
#data <- as.vector(df[,response]) #save response data
covariates <- df[,variables_to_keep] #remove response (and useless) covariates

#tot_obs_patient <- table(df$CAI)
tot_obs_patient <- as.vector(table(df$CAI)[match(patients_mat, names(table(df$CAI)))])
subj <- rep(seq_along(tot_obs_patient), times = tot_obs_patient)

data_model = list(N = N, P = dim(covariates)[2], T = T, 
                     subj = subj, y = response, X = covariates)

return(data_model)
}

#keep just significant covariates and cut responses
chosen_columns <- c(
  "Altezza",
  #"Colesterolo_Hdl",
  "Colesterolo_totale",
  "Distribuzione_di_volume",
  "Ematocrito_hct",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_hb",
  "Emoglobina_massa_media_mch",
  "Eosinofili_perc",
  "Eritrociti_rbc",
  #"Glucosio",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Monociti_perc",
  #"PMAX",
  "Peso",
  "Piastrine",
  "Polso",
  #"Trigliceridi",
  "Volume_medio",
    "Alcool",
    "Attivita_fisica",
  #"Circonferenza_vita",
  "Fumo",
  "Rh",
  "SESSO",
  "AB0"
  #,"delta_date"
)

response <- as.vector(df_stan$Glucosio) #CHANGE HERE!
data_for_model <- data_stan(df_stan,chosen_columns,response)

parallel::detectCores()
#change values for simulation
fit = stan(file = 'model_base.stan', 
                    data = data_for_model, 
                    chains = 2, 
                    iter = 2000, 
                    warmup = 1000, 
                     cores = 2,
                    thin = 1, 
                    seed = 19)

summary(fit, pars = c("beta", "sigma_e"))
traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
print(fit, pars = c("beta", "sigma_e"))
#add something to see convergence of chains??

#RESULT ANALYSIS -------------

#extract results
posterior_samples <- extract(fit)
#access to y_hat
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE) #can be done mean, median or mode. matrix is number of iterations * number of observations
#obtain residuals
res <- response - y_hat_point

#stuff to assess how the model performed (need to add something and eventually fix this part)
rss <- sum(res^2)
rmse <- sqrt(rss/dim(df_stan)[1])
sd_y <- sd(response)
RMSE_to_sd_ratio <- rmse / sd_y

#RESIDUAL ANALYSIS -----------

#group residuals and y_hat for each patient
patient_id <- df_stan$CAI
res_grouped <- split(res, patient_id)
y_hat_grouped <- split(y_hat_point,patient_id)

#study autocorrelation of residuals

#filter patients ith at least 10 observations
res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 10]
#associates timestamps
res_grouped_with_timestamp <- lapply(names(res_grouped_filtered), function(patient_id) {
  patient_data <- res_grouped_filtered[[patient_id]]  # Dati del paziente
  # Ottieni i timestamp per quel paziente
  patient_timestamps <- df_stan[df_stan$CAI == patient_id, "delta_date"]
  
  # Associa timestamp ai dati del paziente
  data.frame(observations = patient_data, timestamp = patient_timestamps)
})
# COMPUTE AUTOCORRELATION (TO DO)




