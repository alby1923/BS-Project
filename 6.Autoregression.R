library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('filled_dataset.rds')

#df_stan <- na.omit(df_stan)

#OBTAIN DATAFRAME ----------------
#for simulation just remove columns with too many NA but you can keep target variables with NA

#STAN MODEL -------
data_stan <- function(df,variables_to_keep,response){
  #data needed for the model
  patients_mat <- as.matrix(unique(df[,1]))
  patients <- as.vector(patients_mat)
  N <- dim(patients_mat)[1]
  T <- dim(df)[1]
  #data <- as.vector(df[,response]) #save response data
  covariates <- df[,variables_to_keep] #remove response (and useless) covariates
  
  df <- df %>%
    group_by(CAI) %>%
    arrange(delta_date, .by_group = TRUE) %>%
    mutate(
      successive_timestamps = if_else(
        delta_date == 0,
        0,               
        delta_date - lag(delta_date, default = 0)
      ),
      ar_flag = if_else(
        successive_timestamps <= 180 & successive_timestamps != 0, 
        1, 
        0
      ))%>%
    ungroup()  
  
  y_trasl <- lag(response, default = 0) * df$ar_flag
  
  #y_trasl_matrix <- matrix(y_trasl, nrow = 1, ncol = length(y_trasl))
  
  #tot_obs_patient <- table(df$CAI)
  tot_obs_patient <- as.vector(table(df$CAI)[match(patients_mat, names(table(df$CAI)))])
  subj <- rep(seq_along(tot_obs_patient), times = tot_obs_patient)
  
  data_model = list(N = N, P = dim(covariates)[2], T = T, 
                    subj = subj, y = response, y_trasl = y_trasl, X = covariates)
  
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
  "AB0",
  'eta_std'
  #,"delta_date"
)
#Colesterolo_Hdl, Circonferenza_vita, Glucosio, PMAX, Trigliceridi, trained one each (no for loop due to problems with saving traceplot)
response <- as.vector(df_stan$Colesterolo_Hdl) #CHANGE HERE!
data_for_model <- data_stan(df_stan,chosen_columns,response)

parallel::detectCores()
#change values for simulation
fit = stan(file = 'autoregressive_model.stan', 
           data = data_for_model, 
           chains = 4, 
           iter = 2000, 
           warmup = 1000, 
           cores = 4,
           thin = 1,
           #control = list(adapt_delta = 0.95), may increase running time but better mix of chains. minimum value = 0.8 (default)
           seed = 19)

summary(fit, pars = c("beta", "sigma_e"))
traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
print(fit, pars = c("beta", "sigma_e"))

#extract results and obtain residuals
posterior_samples <- extract(fit)
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE) #matrix is number of iterations * number of observations
#y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) #try with the median
res <- t(response - y_hat_point)

#analyze performance
rss <- sum(res^2)
cat('\nRSS: ',rss)
rmse <- sqrt(rss/dim(df_stan)[1]) 
cat('\nRMSE: ',rmse)
RMSE_to_sd_ratio <- rmse / sd(t(response))
cat('\nRMSE/SD: ',RMSE_to_sd_ratio)

# Calcola R² a posteriori
response <- t(response)
SS_tot <- sum((response - mean(response))^2)  
SS_res <- apply(y_hat, 1, function(y_h) sum((response - y_h)^2))
R2_posterior <- 1 - SS_res / SS_tot

# Calcola la media di R² e l'intervallo di credibilità
mean_R2 <- mean(R2_posterior)
R2_ci <- quantile(R2_posterior, probs = c(0.025, 0.975))

cat("\nMedia di R² a posteriori: ", round(mean_R2, 3), "\n")
cat("Intervallo di credibilità del 95% per R²: ", round(R2_ci, 3), "\n")