library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('filled_dataset.rds')

df_stan <- na.omit(df_stan) #4 patients had no DATA_NASCITA so we couldn't compute their age

#STAN MODEL -------
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
response <- as.vector(df_stan$Trigliceridi) #CHANGE HERE!
data_for_model <- data_stan(df_stan,chosen_columns,response)

parallel::detectCores()
#change values for simulation
fit = stan(file = 'model_base.stan', 
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

folder_name <- paste0("STAN_Model")
#dir.create(folder_name)
file_path <- file.path(folder_name, paste0("fit_Trigliceridi_new.Rdata"))
save(fit, file = file_path)

#RESULT ANALYSIS -------------
target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
res_list <- list()
for(target in target_variables){
  dir.create(paste0('PLOT_',target))
  folder_name <- paste0('PLOT_',target)
  cat('target: ',target)
  load(paste0("fit_",target,"_new.RData"))
  response <- df_stan[,target]
  #extract results and obtain residuals
  posterior_samples <- extract(fit)
  y_hat <- posterior_samples$y_hat
  y_hat_point <- colMeans(y_hat, na.rm = TRUE) #matrix is number of iterations * number of observations
  #y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) #try with the median
  res <- t(response - y_hat_point)
  res_list[[target]] <- res
  #analyze performance
  rss <- sum(res^2)
  cat('\nRSS: ',rss)
  rmse <- sqrt(rss/dim(df_stan)[1]) 
  cat('\nRMSE: ',rmse)
  RMSE_to_sd_ratio <- rmse / sd(t(response))
  cat('\nRMSE/SD: ',RMSE_to_sd_ratio)
  
  # Posterior di sigma
  sigma_posterior <- posterior_samples$sigma
  
  # Plotta la distribuzione a posteriori di sigma e colorala di azzurro
  q <- stan_dens(fit, pars = "sigma_e")
  q + 
    ggtitle(paste("Sigma posterior distribution of", target)) +  # Titolo personalizzato
    theme_minimal() +  # Usa un tema minimale per il grafico
    geom_density(fill = "skyblue", color = "blue")  # Colora la densità in azzurro
  ggsave(paste0(folder_name, "/Sigma posterior distribution of", target, ".png"), plot = q, width = 6, height = 4)
  
  # Plotta la distribuzione a posteriori dei beta e colorala di azzurro
  h <- stan_dens(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]")))
  h + 
    ggtitle(paste("Beta posterior distribution of", target)) +  # Titolo personalizzato
    theme_minimal() +  # Tema minimale
    geom_density(fill = "skyblue", color = "blue")  # Colora la densità in azzurro
  ggsave(paste0(folder_name, "/Beta posterior distribution of", target, ".png"), plot = h, width = 6, height = 4)
  
  p <- ggplot(data.frame(t(res)), aes(x = t(res))) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1.5) +  # Aggiungi la curva di densità
    #stat_function(fun = dnorm, args = list(mean = mean(res), sd = sd(res)), color = "green", size = 1.5) +  # Aggiungi una curva normale teorica
    labs(title = paste("Residuals Histogram", target), x = "Valore", y = "Densità") +
    theme_minimal()
  plot(p)
  # Salva il grafico come immagine PNG
  ggsave(paste0(folder_name, "/residuals_histogram_", target, ".png"), plot = p, width = 6, height = 4)
  
  #print(summary(fit, pars = c("beta", "sigma_e")))
  #q <- traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
  #plot(q)
  #print(fit, pars = c("beta", "sigma_e"))

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
  
  cat('\n-------------------------------\n\n\n')
}