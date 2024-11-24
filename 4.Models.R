library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('filled_dataset.rds')
# counts how many dates are before the reference date
reference_date <- as.Date("2019-01-01")
dates_before <- sum(df_stan$Date < reference_date, na.rm = TRUE)
dates_before #many variables started being tracked from 2019 on
#compute the delta time from first osbervation for each patient
df_stan <- df_stan[order(df_stan$CAI), ]
calculate_delta_date <- function(df, patient_col, date_col) {
  
  df$delta_date <- unlist(by(
    df, 
    df[[patient_col]], 
    function(sub_df) {
      
      first_date <- sub_df[[date_col]][1]
      
      delta <- as.numeric(difftime(sub_df[[date_col]], first_date, units = "days"))
      return(delta)
    }
  ))
  
  return(df)
}
df_stan <- calculate_delta_date(df_stan, "CAI", "Date")

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
#Colesterolo_Hdl, Circonferenza_vita, Glucosio, PMAX, Trigliceridi, trained one each (no for loop due to problems with saving traceplot)
response <- as.vector(df_stan$PMAX) #CHANGE HERE!
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
                    seed = 19)

summary(fit, pars = c("beta", "sigma_e"))
traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
print(fit, pars = c("beta", "sigma_e"))

folder_name <- paste0("STAN_Model")
#dir.create(folder_name)
file_path <- file.path(folder_name, paste0("fit_PMAX_new.Rdata"))
save(fit, file = file_path)

#RESULT ANALYSIS -------------
target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
res_list <- list()
folder_name <- paste0("STAN_Model")
for(target in target_variables){
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
  
  p <- ggplot(data.frame(t(res)), aes(x = t(res))) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1.5) +  # Aggiungi la curva di densità
    #stat_function(fun = dnorm, args = list(mean = mean(res), sd = sd(res)), color = "green", size = 1.5) +  # Aggiungi una curva normale teorica
    labs(title = paste("Residuals Histogram", target), x = "Valore", y = "Densità") +
    theme_minimal()
  # Salva il grafico come immagine PNG
  ggsave(paste0(folder_name, "/residuals_histogram_", target, ".png"), plot = p, width = 6, height = 4)
  cat('\n-------------------------------\n\n\n')
}
#AUTOCORRELATION ---------
fake_pacf <- function(res_grouped_with_timestamps) {
  # Creare una lista vuota per raccogliere i valori per ciascun lag
  pacf_results <- list()
  
  # Loop attraverso ciascun gruppo (time series + timestamps)
  for (CAI in names(res_grouped_with_timestamps)) {
    group <- res_grouped_with_timestamps[[CAI]]
    ts <- group$results
    ts_time <- group$timestamps
    
    # Calcolare la media della serie temporale
    mu <- mean(ts)
    
    # Calcolare il massimo lag possibile
    max_lag <- max(ts_time) - min(ts_time)
    
    for (lag in 1:max_lag) {
      # Trova tutte le coppie di timestamp distanziate da 'lag'
      pairs <- combn(seq_along(ts_time), 2, simplify = FALSE)
      valid_pairs <- Filter(function(p) abs(ts_time[p[2]] - ts_time[p[1]]) == lag, pairs)
      
      if (length(valid_pairs) > 0) {
        # Calcolare la covarianza per questo lag
        cov <- sum(sapply(valid_pairs, function(p) {
          (ts[p[1]] - mu) * (ts[p[2]] - mu)
        })) / length(valid_pairs)
        
        # Calcolare la varianza della serie temporale
        var_ts <- sum((ts - mu)^2) / (length(ts) - 1)
        
        # Calcolare la correlazione (covarianza / varianza)
        correlation <- cov / var_ts
        
        # Aggiungere il valore alla lista del lag corrispondente
        lag_key <- paste0("lag", lag)
        pacf_results[[lag_key]] <- c(pacf_results[[lag_key]], correlation)
      }
    }
  }
  
  pacf_df <- data.frame(
    lag = sprintf("%02d", as.numeric(gsub("lag", "", names(pacf_results)))),  # Assegna i lag con il formato "01", "02", ...
    mean_value = sapply(pacf_results, mean, na.rm = TRUE),
    count = sapply(pacf_results, length)
  )
  
  return(pacf_df)
}
#loop for on residuals
target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
df_list <- list()
for(target in target_variables){
  cat('\nstarting: ',target)
  res <- res_list[[target]]
  #group residuals to each patient and take only time series with at least 3 obs
  patient_id <- df_stan$CAI
  res_grouped <- split(res, patient_id)
  res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 3]
  
  #associate timestamp to each time series
  res_grouped_with_filtered_timestamps <- lapply(names(res_grouped_filtered), function(CAI) {
    group <- res_grouped_filtered[[CAI]]
    relevant_timestamps <- unname(df_stan$delta_date[df_stan$CAI == CAI]%/%90)
    list(
      results = group,         
      timestamps = relevant_timestamps  
    )
  })
  names(res_grouped_with_filtered_timestamps) <- names(res_grouped_filtered)
  
  #compute the fake PACF for irregular timestamps
  df_model <- fake_pacf(res_grouped_with_filtered_timestamps)
  df_model <- rbind(data.frame(lag = "0", mean_value = 1, count = 0), df_model)
  #df$lag <- factor(df$lag, levels = c("lag0", paste0("lag", 1:max(as.numeric(gsub("lag", "", df$lag))))))
  df_list[[target]] <- df_model
}

df_vita <- df_list[['Circonferenza_vita']]
df_col <- df_list[["Colesterolo_Hdl"]]
df_gluc <- df_list[["Glucosio"]]
df_pmax <- df_list[["PMAX"]]
df_trig <- df_list[["Trigliceridi"]]

merged_data <- full_join(df_vita, df_col, by = "lag", suffix = c("_vita", "_col")) %>%
  full_join(df_gluc, by = "lag", suffix = c("_col", "_gluc")) %>%
  full_join(df_pmax, by = "lag", suffix = c("_gluc", "_pmax")) %>%
  full_join(df_trig, by = "lag", suffix = c("_pmax", "_trig"))

#substitute NA with zero
merged_data <- merged_data %>%
  mutate(
    mean_value_vita = ifelse(is.na(mean_value_vita), 0, mean_value_vita),
    mean_value_col = ifelse(is.na(mean_value_col), 0, mean_value_col),
    mean_value_gluc = ifelse(is.na(mean_value_gluc), 0, mean_value_gluc),
    mean_value_pmax = ifelse(is.na(mean_value_pmax), 0, mean_value_pmax),
    mean_value = ifelse(is.na(mean_value), 0, mean_value),
    count_vita = ifelse(is.na(count_vita), 0, count_vita),
    count_col = ifelse(is.na(count_col), 0, count_col),
    count_gluc = ifelse(is.na(count_gluc), 0, count_gluc),
    count_pmax = ifelse(is.na(count_pmax), 0, count_pmax),
    count = ifelse(is.na(count), 0, count)
  )

#compute weighted values
result <- merged_data %>%
  mutate(
    weighted_value = (mean_value_vita * count_vita + mean_value_col * count_col + mean_value_gluc * count_gluc + 
                        mean_value_pmax * count_pmax + mean_value * count) / 
      (count_vita + count_col + count_gluc + count_pmax + count),
    total_count = count_vita + count_col + count_gluc + count_pmax +count
  ) %>%
  select(lag, weighted_value, total_count)

#pacf lag0 always equal to 1
result[1, "weighted_value"] = 1

# Calcolare max e min value escluso 1 per le linee tratteggiate nere
max_value <- max(result$weighted_value[result$weighted_value != 1])
min_value <- min(result$weighted_value[result$weighted_value != 1])

# Crea un grafico a dispersione per visualizzare i valori PACF medi
ggplot(result, aes(x = lag, y = weighted_value)) +
  geom_point(color = "blue", size = 3) +  # Punti per ogni lag
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Linea rossa a zerohttp://127.0.0.1:26863/graphics/plot_zoom_png?width=980&height=712
  geom_hline(yintercept = max_value, color = "black", linetype = "dashed") +  # Linea nera al max valore
  geom_hline(yintercept = min_value, color = "black", linetype = "dashed") +  # Linea nera al min valore
  # Aggiungere le linee perpendicolari dai punti alla linea rossa
  geom_segment(aes(x = lag, xend = lag, y = weighted_value, yend = 0), color = "blue") +  # Linea perpendicolare da ogni punto a y = 0
  # Aggiungere le etichette per il massimo e il minimo
  geom_text(aes(x = lag[weighted_value == max_value], y = max_value, label = paste("Max:", round(max_value, 2))),
            color = "black", vjust = -1, size = 3) +  # Etichetta per il massimo
  geom_text(aes(x = lag[weighted_value == min_value], y = min_value, label = paste("Min:", round(min_value, 2))),
            color = "black", vjust = 1.5, size = 3) +  # Etichetta per il minimo
  theme_minimal() +  # Tema minimalista
  labs(title = "PACF: Mean Values by Lag",
       x = "Lag",
       y = "Mean PACF Value") +
  theme(axis.text.x = element_blank()) +  # Ruota le etichette sull'asse x
  ylim(-1, 1)  # Imposta i limiti dell'asse y da -1 a 1

file_path <- file.path(folder_name, paste0("df_pacf_combined.Rdata"))
save(result, file = file_path)



