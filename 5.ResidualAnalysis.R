library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('df_filled_raw.rds')

#FUNCTIONS --------
subset_dataset <- function(df_wide,names){
  rows_to_remove = c()
  target_cols = which(colnames(df_wide) %in% names)
  pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
  for(i in 1:nrow(df_wide)){
    flag <- 0
    pb$tick()
    for(col in target_cols){
      if(col == 'Date'){
        next
      }
      if(flag == 1)
        break
      if(is.na(df_wide[i, col])){
        rows_to_remove <- c(rows_to_remove, i)
        flag <- 1
      }
    }
  }
  length(rows_to_remove)
  df <- df_wide[-rows_to_remove,]
  return(df)
}
plot_na_percentages <- function(dataset) {
  # Calcola la percentuale di NA per ogni colonna
  na_percentages <- sapply(dataset, function(col) mean(is.na(col)) * 100)
  
  # Crea un dataframe per ggplot
  na_data <- data.frame(
    Column = names(na_percentages),
    NA_Percentage = na_percentages
  )
  
  # Ordina per percentuale decrescente
  na_data <- na_data[order(-na_data$NA_Percentage), ]
  
  na_data$Highlight <- ifelse(na_data$NA_Percentage > 5, "Remove", "Keep")
  
  # Barplot
  plot <- ggplot(na_data, aes(x = reorder(Column, -NA_Percentage), y = NA_Percentage)) +
    geom_bar(
      stat = "identity",
      aes(fill = Highlight),
      width = .85,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("Remove" = "red", "Keep" = "steelblue")) +
    geom_text(
      aes(label = paste0(round(NA_Percentage, 1), "%")),
      hjust = -0.1,
      size = 3.5,
      color = "black"
    ) +
    coord_flip() +
    labs(
      x = "Columns",
      y = "Percentage of NA (%)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(size = 12, color = "darkgray"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()
    ) +
    ylim(0, max(na_data$NA_Percentage) + 5)  
  
  #print(plot)
  return(list(na_percentages = na_percentages, plot = plot))
}
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
choose_columns <- function(df,response,threshold){
  folder_name <- paste0("STAN_", response)
  dir.create(folder_name)
  #Ttake only rows with no NA in a given column
  df_responses_filled <- subset_dataset(df,response)
  #compute percentage of NA for ny covariate
  values <- plot_na_percentages(df_responses_filled)
  percentages <- values$na_percentages
  ggsave(file.path(folder_name, paste0("NApercentage_", response, ".pdf")), values$plot, width = 10, height = 8, dpi = 300)
  #we keep only covariates with maximum 10% of NA and take out other responses
  target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
  other_responses <- setdiff(target_variables,response)
  columns_to_keep <- setdiff(names(percentages[percentages <= 10]),other_responses)
  #remove rows with at least one NA in any of column to keep
  df_opt1 <- subset_dataset(df_responses_filled,columns_to_keep)
  df_opt1 <- df_opt1[,columns_to_keep]
  #just change name of df
  df_stan <- df_opt1
  #take only observations with at least a given number of observations
  df_stan_obs <- minimum_observations(df_stan,threshold)
  #build a linear model to discard non-significant covariates to speed up computations
  cols_to_include <- setdiff(names(df_stan_obs), c("CAI", "Date", response,'delta_date'))
  model_formula <- as.formula(paste(response, "~", paste(cols_to_include, collapse = " + ")))
  model <- lm(model_formula, data = df_stan_obs)
  coeffs <- summary(model)$coefficients
  chosen_columns <- rownames(coeffs)[coeffs[, "Pr(>|t|)"] < 0.1 & rownames(coeffs) != "(Intercept)"]
  df_stan_obs <-df_stan_obs[,c("CAI","Date",chosen_columns,response,'delta_date')]
  return(df_stan_obs)
}
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

#PLOT FOR OCCURRENCIES ------------
# Contare le occorrenze per ciascun CAI
cai_counts <- df_stan %>%
  group_by(CAI) %>%
  summarise(occurrences = n()) %>%
  arrange(desc(occurrences))  # Ordinare in ordine decrescente per comodità

# Plot delle occorrenze
ggplot(cai_counts, aes(x = reorder(CAI, -occurrences), y = occurrences)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_minimal() +
  labs(
    title = "Occorrenze di CAI nel dataset",
    x = "CAI",
    y = "Numero di occorrenze"
  ) +
  theme(axis.text.x = element_blank())  # Ruota le etichette per leggibilità

#BUILD MODEL -----------
# counts how many dates are before the reference date
reference_date <- as.Date("2019-01-01")
dates_before <- sum(df_stan$Date < reference_date, na.rm = TRUE)
dates_before #many variables started being tracked from 2019 on
#compute the delta time from first osbervation for each patient
df_stan <- df_stan[order(df_stan$CAI), ]
df_stan <- calculate_delta_date(df_stan, "CAI", "Date")
target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Trigliceridi","Glucosio","PMAX")
for (response in target_variables){
  folder_name <- paste0("STAN_", response)
  cat('\nStarting ',response,'\n')
  #create dataset and save it in the subfolder
  df_stan_obs <- choose_columns(df_stan,response,1)
  file_path <- file.path(folder_name, paste0("df_stan_obs_", response, "_new.rds"))
  saveRDS(df_stan_obs, file = file_path)
  #predictors + target
  chosen_columns <- setdiff(names(df_stan_obs), c("CAI", "Date", response,'delta_date'))
  target <- as.vector(df_stan_obs[[response]])
  #create data for model
  data_for_model <- data_stan(df_stan_obs,chosen_columns,target)
  #fit the model
  cat('\nsi parte\n')
  fit = stan(file = 'model_base.stan', 
             data = data_for_model, 
             chains = 2, 
             iter = 1000, 
             warmup = 500, 
             cores = 2,
             thin = 1, 
             seed = 19)
  #save fit
  file_path <- file.path(folder_name, paste0("fit_marginal_", response, "_new.Rdata"))
  save(fit, file = file_path)
  #save names of columns
  file_path <- file.path(folder_name, paste0("chosen_columns_", response, ".txt"))
  writeLines(chosen_columns, file_path)
  cat('\nEND RUNNING MCMC\n---------------------------------\n')
  rm(fit)
}

#sorry for the bad programming here, could have been done with a for loop instead of copy-pasting
#CIRCONFERENZA VITA -------------
my_df <- list()

df_vita <- readRDS('df_stan_obs_Circonferenza_vita_new.rds')
response <- df_vita$Circonferenza_vita
chosen_columns <- c(
  "Altezza",
  "Distribuzione_di_volume",
  "Ematocrito_hct",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_hb",
  "Eritrociti_rbc",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Peso",
  "Piastrine",
  "Polso",
  "Alcool",
  "Attivita_fisica",
  "Fumo",
  "Rh",
  "AB0",
  "SESSO"
)
#results of the STAN model
fit_vita <- load("fit_marginal_Circonferenza_vita_new.RData")
fit_vita <- fit
rm(fit)
traceplot(fit_vita, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
print(fit_vita, pars = c("beta", "sigma_e"))
summary(fit_vita, pars =  c("beta", "sigma_e"))

#extract results and obtain residuals
posterior_samples <- extract(fit_vita)
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE) #matrix is number of iterations * number of observations
#y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) #try with the median
res <- response - y_hat_point

#analyze performance
rss <- sum(res^2)
rmse <- sqrt(rss/dim(df_vita)[1]) #0.02085122
RMSE_to_sd_ratio <- rmse / sd(response) #0.1807194
hist(res)

# counts how many dates are before the reference date
reference_date <- as.Date("2019-01-01")
dates_before <- sum(df_vita$Date < reference_date, na.rm = TRUE)
dates_before #many variables started being tracked from 2019 on
#compute the delta time from first osbervation for each patient
df_vita <- df_vita[order(df_vita$CAI), ]
df_vita <- calculate_delta_date(df_vita, "CAI", "Date")

#group residuals to each patient and take only time series with at least 10 obs
patient_id <- df_vita$CAI
res_grouped <- split(res, patient_id)
res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 5]

#associate timestamp to each time series
res_grouped_with_filtered_timestamps <- lapply(names(res_grouped_filtered), function(CAI) {
  group <- res_grouped_filtered[[CAI]]
  relevant_timestamps <- unname(df_vita$delta_date[df_vita$CAI == CAI]%/%90)
  list(
    results = group,         
    timestamps = relevant_timestamps  
  )
})
names(res_grouped_with_filtered_timestamps) <- names(res_grouped_filtered)

#compute the fake PACF for irregular timestamps
df <- fake_pacf(res_grouped_with_filtered_timestamps)
df <- rbind(data.frame(lag = "0", mean_value = 1, count = 0), df)
#df$lag <- factor(df$lag, levels = c("lag0", paste0("lag", 1:max(as.numeric(gsub("lag", "", df$lag))))))

# Calcolare max e min value escluso 1 per le linee tratteggiate nere
max_value <- max(df$mean_value[df$mean_value != 1])
min_value <- min(df$mean_value[df$mean_value != 1])

# Crea un grafico a dispersione per visualizzare i valori PACF medi
ggplot(df, aes(x = lag, y = mean_value)) +
  geom_point(color = "blue", size = 3) +  # Punti per ogni lag
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Linea rossa a zero
  geom_hline(yintercept = max_value, color = "black", linetype = "dashed") +  # Linea nera al max valore
  geom_hline(yintercept = min_value, color = "black", linetype = "dashed") +  # Linea nera al min valore
  # Aggiungere le linee perpendicolari dai punti alla linea rossa
  geom_segment(aes(x = lag, xend = lag, y = mean_value, yend = 0), color = "blue") +  # Linea perpendicolare da ogni punto a y = 0
  # Aggiungere le etichette per il massimo e il minimo
  geom_text(aes(x = lag[mean_value == max_value], y = max_value, label = paste("Max:", round(max_value, 2))),
            color = "black", vjust = -1, size = 3) +  # Etichetta per il massimo
  geom_text(aes(x = lag[mean_value == min_value], y = min_value, label = paste("Min:", round(min_value, 2))),
            color = "black", vjust = 1.5, size = 3) +  # Etichetta per il minimo
  theme_minimal() +  # Tema minimalista
  labs(title = "PACF: Mean Values by Lag for Circonferenza Vita",
       x = "Lag",
       y = "Mean PACF Value") +
  theme(axis.text.x = element_blank()) +  # Ruota le etichette sull'asse x
  ylim(-1, 1)  # Imposta i limiti dell'asse y da -1 a 1

my_df[['Circonferenza_vita']] <- df
rm(fit_vita)
#COLESTEROLO ----------------
df_colesterolo <- readRDS('df_stan_obs_Colesterolo_Hdl_new.rds')
response <- df_colesterolo$Colesterolo_Hdl
chosen_columns <- c(
  "Colesterolo_totale",
  "Distribuzione_di_volume",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_hb",
  "Eritrociti_rbc",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Peso",
  "Piastrine",
  "Polso",
  "Rh",
  "SESSO"
)

#results of the STAN model
fit_colesterolo <- load("fit_marginal_Colesterolo_Hdl_new.RData")
fit_colesterolo <- fit
rm(fit)
traceplot(fit_colesterolo, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
print(fit_colesterolo, pars =  c("beta", "sigma_e"))
summary(fit_colesterolo, pars =  c("beta", "sigma_e"))

#extract results and obtain residuals
posterior_samples <- extract(fit_colesterolo)
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE) #matrix is number of iterations * number of observations
#y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) #try with the median
res <- response - y_hat_point

#analyze performance
rss <- sum(res^2)
rmse <- sqrt(rss/dim(df_colesterolo)[1]) #0.088
RMSE_to_sd_ratio <- rmse / sd(response) #0.363
hist(res)

# counts how many dates are before the reference date
reference_date <- as.Date("2019-01-01")
dates_before <- sum(df_colesterolo$Date < reference_date, na.rm = TRUE)
dates_before #many variables started being tracked from 2019 on
#compute the delta time from first osbervation for each patient
df_colesterolo <- df_colesterolo[order(df_colesterolo$CAI), ]
df_colesterolo <- calculate_delta_date(df_colesterolo, "CAI", "Date")

#group residuals to each patient and take only time series with at least 10 obs
patient_id <- df_colesterolo$CAI
res_grouped <- split(res, patient_id)
res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 10]

#associate timestamp to each time series
res_grouped_with_filtered_timestamps <- lapply(names(res_grouped_filtered), function(CAI) {
  group <- res_grouped_filtered[[CAI]]
  relevant_timestamps <- unname(df_colesterolo$delta_date[df_colesterolo$CAI == CAI]%/%90)
  list(
    results = group,         
    timestamps = relevant_timestamps  
  )
})
names(res_grouped_with_filtered_timestamps) <- names(res_grouped_filtered)

df_col <- fake_pacf(res_grouped_with_filtered_timestamps)
df_col <- rbind(data.frame(lag = "0", mean_value = 1, count = 0), df_col)
#df$lag <- factor(df$lag, levels = c("lag0", paste0("lag", 1:max(as.numeric(gsub("lag", "", df$lag))))))

# Calcolare max e min value escluso 1 per le linee tratteggiate nere
max_value <- max(df_col$mean_value[df_col$mean_value != 1])
min_value <- min(df_col$mean_value[df_col$mean_value != 1])

# Crea un grafico a dispersione per visualizzare i valori PACF medi
ggplot(df_col, aes(x = lag, y = mean_value)) +
  geom_point(color = "blue", size = 3) +  # Punti per ogni lag
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Linea rossa a zero
  geom_hline(yintercept = max_value, color = "black", linetype = "dashed") +  # Linea nera al max valore
  geom_hline(yintercept = min_value, color = "black", linetype = "dashed") +  # Linea nera al min valore
  # Aggiungere le linee perpendicolari dai punti alla linea rossa
  geom_segment(aes(x = lag, xend = lag, y = mean_value, yend = 0), color = "blue") +  # Linea perpendicolare da ogni punto a y = 0
  # Aggiungere le etichette per il massimo e il minimo
  geom_text(aes(x = lag[mean_value == max_value], y = max_value, label = paste("Max:", round(max_value, 2))),
            color = "black", vjust = -1, size = 3) +  # Etichetta per il massimo
  geom_text(aes(x = lag[mean_value == min_value], y = min_value, label = paste("Min:", round(min_value, 2))),
            color = "black", vjust = 1.5, size = 3) +  # Etichetta per il minimo
  theme_minimal() +  # Tema minimalista
  labs(title = "PACF: Mean Values by Lag for Colesterolo",
       x = "Lag",
       y = "Mean PACF Value") +
  theme(axis.text.x = element_blank()) +  # Ruota le etichette sull'asse x
  ylim(-1, 1)  # Imposta i limiti dell'asse y da -1 a 1

my_df[['Colesterolo_Hdl']] <- df_col
rm(fit_colesterolo)
#TRIGLICERIDI --------------
df_trigliceridi <- readRDS('df_stan_obs_Trigliceridi_new.rds')
response <- df_trigliceridi$Trigliceridi
chosen_columns <- c(
  "Altezza",
  "Colesterolo_totale",
  "Distribuzione_di_volume",
  "Eosinofili_perc",
  "Eritrociti_rbc",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Peso",
  "Rh",
  "AB0",
  "SESSO"
)

#results of the STAN model
fit_trigliceridi <- load("fit_marginal_Trigliceridi_new.RData")
fit_trigliceridi <- fit
rm(fit)
traceplot(fit_trigliceridi, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
print(fit_trigliceridi, pars = c('beta', "sigma_e"))
summary(fit_trigliceridi, pars = c('beta', "sigma_e"))

#extract results and obtain residuals
posterior_samples <- extract(fit_trigliceridi)
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE) #matrix is number of iterations * number of observations
#y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) #try with the median
res <- response - y_hat_point

#analyze performance
rss <- sum(res^2)
rmse <- sqrt(rss/dim(df_trigliceridi)[1]) #0.2405322
RMSE_to_sd_ratio <- rmse / sd(response) #0.5189739
hist(res)

# counts how many dates are before the reference date
reference_date <- as.Date("2019-01-01")
dates_before <- sum(df_trigliceridi$Date < reference_date, na.rm = TRUE)
dates_before #many variables started being tracked from 2019 on
#compute the delta time from first osbervation for each patient
df_trigliceridi <- df_trigliceridi[order(df_trigliceridi$CAI), ]
df_trigliceridi <- calculate_delta_date(df_trigliceridi, "CAI", "Date")

#group residuals to each patient and take only time series with at least 10 obs
patient_id <- df_trigliceridi$CAI
res_grouped <- split(res, patient_id)
res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 10]

#associate timestamp to each time series
res_grouped_with_filtered_timestamps <- lapply(names(res_grouped_filtered), function(CAI) {
  group <- res_grouped_filtered[[CAI]]
  relevant_timestamps <- unname(df_trigliceridi$delta_date[df_trigliceridi$CAI == CAI]%/%90)
  list(
    results = group,         
    timestamps = relevant_timestamps  
  )
})
names(res_grouped_with_filtered_timestamps) <- names(res_grouped_filtered)

#compute the fake PACF for irregular timestamps
df_trig <- fake_pacf(res_grouped_with_filtered_timestamps)
df_trig <- rbind(data.frame(lag = "0", mean_value = 1, count = 0), df_trig)
#df$lag <- factor(df$lag, levels = c("lag0", paste0("lag", 1:max(as.numeric(gsub("lag", "", df$lag))))))

# Calcolare max e min value escluso 1 per le linee tratteggiate nere
max_value <- max(df_trig$mean_value[df_trig$mean_value != 1])
min_value <- min(df_trig$mean_value[df_trig$mean_value != 1])

# Crea un grafico a dispersione per visualizzare i valori PACF medi
ggplot(df_trig, aes(x = lag, y = mean_value)) +
  geom_point(color = "blue", size = 3) +  # Punti per ogni lag
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Linea rossa a zero
  geom_hline(yintercept = max_value, color = "black", linetype = "dashed") +  # Linea nera al max valore
  geom_hline(yintercept = min_value, color = "black", linetype = "dashed") +  # Linea nera al min valore
  # Aggiungere le linee perpendicolari dai punti alla linea rossa
  geom_segment(aes(x = lag, xend = lag, y = mean_value, yend = 0), color = "blue") +  # Linea perpendicolare da ogni punto a y = 0
  # Aggiungere le etichette per il massimo e il minimo
  geom_text(aes(x = lag[mean_value == max_value], y = max_value, label = paste("Max:", round(max_value, 2))),
            color = "black", vjust = -1, size = 3) +  # Etichetta per il massimo
  geom_text(aes(x = lag[mean_value == min_value], y = min_value, label = paste("Min:", round(min_value, 2))),
            color = "black", vjust = 1.5, size = 3) +  # Etichetta per il minimo
  theme_minimal() +  # Tema minimalista
  labs(title = "PACF: Mean Values by Lag for Trigliceridi",
       x = "Lag",
       y = "Mean PACF Value") +
  theme(axis.text.x = element_blank()) +  # Ruota le etichette sull'asse x
  ylim(-1, 1)  # Imposta i limiti dell'asse y da -1 a 1

my_df[['Trigliceridi']] <- df_trig
rm(fit_trigliceridi)
#GLUCOSIO -------------
df_glucosio <- readRDS('df_stan_obs_Glucosio_new.rds')
response <- df_glucosio$Glucosio
chosen_columns <- c(
  "Altezza",
  "Distribuzione_di_volume",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_massa_media_mch",
  "Eosinofili_perc",
  "Eritrociti_rbc",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Monociti_perc",
  "Peso",
  "Piastrine",
  "Polso",
  "Volume_medio",
  "AB0",
  "SESSO"
)

# Results of the STAN model
fit_glucosio <- load("fit_marginal_Glucosio_new.RData")
fit_glucosio <- fit
rm(fit)
traceplot(fit_glucosio, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
print(fit_glucosio, pars = c("beta", "sigma_e"))
summary(fit_glucosio, pars = c("beta", "sigma_e"))

# Extract results and obtain residuals
posterior_samples <- extract(fit_glucosio)
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE) # Matrix is number of iterations * number of observations
# y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) # Try with the median
res <- response - y_hat_point

# Analyze performance
rss <- sum(res^2)
rmse <- sqrt(rss/dim(df_glucosio)[1]) # 0.093
RMSE_to_sd_ratio <- rmse / sd(response) # 0.79
hist(res)

# Counts how many dates are before the reference date
reference_date <- as.Date("2019-01-01")
dates_before <- sum(df_glucosio$Date < reference_date, na.rm = TRUE)
dates_before # Many variables started being tracked from 2019 on
# Compute the delta time from first observation for each patient
df_glucosio <- df_glucosio[order(df_glucosio$CAI), ]
df_glucosio <- calculate_delta_date(df_glucosio, "CAI", "Date")

# Group residuals to each patient and take only time series with at least 10 obs
patient_id <- df_glucosio$CAI
res_grouped <- split(res, patient_id)
res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 10]

# Associate timestamp to each time series
res_grouped_with_filtered_timestamps <- lapply(names(res_grouped_filtered), function(CAI) {
  group <- res_grouped_filtered[[CAI]]
  relevant_timestamps <- unname(df_glucosio$delta_date[df_glucosio$CAI == CAI] %/% 90)
  list(
    results = group,         
    timestamps = relevant_timestamps  
  )
})
names(res_grouped_with_filtered_timestamps) <- names(res_grouped_filtered)

# Compute the fake PACF for irregular timestamps
df <- fake_pacf(res_grouped_with_filtered_timestamps)
df <- rbind(data.frame(lag = "0", mean_value = 1, count = 0), df)

# Calcolare max e min value escluso 1 per le linee tratteggiate nere
max_value <- max(df$mean_value[df$mean_value != 1])
min_value <- min(df$mean_value[df$mean_value != 1])

# Create a scatter plot to visualize the mean PACF values
ggplot(df, aes(x = lag, y = mean_value)) +
  geom_point(color = "blue", size = 3) +  # Points for each lag
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Red line at zero
  geom_hline(yintercept = max_value, color = "black", linetype = "dashed") +  # Black line at max value
  geom_hline(yintercept = min_value, color = "black", linetype = "dashed") +  # Black line at min value
  # Add perpendicular lines from points to the red dashed line
  geom_segment(aes(x = lag, xend = lag, y = mean_value, yend = 0), color = "blue") +  # Perpendicular line from each point to y = 0
  # Add labels for max and min
  geom_text(aes(x = lag[mean_value == max_value], y = max_value, label = paste("Max:", round(max_value, 2))),
            color = "black", vjust = -1, size = 3) +  # Label for max
  geom_text(aes(x = lag[mean_value == min_value], y = min_value, label = paste("Min:", round(min_value, 2))),
            color = "black", vjust = 1.5, size = 3) +  # Label for min
  theme_minimal() +  # Minimal theme
  labs(title = "PACF: Mean Values by Lag for Glucosio",
       x = "Lag",
       y = "Mean PACF Value") +
  theme(axis.text.x = element_blank()) +  # Remove x-axis text
  ylim(-1, 1)  # Set y-axis limits to -1 and 1

my_df[['Glucosio']] <- df
rm(fit_glucosio)
#PMAX -----------
df_p <- readRDS('df_stan_obs_PMAX_new.rds')
response <- df_p$PMAX
chosen_columns <- c(
  "Altezza",
  "Distribuzione_di_volume",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_massa_media_mch",
  "Eosinofili_perc",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Monociti_perc",
  "Peso",
  "Piastrine",
  "Polso",
  "Volume_medio",
  "SESSO"
)

# Results of the STAN model
fit_p <- load("fit_marginal_PMAX_new.RData")
fit_p <- fit
rm(fit)
traceplot(fit_p, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
#print(fit_p, pars = c("beta", "sigma_e"))
#summary(fit_p, pars = c("beta", "sigma_e"))

# Extract results and obtain residuals
posterior_samples <- extract(fit_p)
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE) # Matrix is number of iterations * number of observations
# y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) # Try with the median
res <- response - y_hat_point

# Analyze performance
rss <- sum(res^2)
rmse <- sqrt(rss/dim(df_p)[1]) # 0.06
RMSE_to_sd_ratio <- rmse / sd(response) # 0.52
hist(res)

# Counts how many dates are before the reference date
reference_date <- as.Date("2019-01-01")
dates_before <- sum(df_p$Date < reference_date, na.rm = TRUE)
dates_before # Many variables started being tracked from 2019 on
# Compute the delta time from first observation for each patient
df_p <- df_p[order(df_p$CAI), ]
df_p <- calculate_delta_date(df_p, "CAI", "Date")

# Group residuals to each patient and take only time series with at least 10 obs
patient_id <- df_p$CAI
res_grouped <- split(res, patient_id)
res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 10]

# Associate timestamp to each time series
res_grouped_with_filtered_timestamps <- lapply(names(res_grouped_filtered), function(CAI) {
  group <- res_grouped_filtered[[CAI]]
  relevant_timestamps <- unname(df_p$delta_date[df_p$CAI == CAI] %/% 90)
  list(
    results = group,         
    timestamps = relevant_timestamps  
  )
})
names(res_grouped_with_filtered_timestamps) <- names(res_grouped_filtered)

# Compute the fake PACF for irregular timestamps
df <- fake_pacf(res_grouped_with_filtered_timestamps)
df <- rbind(data.frame(lag = "0", mean_value = 1, count = 0), df)

# Calcolare max e min value escluso 1 per le linee tratteggiate nere
max_value <- max(df$mean_value[df$mean_value != 1])
min_value <- min(df$mean_value[df$mean_value != 1])

# Create a scatter plot to visualize the mean PACF values
ggplot(df, aes(x = lag, y = mean_value)) +
  geom_point(color = "blue", size = 3) +  # Points for each lag
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Red line at zero
  geom_hline(yintercept = max_value, color = "black", linetype = "dashed") +  # Black line at max value
  geom_hline(yintercept = min_value, color = "black", linetype = "dashed") +  # Black line at min value
  # Add perpendicular lines from points to the red dashed line
  geom_segment(aes(x = lag, xend = lag, y = mean_value, yend = 0), color = "blue") +  # Perpendicular line from each point to y = 0
  # Add labels for max and min
  geom_text(aes(x = lag[mean_value == max_value], y = max_value, label = paste("Max:", round(max_value, 2))),
            color = "black", vjust = -1, size = 3) +  # Label for max
  geom_text(aes(x = lag[mean_value == min_value], y = min_value, label = paste("Min:", round(min_value, 2))),
            color = "black", vjust = 1.5, size = 3) +  # Label for min
  theme_minimal() +  # Minimal theme
  labs(title = "PACF: Mean Values by Lag for PMAX",
       x = "Lag",
       y = "Mean PACF Value") +
  theme(axis.text.x = element_blank()) +  # Remove x-axis text
  ylim(-1, 1)  # Set y-axis limits to -1 and 1

my_df[['PMAX']] <- df
rm(fit_p)

