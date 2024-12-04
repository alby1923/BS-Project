library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('df_filled_raw.rds')

#FUNCTIONS -----------
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
choose_columns <- function(df,response,threshold){
  folder_name <- paste0("STAN_", response)
  dir.create(folder_name)
  #Take only rows with no NA in a given column
  df_responses_filled <- subset_dataset(df,response)
  #compute percentage of NA for any covariate
  values <- plot_na_percentages(df_responses_filled)
  percentages <- values$na_percentages
  #ggsave(file.path(folder_name, paste0("NApercentage_", response, ".pdf")), values$plot, width = 10, height = 8, dpi = 300)
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
  cols_to_include <- setdiff(names(df_stan_obs), c("CAI", "Date", response,'delta_date','eta'))
  model_formula <- as.formula(paste(response, "~", paste(cols_to_include, collapse = " + ")))
  model <- lm(model_formula, data = df_stan_obs)
  coeffs <- summary(model)$coefficients
  chosen_columns <- rownames(coeffs)[coeffs[, "Pr(>|t|)"] < 0.1 & rownames(coeffs) != "(Intercept)"]
  df_stan_obs <-df_stan_obs[,c("CAI","Date",chosen_columns,response,'delta_date','eta')]
  return(df_stan_obs)
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

#MODELS ----------------
target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Trigliceridi","Glucosio","PMAX")
res_list <- list()
for (response in target_variables){
  folder_name <- paste0("STAN_", response)
  cat('\nStarting ',response,'\n')
  #create dataset and save it in the subfolder
  df_stan_obs <- choose_columns(df_stan,response,1)
  file_path <- file.path(folder_name, paste0("df_stan_obs_", response, "_new.rds"))
  saveRDS(df_stan_obs, file = file_path)
  #predictors + target
  chosen_columns <- setdiff(names(df_stan_obs), c("CAI", "Date", response,'delta_date','eta'))
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
  
  #extract results and obtain residuals
  posterior_samples <- extract(fit)
  y_hat <- posterior_samples$y_hat
  y_hat_point <- colMeans(y_hat, na.rm = TRUE) #matrix is number of iterations * number of observations
  #y_hat_point <- apply(y_hat, 2, median, na.rm = TRUE) #try with the median
  res <- t(target - y_hat_point)
  
  #group residuals to each patient
  patient_id <- df_stan_obs$CAI
  res_grouped <- split(res, patient_id)
  res_grouped_filtered <- res_grouped[sapply(res_grouped, length) >= 1]
  
  #associate timestamp to each time series
  res_grouped_with_filtered_timestamps <- lapply(names(res_grouped_filtered), function(CAI) {
    group <- res_grouped_filtered[[CAI]]
    relevant_timestamps <- unname(df_stan_obs$delta_date[df_stan_obs$CAI == CAI]/1) #if you want to create windows: %/%90 instead of /1
    list(
      results = group,         
      timestamps = relevant_timestamps  
    )
  })
  names(res_grouped_with_filtered_timestamps) <- names(res_grouped_filtered)
  res_list[[response]] <- res_grouped_with_filtered_timestamps
  file_path <- file.path(folder_name, paste0("residuals_", response, ".Rdata"))
  save(res_grouped_with_filtered_timestamps, file = file_path)
  
  #analyze performance
  rss <- sum(res^2)
  cat('\nRSS: ',rss)
  rmse <- sqrt(rss/dim(df_stan)[1]) 
  cat('\nRMSE: ',rmse)
  RMSE_to_sd_ratio <- rmse / sd(t(target))
  cat('\nRMSE/SD: ',RMSE_to_sd_ratio)
  
  # Posterior di sigma
  sigma_posterior <- posterior_samples$sigma
  
  # Plotta la distribuzione a posteriori di sigma e colorala di azzurro
  q <- stan_dens(fit, pars = "sigma_e")
  q + 
    ggtitle(paste("Sigma posterior distribution of ", response)) +  # Titolo personalizzato
    theme_minimal() +  # Usa un tema minimale per il grafico
    geom_density(fill = "skyblue", color = "blue")  # Colora la densità in azzurro
  plot(q)
  ggsave(paste0(folder_name, "/Sigma posterior distribution of ", response, ".png"), plot = q, width = 6, height = 4)
  
  # Plotta la distribuzione a posteriori dei beta e colorala di azzurro
  h <- stan_dens(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]")))
  h + 
    ggtitle(paste("Beta posterior distribution of ", response)) +  # Titolo personalizzato
    theme_minimal() +  # Tema minimale
    geom_density(fill = "skyblue", color = "blue")  # Colora la densità in azzurro
  plot(h)
  ggsave(paste0(folder_name, "/Beta posterior distribution of ", response, ".png"), plot = h, width = 6, height = 4)
  
  p <- ggplot(data.frame(t(res)), aes(x = t(res))) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1.5) +  # Aggiungi la curva di densità
    #stat_function(fun = dnorm, args = list(mean = mean(res), sd = sd(res)), color = "green", size = 1.5) +  # Aggiungi una curva normale teorica
    labs(title = paste("Residuals Histogram", response), x = "Valore", y = "Densità") +
    theme_minimal()
  plot(p)
  # Salva il grafico come immagine PNG
  ggsave(paste0(folder_name, "/Residuals histogram of ", response, ".png"), plot = p, width = 6, height = 4)
  
  print(summary(fit, pars = c("beta", "sigma_e")))
  q <- traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
  plot(q)
  ggsave(paste0(folder_name, "/Traceplots of ", response, ".png"), plot = q, width = 6, height = 4)
  #print(fit, pars = c("beta", "sigma_e"))
  
  # Calcola R² a posteriori
  target <- t(target)
  SS_tot <- sum((target - mean(target))^2)  
  SS_res <- apply(y_hat, 1, function(y_h) sum((target - y_h)^2))
  R2_posterior <- 1 - SS_res / SS_tot
  
  # Calcola la media di R² e l'intervallo di credibilità
  mean_R2 <- mean(R2_posterior)
  R2_ci <- quantile(R2_posterior, probs = c(0.025, 0.975))
  
  cat("\nMedia di R² a posteriori: ", round(mean_R2, 3), "\n")
  cat("Intervallo di credibilità del 95% per R²: ", round(R2_ci, 3), "\n")

  cat('\nEND\n---------------------------------\n')
  rm(fit)
  rm(posterior_samples)
}

#RESIDUALS ANALYSIS -----------
#Filter patients with at least threshold observations
create_residuals_dataset <- function(target_data, patient_info, threshold,target) {
  #list of results
  results <- list()
  total_steps <- sum(sapply(target_data, function(targets_data) length(targets_data)))
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  step_counter <- 0
    
    #Loops in CAI
    for (CAI in names(target_data)) {
      
      step_counter <- step_counter + 1
      setTxtProgressBar(pb, step_counter)
      
      patient_data <- target_data[[CAI]]
      residui <- patient_data$results
      delta_date <- patient_data$timestamps
      
      #Filter only patients with at least threshold observations
      if (length(residui) >= threshold) {
        patient_row <- patient_info[patient_info$CAI == CAI, ]
        if (nrow(patient_row) == 0) next  #just a check patient exists
        sesso <- patient_row$SESSO
        eta <- patient_row$eta
        
        #build dataframe for the patient
        patient_df <- data.frame(
          CAI = rep(CAI, length(residui)),
          residuo = residui,
          delta_date = delta_date,
          sesso = as.factor(rep(sesso, length(residui))),
          eta = eta
        )
        
        #add dataframe to list of results
        results[[paste(target, CAI, sep = "_")]] <- patient_df
      }
    }
  
  close(pb)
  #combine all dataframes
  final_dataset <- do.call(rbind, results)
  return(final_dataset)
}

target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Trigliceridi","Glucosio","PMAX")
res_list <- list()

for(response in target_variables){
  load(file.path(paste0("STAN_", response), paste0("residuals_", response, ".Rdata")))
  res_list[[response]] <- res_grouped_with_filtered_timestamps
}

residuals_dataset <- list()
for(response in target_variables){
cat('\nstarting: ', response)
df <- readRDS(paste0(file.path(paste0("STAN_", response), paste0("df_stan_obs_", response, "_new.rds"))))
patient_info <- df[, c("CAI", "SESSO", "eta")] %>% distinct(CAI, .keep_all = TRUE)
final_dataset <- create_residuals_dataset(res_list[[response]], patient_info, threshold = 10,response)
quantiles <- quantile(final_dataset$eta, probs = seq(0, 1, 0.25), na.rm = TRUE)
cat("\nQuantiles per age (", response, "):\n")
print(quantiles)
#add levels for ages based on quantiles
final_dataset$age_level <- cut(
  final_dataset$eta,
  breaks = quantiles,
  include.lowest = TRUE,
  labels = 1:4        
)
residuals_dataset[[response]] <- final_dataset
}

#add difference in time from previous observation
for (response in target_variables){
  df <- residuals_dataset[[response]]
  df <- df %>%
    group_by(CAI) %>% # Raggruppa per paziente
    arrange(delta_date, .by_group = TRUE) %>% # Ordina le osservazioni per delta_date
    mutate(
      successive_timestamps = if_else(
        delta_date == 0, # Se delta_date è 0
        0,               # Imposta a 0
        delta_date - lag(delta_date, default = 0) # Altrimenti calcola la differenza
      )
    ) %>%
    ungroup() # Rimuovi il raggruppamento
  residuals_dataset[[response]] <- df
}

working_directory <- getwd()
file_path <- file.path(working_directory, paste0("residuals_dataset.Rdata"))
save(residuals_dataset, file = file_path)

load('residuals_dataset.Rdata')

#build a total linear model
for(response in target_variables){
  cat('\n',response,'\n')
  df <- residuals_dataset[[response]]
  my_lm <- lm(residuo ~ delta_date + I(delta_date^2) + I(delta_date^3) + age_level*sesso,data = df)
  print(summary(my_lm))
  cat('--------------------------------------')
}

#binning version for residuals
# Loop sulle target variables
target_variables <- c("Colesterolo_Hdl","Trigliceridi","Glucosio","PMAX") #"Circonferenza_vita", too few observations
for (response in target_variables) {
  cat("\nTarget Variable: ", response, "\n")
  
  #extract dataset
  df <- residuals_dataset[[response]]
  
  #fixed intervals of bins 
  custom_breaks <- c(0, 180, 365, max(df$successive_timestamps))
  num_bins <- 3
  df$delta_date_bin <- cut(
    df$successive_timestamps,
    breaks = custom_breaks,
    include.lowest = TRUE,
    labels = paste0("Range_", seq_along(custom_breaks[-1]))
  )
  
  #all bins have the same amount of data
  #num_bins <- 4
  #df <- df %>%
  # mutate(delta_date_bin = ntile(delta_date, num_bins))
  #df$delta_date_bin <- paste0("Bin_", df$delta_date_bin)
  
  #print bin ranges
  cat("\nIntervalli dei bin per delta_date:\n")
  print(table(df$delta_date_bin))
  
  bin_ranges <- df %>%
    group_by(delta_date_bin) %>%
    summarise(
      min_delta_date = min(successive_timestamps, na.rm = TRUE),
      max_delta_date = max(successive_timestamps, na.rm = TRUE)
    )
  
  cat("\nIntervalli dei bin per delta_date:\n")
  print(bin_ranges)
  
  
  #linear model for each bin
  for (bin in unique(df$delta_date_bin)) {
    cat("\nBin: ", bin, "\n")
   #filter data for current bin
    bin_df <- df[df$delta_date_bin == bin, ]
    my_lm <- lm(residuo ~ successive_timestamps + I(successive_timestamps^2) + I(successive_timestamps^3) + age_level + sesso + age_level:sesso, data = bin_df)
    cat("\nSummary linear model for ", bin, ":\n")
    print(summary(my_lm))
  }
  cat('----------------------------')
}


library(BDgraph)

#mcmc filling (data simulation)
#residuals on filled dataset 
# (tentative) graph
#add new R2
#multivariate model 
