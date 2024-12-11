library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

#load data
df_stan <- readRDS('to_simulate.rds')

#NA VISUALIZATION --------------
#visualization of how many NA
plot_na_counts <- function(df, variables, ylim_max) {
  # Create dataset with NA counts
  na_counts <- data.frame(
    Variable = variables,
    NA_Count = sapply(variables, function(var) sum(is.na(df[[var]])))
  )
  
  # Barplot
  plot <- ggplot(na_counts, aes(x = reorder(Variable, NA_Count), y = NA_Count)) +
    geom_bar(stat = "identity", aes(fill = NA_Count), width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = NA_Count), vjust = -0.5, size = 4, color = "black") +
    scale_fill_gradient(low = "skyblue", high = "darkblue") +
    labs(
      x = "Variable",
      y = "Count of NA"
    ) +
    ylim(0, ylim_max) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(),  # Rimuove le linee di griglia principali
      panel.grid.minor = element_blank(),  # Rimuove le linee di griglia minori
      panel.background = element_blank(),  # Imposta lo sfondo del pannello su bianco
      plot.background = element_blank()    # Imposta lo sfondo del grafico su bianco
    )
  
  ggsave("NAresponses.pdf", plot, width = 10, height = 8, dpi = 300)
  # Print the plot
  print(plot)
}

target_variables <- c("Glucosio","PMAX","Colesterolo_Hdl","Trigliceridi","Circonferenza_vita")

plot_na_counts(
  df = df_stan, 
  variables = target_variables, 
  ylim_max = 90000
)
#remove rows with NA in PMAX and Glucosio since they are just a few (meaningless to run a 3 hours simulation for less than 4000 observations out 90k)
df_stan <- df_stan[!is.na(df_stan$PMAX) & !is.na(df_stan$Glucosio), ]
#DATA SIMULATION ----------------
data_stan <- function(df, variables_to_keep, response) {

  patients_mat <- as.matrix(unique(df[,1]))
  patients <- as.vector(patients_mat)
  N <- length(patients)
  T <- nrow(df)
  
  covariates <- df[, variables_to_keep]
  
  tot_obs_patient <- as.vector(table(df$CAI)[match(patients_mat, names(table(df$CAI)))])
  subj <- rep(seq_along(tot_obs_patient), times = tot_obs_patient)
  
  observed_indices <- which(!is.na(response))
  missing_indices <- which(is.na(response))
  y_obs <- response[observed_indices]
  
  X1 <- as.matrix(covariates[observed_indices, ])
  X2 <- as.matrix(covariates[missing_indices, ])
  
  data_model <- list(
    N = N, 
    P = ncol(covariates), 
    T = T, 
    subj = subj, 
    X1 = X1, 
    X2 = X2,
    n_obs = length(observed_indices),
    n_miss = length(missing_indices),
    obs_indices = observed_indices,
    miss_indices = missing_indices,
    y_obs = y_obs
  )
  
  return(data_model)
}

#cut responses
chosen_columns <- c(
  "Altezza",
  #"Colesterolo_Hdl",
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
  #"Circonferenza_vita",
  "Rh",
  "SESSO",
  "AB0",
  'eta_std'
)

target_variables <- c("Circonferenza_vita","Colesterolo_Hdl","Trigliceridi")
#function to compute bayesian R2
bayes_R2 <- function(posterior_samples) {
  y_pred <- posterior_samples$y_obs_hat
  var_fit <- apply(y_pred, 1, var)
  var_res <- posterior_samples$sigma_e
  var_fit / (var_fit + var_res)
}
df_stan_filtered <- df_stan
rm(df_stan)

#remove donor who dont have any true value in a covariate (you cant estimate the patient-specific intercept)
for (target in target_variables){
df_stan_filtered <- df_stan_filtered %>%
  group_by(CAI) %>%
  filter(any(!is.na(.data[[target]]))) %>%
  ungroup()
}

target_variables <- c("Circonferenza_vita","Colesterolo_Hdl","Trigliceridi") #
for (target in target_variables){
  
  cat('\n Starting: ',target)
  
  dir.create(paste0('SIMULATED_',target))
  folder_name <- paste0('SIMULATED_',target)
  
  response <- as.vector(t(df_stan_filtered[,target]))
  data_for_model <- data_stan(df_stan_filtered,chosen_columns,response)
  
  #change values for simulation
  fit = stan(file = 'model_simulate.stan', 
             data = data_for_model, 
             chains = 2, 
             iter = 1000, 
             warmup = 500, 
             cores = 2,
             thin = 1,
             seed = 19)
  
  p <- traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
  ggsave(paste0(folder_name, "/Traceplots of ", target, ".png"), plot = p, width = 6, height = 4)
  
  q <- stan_dens(fit, pars = "sigma_e")
  ggsave(paste0(folder_name, "/Sigma posterior distribution of ", target, ".png"), plot = q, width = 6, height = 4)
  
  file_path <- file.path(folder_name, paste0("fit_",target,"_simulated.Rdata"))
  save(fit, file = file_path)
  
  summary_stats <- summary(fit, pars = c("beta", "sigma_e"))
  summary_text <- capture.output(print(summary_stats))
  writeLines(summary_text, con = paste0('summary ',target))
  
  posterior_samples <- extract(fit)
  rm(fit)
  
  r2 <- bayes_R2(posterior_samples)
  k <- ggplot(data.frame(r2 = r2), aes(x = r2)) +
    geom_histogram(binwidth = 0.001, fill = "blue", color = "black", alpha = 0.7) +
    labs(
      title = paste("Histogram of ", target),
      x = target,
      y = "Density"
    ) +
    theme_minimal()
  ggsave(paste0(folder_name, "/Histogram of Bayesian R2 of ", target, ".png"), plot = k, width = 6, height = 4)
  #print(mean(r2))
  
  simulated <- colMeans(posterior_samples$y_miss_hat)
  missing_indices <- which(is.na(df_stan_filtered[,target]))
  if (length(missing_indices) == length(simulated)) {
    for (i in 1:length(missing_indices)) {
      df_stan_filtered[missing_indices[i], target] <- simulated[i]
    }
  } else {
    stop("\nMismatch between missing indices and y_miss_hat length!")
  }
  cat('\n---------------------------------------\n')
  
  rm(posterior_samples)
  rm(p)
  rm(q)
  rm(k)
  rm(summary_text)
  rm(summary_stats)
}

saveRDS(df_stan_filtered,'ultimate_dataset.rds')
