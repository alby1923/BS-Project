library(rstan)
library(dplyr)
library(ggplot2)

df_stan <- readRDS('ultimate_dataset.rds')

data_stan <- function(df,variables_to_keep,target_variables,lags){
  
  df_trans_list <- list()
  for (i in seq_along(lags)) {
    df1 <- df %>%
      group_by(CAI) %>%
      arrange(delta_date, .by_group = TRUE) %>%
      mutate(
        successive_timestamps = if_else(
          delta_date == 0,
          0,               
          delta_date - lag(delta_date, default = 0)
        ),
        ar_flag = if_else(
          successive_timestamps <= lags[i] & successive_timestamps != 0, #if zero it is the first observation of the donor
          1, 
          0 #if less than zero it means that in delta_date - lag(delta_date, default = 0) it took the last observation of the previous donor
        ))%>%
      ungroup()  
    df_trans_list[[i]] <- df1
  }
  
  y_trasl <- matrix(0, nrow = length(target_variables), ncol = dim(df)[1])
  
  k = 0
  for (target in target_variables){
    k = k+1
    df2 <- df_trans_list[[k]]
    response <- df[[target]]
    aux <- df2[['ar_flag']] * lag(response, default = 0)
    for(j in 1:length(aux)){
      y_trasl[k,j] <- aux[j]
    }
  }
  
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
  
  y <- df[,target_variables]
  
  data_model = list(N = N, 
                    P = dim(covariates)[2], 
                    T = T, 
                    subj = subj,
                    K = length(target_variables),
                    y = y,
                    y_trasl = t(y_trasl),
                    X = covariates)
  
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

target_variables <- c("Glucosio","Colesterolo_Hdl","PMAX","Trigliceridi") #PUT THEM IN ORDER OF DEPENDANCY! "Circonferenza_vita"

df_stan_filtrato <- df_stan %>%
  group_by(CAI) %>%                  # Raggruppa per paziente
  filter(n() >= 20) %>%             # Mantieni solo i gruppi con almeno 20 osservazioni
  ungroup()                         # Rimuovi il raggruppamento

lags <- c(0,180,0,0) #max window allowed for previous observation in autoregressive component, if 0 it means no autoregressive component
data_model <- data_stan(df_stan_filtrato,chosen_columns,target_variables,lags)

fit = stan(file = 'model_final.stan', 
           data = data_model, 
           chains = 2, 
           iter = 1000, 
           warmup = 500, 
           cores = 2,
           thin = 1,
           seed = 19)

#TRACEPLOTS --------------
dir.create('STAN_Final')
folder_name <- 'STAN_Final'
K <- 4 #number of target variables
P <- length(chosen_columns)

for (k in 1:K) {
  # Traceplot di beta e sigma_e per ogni target
  beta_param <- grep(paste0("^beta\\[.*,", k, "\\]"), names(fit@sim$samples[[1]]), value = TRUE)
  traceplot_beta_sigma_e <- traceplot(fit, pars = c(beta_param, paste0('sigma_e[', k, ']')))
  ggsave(paste0(folder_name, "/Traceplot_beta_sigma_e_target_", target_variables[k], ".png"), plot = traceplot_beta_sigma_e)
  
  # Traceplot per phi (intercetta target-specifica)
  traceplot_phi <- traceplot(fit, pars = paste0("phi[", k, "]"))
  ggsave(paste0(folder_name, "/Traceplot_phi_target_", target_variables[k], ".png"), plot = traceplot_phi)
  
  # Traceplot per autoregressive coefficient per target (solo per target specifico)
  traceplot_autoreg_coef <- traceplot(fit, pars = paste0("autoreg_coef[", k, "]"))
  ggsave(paste0(folder_name, "/Traceplot_autoreg_coef_target_", target_variables[k], ".png"), plot = traceplot_autoreg_coef)
  
  # Traceplot di gamma (elementi fuori-diagonale) per ogni target
  params_gamma <- c()
  if(k>1){
    for (i in 1:(k - 1)) {
      params_gamma <- c(params_gamma, paste0("gamma[", i, ",", k, "]"))
    }
    traceplot_gamma <- traceplot(fit, pars = params_gamma)
    ggsave(paste0(folder_name, "/Traceplot_gamma_target_", target_variables[k], ".png"), plot = traceplot_gamma)
  }
}

# Traceplot di alpha per il primo e il secondo paziente
traceplot_alpha <- traceplot(fit, pars = c("alpha[1]","alpha[2]","alpha[3]","alpha[4]"))
ggsave(paste0(folder_name, "/Traceplot_alpha_patients.png"), plot = traceplot_alpha)

#RESULTS ANALYSIS -------------------
bayes_R2 <- function(posterior_samples, target_index) {
  y_pred <- posterior_samples$y_hat[, , target_index]
  var_fit <- apply(y_pred, 1, var)
  var_res <- posterior_samples$sigma_e[, target_index]
  var_fit / (var_fit + var_res)
}

posterior_samples <- extract(fit)

for (k in 1:K) {
  target_name <- target_variables[k]
  r2 <- bayes_R2(posterior_samples, target_index = k)
  r2_plot <- ggplot(data.frame(r2 = r2), aes(x = r2)) +
    geom_histogram(binwidth = 0.0001, fill = "blue", color = "black", alpha = 0.7) +
    labs(
      title = paste("Histogram of R2 for ", target_name),
      x = paste("R2 (", target_name, ")"),
      y = "Density"
    ) +
    theme_minimal()
  ggsave(paste0(folder_name, "/Bayesian R2 of ", target_name, ".png"), plot = r2_plot, width = 6, height = 4)
  
  y_hat <- posterior_samples$y_hat[, , k]
  y_hat_point <- colMeans(y_hat, na.rm = TRUE)
  res <- t(df_stan[, target_name] - y_hat_point)
  rss <- sum(res^2)
  cat("\nRSS for", target_name, ": ", rss)
  rmse <- sqrt(rss / nrow(df_stan))
  cat("\nRMSE for", target_name, ": ", rmse)
  RMSE_to_sd_ratio <- rmse / sd(t(df_stan[, target_name]))
  cat("\nRMSE/SD for", target_name, ": ", RMSE_to_sd_ratio)
  
  sigma_plot <- stan_dens(fit, pars = paste0("sigma_e[", k, "]")) +
    labs(
      title = paste("Posterior Distribution of Sigma for ", target_name),
      x = paste("Sigma (", target_name, ")"),
      y = "Density"
    )
  ggsave(paste0(folder_name, "/Sigma posterior distribution of ", target_name, ".png"), plot = sigma_plot, width = 6, height = 4)
  
  gamma_param <- c()
  if(k>1){
  for (i in 1:(k-1)) {
    gamma_param <- c(gamma_param, paste0("gamma[", i, ",", k, "]"))
  }
  }
  
  summary_stats <- summary(fit, pars = c(paste0("beta[", 1:P, ",", k, "]") ,
                                         paste0("sigma_e[", k, "]"), 
                                         paste0("phi[", k, "]"),
                                         paste0("autoreg_coef[", k, "]"),
                                         gamma_param))
  
  summary_text <- capture.output(print(summary_stats))
  writeLines(summary_text, con = paste0(folder_name, "/Summary of ", target_name, ".txt"))
  cat('\n--------------------------\n')
}

#RISK FACTORS -----------------

target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
df_stan_clusters <- df_stan

for (target in target_variables) {
  df_stan_clusters[[target]] <- exp(df_stan_clusters[[target]])
}

# Definizione delle soglie per ogni variabile
glucosio_threshold <- 100
trigliceridi_threshold <- 150
pmax_threshold <- 130
circonferenza_vita_maschi_threshold <- 102
circonferenza_vita_femmine_threshold <- 88
colesterolo_maschi_threshold <- 40
colesterolo_femmine_threshold <- 50

# Creazione della variabile risk_factor
df_stan_clusters$risk_factor <- NA  # inizializza la variabile

# Funzione per calcolare il rischio
for (i in 1:nrow(df_stan_clusters)) {
  # Estrai i valori per ogni osservazione
  glucosio <- df_stan_clusters$Glucosio[i]
  trigliceridi <- df_stan_clusters$Trigliceridi[i]
  pmax <- df_stan_clusters$PMAX[i]
  circonferenza_vita <- df_stan_clusters$Circonferenza_vita[i]
  colesterolo <- df_stan_clusters$Colesterolo_Hdl[i]
  sesso <- df_stan_clusters$SESSO[i]
  
  # Determina il sesso: se negativo, è maschio, altrimenti femmina
  if (sesso < 0) {
    sesso <- "maschio"
    circonferenza_threshold <- circonferenza_vita_maschi_threshold
    colesterolo_threshold <- colesterolo_maschi_threshold
  } else {
    sesso <- "femmina"
    circonferenza_threshold <- circonferenza_vita_femmine_threshold
    colesterolo_threshold <- colesterolo_femmine_threshold
  }
  
  # Verifica le condizioni per ciascuna variabile
  risk_count <- 0
  
  # Condizioni per glucosio
  if (glucosio > glucosio_threshold) risk_count <- risk_count + 1
  
  # Condizioni per trigliceridi
  if (trigliceridi > trigliceridi_threshold) risk_count <- risk_count + 1
  
  # Condizioni per pmax
  if (pmax > pmax_threshold) risk_count <- risk_count + 1
  
  # Condizioni per circonferenza_vita
  if (circonferenza_vita > circonferenza_threshold) risk_count <- risk_count + 1
  
  # Condizioni per colesterolo
  if (colesterolo < colesterolo_threshold) risk_count <- risk_count + 1
  
  # Assegna il rischio in base al numero di condizioni soddisfatte
  if (risk_count >= 3) {
    df_stan_clusters$risk_factor[i] <- "HIGH"
  } else if (risk_count == 2) {
    df_stan_clusters$risk_factor[i] <- "MEDIUM"
  } else {
    df_stan_clusters$risk_factor[i] <- "LOW"
  }
}

# Carica ggplot2 se non è già caricato
library(ggplot2)

# Creiamo un grafico a barre che mostra la distribuzione della variabile risk_factor
ggplot(df_stan_clusters, aes(x = risk_factor, fill = risk_factor)) +
  geom_bar(stat = "count", show.legend = FALSE) +  # Conta le osservazioni per ogni categoria
  scale_fill_manual(values = c("HIGH" = "red", "MEDIUM" = "yellow", "LOW" = "green")) +  # Colori personalizzati per ogni categoria
  labs(
    title = "Distribuzione del fattore di rischio",  # Titolo del grafico
    x = "Fattore di rischio",  # Etichetta asse X
    y = "Numero di osservazioni"  # Etichetta asse Y
  ) +
  theme_minimal() +  # Tema del grafico
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotazione delle etichette sull'asse X
    plot.title = element_text(hjust = 0.5),  # Centra il titolo
    axis.title = element_text(size = 12),  # Imposta la dimensione del titolo dell'asse
    axis.text = element_text(size = 10)  # Imposta la dimensione del testo sugli assi
  )

#CLUSTERS --------

data_for_model_clusters <- list(
  N = nrow(df_stan),
  D = length(chosen_columns),
  X = df_stan[,chosen_columns]
)

fit = stan(file = 'model_clusters.stan', 
           data = data_for_model_clusters, 
           chains = 2, 
           iter = 1000, 
           warmup = 500, 
           cores = 2,
           thin = 1,
           seed = 19)
