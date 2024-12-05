library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('ultimate_dataset.rds')

dir.create('STAN_...') #add your target instead of ...
folder_name <- paste0("STAN_...") #add your target instead of ...

#STAN MODEL -------
data_stan <- function(df,variables_to_keep,response,threshold){
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
        successive_timestamps <= threshold & successive_timestamps != 0, #put instead of AAA the lag for your response
        1, 
        0
      ))%>%
    ungroup()  
  
  y_trasl <- lag(response, default = 0) * df$ar_flag
  
  tot_obs_patient <- as.vector(table(df$CAI)[match(patients_mat, names(table(df$CAI)))])
  subj <- rep(seq_along(tot_obs_patient), times = tot_obs_patient)
  
  data_model = list(N = N, P = dim(covariates)[2], T = T, 
                    subj = subj, y = response, y_trasl = y_trasl, X = covariates)
  
  return(data_model)
}

bayes_R2 <- function(posterior_samples) {
  y_pred <- posterior_samples$y_obs_hat
  var_fit <- apply(y_pred, 1, var)
  var_res <- posterior_samples$sigma_e
  var_fit / (var_fit + var_res)
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

#Colesterolo_Hdl, Circonferenza_vita, Glucosio, PMAX, Trigliceridi, trained one each

target <- as.vector(df_stan$...) #add your target instead of ...
data_for_model <- data_stan(df_stan,chosen_columns,target,...) #ADD HERE YOUR THRESHOLD

#change values for simulation
fit = stan(file = 'model_autoregressive.stan', 
           data = data_for_model, 
           chains = 2, 
           iter = 1000, 
           warmup = 500, 
           cores = 2,
           thin = 1,
           seed = 19)

q <- traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
ggsave(paste0(folder_name, "/Traceplots of ....png"), plot = q, width = 6, height = 4) #add your target instead of ... but leave last . before png

#extract results and obtain residuals
posterior_samples <- extract(fit)

r2 <- bayes_R2(posterior_samples)
k <- ggplot(data.frame(r2 = r2), aes(x = r2)) +
  geom_histogram(binwidth = 0.0001, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = paste("Histogram of..."), #add your target instead of ...
    x = '...', #add your target instead of ...
    y = "Density"
  ) +
  theme_minimal()
ggsave(paste0(folder_name, "/Bayesian R2 of ....png"), plot = k, width = 6, height = 4) #add your target instead of ... but leave last . before png
print(mean(r2)) #save this value in the txt file created at the end of this code

#analyze performance
y_hat <- posterior_samples$y_hat
y_hat_point <- colMeans(y_hat, na.rm = TRUE)
res <- t(target - y_hat_point)

#save this values in the txt file created at the end of this code
rss <- sum(res^2)
cat('\nRSS: ',rss)
rmse <- sqrt(rss/dim(df_stan)[1]) 
cat('\nRMSE: ',rmse)
RMSE_to_sd_ratio <- rmse / sd(t(target))
cat('\nRMSE/SD: ',RMSE_to_sd_ratio)

#plot posterior distribution of sigma
p <- stan_dens(fit, pars = "sigma_e")
ggsave(paste0(folder_name, "/Sigma posterior distribution of ....png"), plot = p, width = 6, height = 4) #add your target instead of ... but leave last . before png

summary_stats <- summary(fit, pars = c("beta", "sigma_e"))
summary_text <- capture.output(print(summary_stats))
writeLines(summary_text, con = paste0(folder_name, "/Summary of ....txt")) #add your target instead of ... but leave last . before txt