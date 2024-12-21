library(rstan)
library(dplyr)
library(ggplot2)

#load data
df_train <- readRDS('df_train_set.rds')
df_test <- readRDS('df_test_set.rds')
classifier_df <- readRDS('df_classifier.rds')
df_stan_last <- readRDS('risk_scores_real.rds')
#STAN MODEL ---------------
data_stan <- function(df,target,variables_to_keep,response,threshold,df_test,previous_target,previous_predictions_train,previous_predictions_test){
  
  #compute translated vector for autoregressive component and set it to zero if previous observation is older than threshold
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
        successive_timestamps <= threshold & successive_timestamps != 0,
        1, 
        0
      ))%>%
    ungroup()  
  
  y_trasl <- lag(response, default = 0) * df$ar_flag
  
  #number of patients and total observations
  patients_mat <- as.matrix(unique(df[,1]))
  patients <- as.vector(patients_mat)
  N <- dim(patients_mat)[1]
  T <- dim(df)[1]
  
  #take covariates and number of covariates for train and test set
  covariates <- df[,variables_to_keep]
  covariates_test <- df_test[,variables_to_keep]
  P <- dim(covariates)[2]
  
  #indexes for each patient (for random intercept)
  tot_obs_patient <- as.vector(table(df$CAI)[match(patients_mat, names(table(df$CAI)))])
  subj <- rep(seq_along(tot_obs_patient), times = tot_obs_patient)
  subj_test <- 1:dim(df_test)[1]

  y_old <- df[,previous_target]
  y_std <- scale(y_old)
  
  #to rescale previous predictions it needs to be with the same mean and sd used for both train and test set and same magnitude of observed values
  col_means <- colMeans(y_old)
  col_sds <- apply(y_old, 2, sd)
  previous_predictions_train_std <- (previous_predictions_train - col_means)/col_sds
  previous_predictions_test_std <- (previous_predictions_test - col_means)/col_sds
  
  #compute if there is autoregressive component for predicted value
  df_tot <- bind_rows(df, df_test) %>%
    arrange(CAI, delta_date)
  
  df_tot <- df_tot %>%
    group_by(CAI) %>%
    arrange(delta_date, .by_group = TRUE) %>%
    mutate(
      successive_timestamps = if_else(
        delta_date == 0,
        0,               
        delta_date - lag(delta_date, default = 0)
      ),
      ar_flag = if_else(
        successive_timestamps <= threshold & successive_timestamps != 0,
        1, 
        0
      ))%>%
    ungroup()  
  
  df_last <- df_tot %>%
    group_by(CAI) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  df_last_value <- df %>%
    group_by(CAI) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  y_trasl_test <- df_last_value[,target] * df_last$ar_flag
  
  data_model = list(N = N, 
                    P = P, 
                    T = T, 
                    subj = subj,
                    y = response,
                    y_std = y_std,
                    y_trasl = y_trasl,
                    lag = threshold,
                    X = covariates,
                    K = length(previous_target),
                    y_pred_train_std = previous_predictions_train_std,
                    y_pred_test_std = previous_predictions_test_std,
                    Xtest = covariates_test,
                    subj_test = subj_test,
                    y_trasl_test = unlist(y_trasl_test))
  
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

#ORDER: CIRCONFERENZA VITA - TRIGLICERIDI - COLESTEROLO - GLUCOSIO - PMAX, besides the first one iteratively substitute the new target where ### occurs

#take response values
response <- df_train$PMAX ###
response_test <- df_test$PMAX ###
#autoregression threshold
threshold <- 180 ######### trigliceridi: 365, colesterolo: 365, glucosio: 0, pmax: 180
#order of targets
previous_target <- c('Circonferenza_vita','Trigliceridi','Colesterolo_Hdl','Glucosio')
#previous predicted values by preceding models
load("Train_set_output_Glucosio.Rdata") ###
y_old_train <- matrix(y_hat_point, ncol = 1)
rm(y_hat_point)
#pp_train <- y_old_train #use this just for trigliceridi since it is the first one for model_B
load('pp_train.Rdata')
pp_train <- cbind(pp_train,y_old_train)
load("Test_set_output_Glucosio.Rdata") ###
y_old_test <- matrix(y_test_point, ncol = 1)
rm(y_test_point)
#pp_test <- y_old_test #use this just for trigliceridi since it is the first one for model_B
load('pp_test.Rdata')
pp_test <- cbind(pp_test,y_old_test)

save(pp_train,file = 'pp_train.Rdata')
save(pp_test,file = 'pp_test.Rdata')

target <- 'PMAX' #put here the target you are running
data_model <- data_stan(df_train,target,chosen_columns,response,threshold,df_test,previous_target,pp_train,pp_test)

#fit the model
fit = stan(file = 'model_final_B.stan', 
           data = data_model, 
           chains = 2, 
           iter = 1000, 
           warmup = 500, 
           cores = 2,
           thin = 1,
           seed = 19)

#TRACEPLOTS AND RESULTS --------------
#create folder to save results
dir.create('STAN_Final_PMAX') ###
folder_name <- 'STAN_Final_PMAX' ###
file_path <- file.path(folder_name, "final_model_PMAX.Rdata")  ###
save(fit,file = file_path)

#save traceplot of beta and sigma
q <- traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
ggsave(paste0(folder_name, "/Traceplots of PMAX.png"), plot = q, width = 6, height = 4) ###

#save plot posterior distribution of sigma
p <- stan_dens(fit, pars = "sigma_e")
ggsave(paste0(folder_name, "/Sigma posterior distribution of PMAX.png"), plot = p, width = 6, height = 4) ###

#save plot for autoregressive component
t <- traceplot(fit, pars = "phi")
ggsave(paste0(folder_name, "/Traceplot of autoregressive coefficient of PMAX.png"), plot = t, width = 6, height = 4) ###

#save plot for gamma
l <- traceplot(fit, pars = paste0("gamma[", 1:length(previous_target), "]"))
ggsave(paste0(folder_name, "/Traceplots of gamma of PMAX.png"), plot = l, width = 6, height = 4) ###

#extract results
output <- extract(fit)
y_train <- output$y_train
y_test <- output$y_test
sigma <- output$sigma_e

bayes_R2 <- function(y_pred,sigma_e) {
  var_fit <- apply(y_pred, 1, var)
  var_res <- sigma_e
  var_fit / (var_fit + var_res)
}

r2 <- bayes_R2(y_train,sigma)
k <- ggplot(data.frame(r2 = r2), aes(x = r2)) +
  geom_histogram(binwidth = 0.0001, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = paste("Histogram of PMAX"), ###
    x = 'PMAX', ###
    y = "Density"
  ) +
  theme_minimal()
ggsave(paste0(folder_name, "/Bayesian R2 of PMAX.png"), plot = k, width = 6, height = 4) ###

#analyze performance on train set
y_hat_point <- colMeans(y_train, na.rm = TRUE)
res <- t(response - y_hat_point)

cat('\nTRAIN SET')
rss <- sum(res^2)
cat('\nRSS: ',rss)
rmse <- sqrt(rss/dim(df_train)[1]) 
cat('\nRMSE: ',rmse)
RMSE_to_sd_ratio <- rmse / sd(t(response))
cat('\nRMSE/SD: ',RMSE_to_sd_ratio)

#analyze performance on test set
y_test_point <- colMeans(y_test, na.rm = TRUE)
res_test <- t(response_test - y_test_point)

cat('\nTEST SET')
rss_test <- sum(res_test^2)
cat('\nRSS: ',rss_test)
rmse_test <- sqrt(rss_test/dim(df_stan_last)[1]) 
cat('\nRMSE: ',rmse_test)
RMSE_to_sd_ratio_test <- rmse_test / sd(t(response_test))
cat('\nRMSE/SD: ',RMSE_to_sd_ratio_test)

#save y_train output for following models
file_path <- file.path(folder_name, "Train_set_output_PMAX.Rdata") ###
save(y_hat_point, file = file_path)

#save y_test output for classifier
file_path <- file.path(folder_name, "Test_set_output_PMAX.Rdata") ###
save(y_test_point, file = file_path)

#save summary of stan model
summary_stats <- summary(fit, pars = c("beta","sigma_e",'phi','gamma'))
summary_text <- capture.output(print(summary_stats))
writeLines(summary_text, con = paste0(folder_name, "/Summary of PMAX.txt")) ###

