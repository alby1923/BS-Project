library(rstan)
library(dplyr)
library(ggplot2)

df_stan <- readRDS('ultimate_dataset.rds')

#TEST SET --------
#we will run our model only on patients with at least 5 observations in the training set (hence, at least 6 in total)
df_stan <- df_stan %>%
  group_by(CAI) %>%         
  filter(n() >= 6) %>%    
  ungroup()   

#take out from each patient its last observation since we will try to predict that one with our model
df_stan_last <- df_stan %>%
  group_by(CAI) %>%          
  slice_tail(n = 1) %>%      
  ungroup()                 

#save this version of test set for the stan model, the other one is for classifier
df_stan_last_model <- df_stan_last

#we remove form the data_stan (train set) the observations in the test set
df_stan <- df_stan %>%
  anti_join(df_stan_last, by = colnames(df_stan_last_model))

#for the test set we will only need the target variables and the sex factor since different thresholds are applied to male and females
target_variables <- c("Circonferenza_vita","Colesterolo_Hdl","Trigliceridi","PMAX","Glucosio")
df_stan_last <- df_stan_last[,c(target_variables,'SESSO')]
df_stan_last <- exp(df_stan_last) #we get back to original values
table(df_stan_last$SESSO) #just to see with which value males and females are saved in the dataset
df_stan_last$SESSO <- ifelse(df_stan_last$SESSO < 1, "M", "F")
df_stan_last$SESSO <- as.factor(df_stan_last$SESSO)

#here set 1 if value for that target we have value in rsnge for metabolic syndrome, zero otherwise
df_stan_last <- df_stan_last %>%
  # Glucosio > 100
  mutate(glucosio_flag = if_else(Glucosio > 100, 1, 0)) %>%
  # Trigliceridi > 150
  mutate(trigliceridi_flag = if_else(Trigliceridi > 150, 1, 0)) %>%
  # Pressione massima (PMAX) > 130
  mutate(pmax_flag = if_else(PMAX > 130, 1, 0)) %>%
  # Circonferenza vita > 102 maschi, > 88 femmine
  mutate(circonferenza_flag = if_else(
    (SESSO == 'M' & Circonferenza_vita > 102) | 
      (SESSO == 'F' & Circonferenza_vita > 88), 
    1, 
    0
  )) %>%
  # Colesterolo < 40 maschi, < 50 femmine
  mutate(colesterolo_flag = if_else(
    (SESSO == 'M' & Colesterolo_Hdl < 40) | 
      (SESSO == 'F' & Colesterolo_Hdl < 50), 
    1, 
    0
  )) %>%
  
  # Calcolo del punteggio complessivo di rischio
  mutate(risk_score = rowSums(across(ends_with("_flag"))))

#set 1 if at least three targets in range, hence metabolic syndrome
num_high_risk <- df_stan_last %>%
  filter(risk_score >= 3) %>%
  nrow()

#wee see how many of them and percentage with respect to total number of patients
num_high_risk
num_high_risk/dim(df_stan_last)[1]

#save datasets for following models
saveRDS(df_stan,'df_train_set.rds')
saveRDS(df_stan_last_model,'df_test_set.rds')
saveRDS(df_stan_last,'risk_scores_real.rds')
#STAN MODEL ---------------
data_stan <- function(df,variables_to_keep,response,df_test){
  #data needed for the model
  
  #total number of patients
  patients_mat <- as.matrix(unique(df[,1]))
  patients <- as.vector(patients_mat)
  N <- dim(patients_mat)[1]
  
  #total number of observations
  T <- dim(df)[1]
  
  #covariates and number of covariates
  covariates <- df[,variables_to_keep] #remove response (and useless) covariates
  covariates_test <- df_test[,variables_to_keep]
  P <- dim(covariates)[2]
  
  #index of patient related to each line of covariates (for random intercept)
  tot_obs_patient <- as.vector(table(df$CAI)[match(patients_mat, names(table(df$CAI)))])
  subj <- rep(seq_along(tot_obs_patient), times = tot_obs_patient)
  subj_test <- 1:dim(df_test)[1]
  
  data_model = list(N = N, P = P, T = T, 
                    subj = subj, y = response, X = covariates, Xtest = covariates_test, subj_test = subj_test)
  
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
  #"Rh",
  #"SESSO",
  #"AB0",
  'eta_std'
)

#for waist circumference we have no autoregressive component
response <- df_stan$Circonferenza_vita
response_test <- df_stan_last_model$Circonferenza_vita
data_model <- data_stan(df_stan,chosen_columns,response,df_stan_last_model)

#fit the model
fit = stan(file = 'model_final_A.stan', 
           data = data_model, 
           chains = 2, 
           iter = 1000, 
           warmup = 500, 
           cores = 2,
           thin = 1,
           seed = 19)

#TRACEPLOTS AND RESULTS --------------
#create folder to save results
dir.create('STAN_Final_Circonferenza_vita')
folder_name <- 'STAN_Final_Circonferenza_vita'
file_path <- file.path(folder_name, "final_model_Circonferenza_vita.Rdata") 
save(fit,file = file_path)

#save traceplot of beta and sigma
q <- traceplot(fit, pars = c(paste0("beta[", 1:length(chosen_columns), "]"), "sigma_e"))
ggsave(paste0(folder_name, "/Traceplots of Circonferenza_vita.png"), plot = q, width = 6, height = 4)

#save plot posterior distribution of sigma
p <- stan_dens(fit, pars = "sigma_e")
ggsave(paste0(folder_name, "/Sigma posterior distribution of Circonferenza_vita.png"), plot = p, width = 6, height = 4)

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
    title = paste("Histogram of Circonferenza_vita"),
    x = 'Circonferenza_vita',
    y = "Density"
  ) +
  theme_minimal()
ggsave(paste0(folder_name, "/Bayesian R2 of Circonferenza_vita.png"), plot = k, width = 6, height = 4)

#analyze performance on train set
y_hat_point <- colMeans(y_train, na.rm = TRUE)
res <- t(response - y_hat_point)

cat('\nTRAIN SET')
rss <- sum(res^2)
cat('\nRSS: ',rss)
rmse <- sqrt(rss/dim(df_stan)[1]) 
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
file_path <- file.path(folder_name, "Train_set_output_Circonferenza_vita.Rdata")
save(y_hat_point, file = file_path)

#save y_test output for classifier
file_path <- file.path(folder_name, "Test_set_output_Circonferenza_vita.Rdata")
save(y_test_point, file = file_path)

#save summary of stan model
summary_stats <- summary(fit, pars = c("beta","sigma_e"))
summary_text <- capture.output(print(summary_stats))
writeLines(summary_text, con = paste0(folder_name, "/Summary of Circonferenza_vita.txt"))



