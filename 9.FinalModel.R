library(rstan)

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
          successive_timestamps <= lags[i] & successive_timestamps > 0, #if zero it is the first observation of the donor
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
                    y_trasl = y_trasl,
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

lags <- c(100,200,300,400) #max window allowed for previous observation in autoregressive component
data_model <- data_stan(df_reduced,chosen_columns,target_variables,lags)

fit = stan(file = 'model_final.stan', 
           data = data_model, 
           chains = 1, 
           iter = 20, 
           warmup = 5, 
           cores = 1,
           thin = 1,
           seed = 19)

