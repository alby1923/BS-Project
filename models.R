## Boring creation of subsets (RUN ONCE TO GET THE FILES) ---------------
library(progress)
library(ggplot2)

df_wide <- readRDS('flattened_dataset.rds')

##glucosio
pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
rows_to_remove = c()
for(i in 1:nrow(df_wide)){
  pb$tick()
  if(is.na(df_wide[i, 18])){
      rows_to_remove <- c(rows_to_remove, i)
  }
}
length(rows_to_remove)
df_glucosio = df_wide[-rows_to_remove,]
saveRDS(df_glucosio,'df_glucosio.rds')

##circonferenza vita
pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
rows_to_remove = c()
for(i in 1:nrow(df_wide)){
  pb$tick()
  if(is.na(df_wide[i, 36])){
    rows_to_remove <- c(rows_to_remove, i)
  }
}
length(rows_to_remove) ## !!!
df_circonferenza_vita = df_wide[-rows_to_remove,]
saveRDS(df_circonferenza_vita,'df_circonferenza_vita.rds')

##trigliceridi
pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
rows_to_remove = c()
for(i in 1:nrow(df_wide)){
  pb$tick()
  if(is.na(df_wide[i, 32])){
    rows_to_remove <- c(rows_to_remove, i)
  }
}
length(rows_to_remove)
df_trigliceridi = df_wide[-rows_to_remove,]
saveRDS(df_trigliceridi,'df_trigliceridi.rds')
##PMAX

pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
rows_to_remove = c()
for(i in 1:nrow(df_wide)){
  pb$tick()
  if(is.na(df_wide[i, 22])){
    rows_to_remove <- c(rows_to_remove, i)
  }
}
length(rows_to_remove)
df_PMAX = df_wide[-rows_to_remove,]
saveRDS(df_PMAX,'df_PMAX.rds')

##colesterolo
pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
rows_to_remove = c()
for(i in 1:nrow(df_wide)){
  pb$tick()
  if(is.na(df_wide[i, 6])){
    rows_to_remove <- c(rows_to_remove, i)
  }
}
length(rows_to_remove)
df_colesterolo = df_wide[-rows_to_remove,]
saveRDS(df_colesterolo,'df_colesterolo.rds')

## models -------------
library(rstan)

#glucosio
df_glucosio = readRDS("df_glucosio.rds")

patients_glucosio = as.matrix(unique(df_glucosio[,1]))

N_glucosio = dim(patients_glucosio)[1]
T_glucosio = dim(df_glucosio)[1]

data = df_glucosio$Glucosio
covariates = df_glucosio[,-c(1,2,18,36,32,22,6)]

data_glucosio = list(N = N_glucosio, P = 30, T = T_glucosio, 
                     subj = patients_glucosio, 
                     y = data, X = covariates)

fit_glucosio = stan(file = 'model_base.stan', data = data_glucosio)
