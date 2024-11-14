#LOAD DATA ------------
df <- readRDS('final_dataset.RDS')

#remove columns with all NULL values
df <- df[,-c(42,43)]
df$SESSO <- as.factor(df$SESSO)

#DATA INSPECTION -----------
#for each column (if numeric) we find the minimum and maximum values among all observation and save also the related id
not_numeric_columns <- c()
for (col in 1:ncol(df)){
  #check if the list is composed of char or numeric
  if(class(df[,col][[1]])=='numeric'){
    obs_min <- 0
    obs_max <- 0
    min <- Inf
    max <- -Inf
    my_col <- df[,col]
      for(i in 1:nrow(df)){
        if(min(na.omit(my_col[[i]])) < min){
           min <- min(my_col[[i]])
          obs_min <- i
           }
  if(max(na.omit(my_col[[i]])) > max){
    max <- max(my_col[[i]])
    obs_max <- i
  }
}
cat('Colonna ', colnames(df)[col],'\n')
cat('minimo ', min, ' del paziente ',obs_min, '. Massimo ',max,' del paziente ',obs_max,'\n')
print('---------------')
  }
  else not_numeric_columns <- c(not_numeric_columns, col)
}
#visual check if the not analyzed columns are dates
for(col in not_numeric_columns){
  cat(colnames(df)[col])
  cat('\n')
}

#we saw some zeros, hence we will check how many 
ferritina <- df$Ferritina
eosinofili <- df$Eosinofili_perc
for(i in 1:nrow(df)){
  vec <- na.omit(ferritina[[i]])
  flag <- 0
  for(j in 1:length(na.omit(ferritina[[i]])))
    if(vec[j]==0)
      flag=1
  if(flag==1)
    cat('obs ',i,'has a null value')
}

for(i in 1:nrow(df)){
  vec <- na.omit(eosinofili[[i]])
  flag <- 0
  for(j in 1:length(na.omit(eosinofili[[i]])))
    if(vec[j]==0)
      flag = flag+1
  if(flag>=1)
    cat('obs ',i,'has ',flag, 'null values\n')
}

#we saw a negative value, so we will check the observations with negative values in that column
alan <- df$Alanina_aminotransferasi_alt
for(i in 1:nrow(df)){
  vec <- na.omit(alan[[i]])
  flag <- 0
  for(j in 1:length(na.omit(alan[[i]])))
    if(vec[j]<0)
      flag=1
    if(flag==1)
      cat('obs ',i,'has a negative value')
}
#we found obs 496 has a negative value, we will just swap its sign
df$Alanina_aminotransferasi_alt[[496]][16] = 9

#check consistency of lengths of variables
test <- df[,c(2,3,8:75)]
flag=TRUE
j <- 1

while(flag){
  cat('variable ', colnames(test)[j+1])
  for(i in nrow(test)){
    if(length(test[i,j]) != length(test[i,j+1])){
      cat('observation ',i,'has length inconsistency\n')
    }
    print('-------------')
  }
  if(j+1 == ncol(test))
    flag = FALSE
  j <- j+2
}
#everything is fine

#OUTLIERS DETECTION ----------------
#to fix the heights the smartest choice is to replace in the vectors all the values with the mode, 
#because it is unlikely that someone significantly dropped or grew

##HEIGHT ----------------
#function mode doesn't esist in R
moda <- function(x) {
  unique_x <- unique(x)                 
  freq_table <- table(x)                 
  moda_value <- unique_x[which.max(freq_table)] 
  return(moda_value)
}

col <- 61
for(i in 1:nrow(df)){
  modav <- moda(df[[col]][[i]])
  for(j in 1:length(df[[col]][[i]])){
    df[[col]][[i]][[j]] <- modav
  }
} 

##PESO+CIRCONFERENZA VITA -----
rileva_e_correggi_spike <- function(serie,soglia){
  flag = 0
  if(length(serie)<3){
    med <- median(serie)
    diff <- abs(serie - med)
    if(any(diff > soglia)){
      flag <- 1
      serie[which(diff>soglia)] <- med
    }
  }
  else for(i in 2:(length(serie)-1)){
    aux <- c(serie[i-1],serie[i],serie[i+1])
    med <- median(aux)
    diff <- abs(aux - med)
    if(any(diff > soglia))
      flag <- 1
      if(diff[1]>soglia)
        serie[i-1] <- med
    if(diff[2]>soglia)
      serie[i] <- med
    if(diff[3]>soglia)
      serie[i+1] <- med
  }
  return(list(serie = serie, spike = flag))
}

#if in a sliding window of three it occurs a peak or a drop of more than 20 (or 30) it is swapped withe the median of that window
for(i in 1:nrow(df)){
to_print <- df$Peso[[i]]
risultato <- rileva_e_correggi_spike(df$Peso[[i]],20)
#print observation in which a change occurred
if(risultato$spike == 1){
  df$Peso[[i]] <- risultato$serie
cat('osservazione ',i,'\n')
cat("Serie originale:", to_print, "\n")
cat("Serie corretta:", risultato$serie, "\n\n\n")
}
}

for(i in 1:nrow(df)){
  to_print <- df$Circonferenza_vita[[i]]
  risultato <- rileva_e_correggi_spike(df$Circonferenza_vita[[i]],30)
  #print observation in which a change occurred
  if(risultato$spike == 1){
    df$Circonferenza_vita[[i]] <- risultato$serie
    cat('osservazione ',i,'\n')
    cat("Serie originale:", to_print, "\n")
    cat("Serie corretta:", risultato$serie, "\n\n\n")
  }
}

##PMAX+POLSO ----------
#for PMAX and POLSO we correct values which are unrealistic with the median o
out_of_bounds <- function(serie,lower_bound,upper_bound){
  flag <- 0
  if(any(serie < lower_bound) | any(serie > upper_bound)){
    serie[which(serie<lower_bound | serie>upper_bound)] <- round(median(serie))
    flag <- 1
  }
  return(list(serie = serie, flag = flag))
}

wrist <- df$Polso
pmax <- df$PMAX

#these bounds are as much conservative as possible. no-sense to accept values out of this considering blood donors are healthy
for(i in 1:nrow(df)){
  serie_polso <- out_of_bounds(wrist[[i]],40,120)
  if(serie_polso$flag == 1){
    cat('POLSO\n')
    cat('osservazione ',i,'\n')
    cat("Serie originale:", wrist[[i]], "\n")
    cat("Serie corretta:", serie_polso$serie, "\n\n\n")
    df$Polso[[i]] <- serie_polso$serie
  }
  serie_pmax <- out_of_bounds(pmax[[i]],40,180)
  if(serie_pmax$flag == 1){
    cat('PMAX\n')
    cat('osservazione ',i,'\n')
    cat("Serie originale:", pmax[[i]], "\n")
    cat("Serie corretta:", serie_pmax$serie, "\n\n\n")
    df$PMAX[[i]] <- serie_pmax$serie
  }
}

#LOG-TRANSFORMATION OF RESPONSES -------------
names <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
target_cols <- which(colnames(df) %in% names)

for(index in target_cols){
  for(i in 1:nrow(df)){
    df[[index]][[i]] <- log(df[[index]][[i]]) 
  }
}

#CREATE FINAL FILE --------------
saveRDS(df,'cleaned_dataset.rds')
