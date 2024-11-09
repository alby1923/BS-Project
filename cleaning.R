#load data
df <- readRDS('final_dataset.RDS')

#remove columns with all NULL values
df <- df[,-c(42,43)]
df$SESSO <- as.factor(df$SESSO)

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

#we compute variances for numeric columns so to spot the vectors with anomalous values
#first we save the indexes of numeric columns
numeric_columns <- setdiff(1:81, not_numeric_columns)
known_columns <- c(13,15,61,63,65)
indexes_to_remove <- c()

for(col in 1:ncol(df)){
  if (col %in% not_numeric_columns | col %in% known_columns | col==1) { #we treat known columns in a different way, and we don't consider CAI variable
    next  #next iteration if column not numeric
  }
  cat('\n--------------------\n')
  cat('colonna ',colnames(df)[col])
  vars <- c() #create an auxiliary vector
  for(j in 1:nrow(df)){
    if(length(df[[col]][[j]])>1) #need at least 2 value to compute variance
    vars <- c(vars,var(df[[col]][[j]])) #compute variance of the i-th vector in j-th column
    else vars <- c(vars,0) #if the length is 1, then we have 1 value and so variance is zero
  }
  #once you have a vector with all variances, we compute z-scores to see if there are outliers (so vectors with really high variances)
  # Calculate mean and standard deviation of variances
  mean_var <- mean(vars)
  sd_var <- sd(vars)
  # Calculate Z-scores
  z_scores <- (vars - mean_var) / sd_var
  # Find observations with high Z-scores
  outliers <- which(abs(z_scores) > 15) #how much to set?? 15 means that the value is 15 standard deviations far from the mean (15 is HUGE, usually as a rule of
  #thumb a value higher than 3 si considered an outlier)
  if(length(outliers)>=1) #if there are no outliers nothing to print
  for(k in 1:length(outliers))
    cat('\n values of the outlier ',outliers[k],'with z-score ',z_scores[outliers[k]],': ',df[[col]][[outliers[k]]]) #print the vector of values of the outlier
  indexes_to_remove <- c(indexes_to_remove,outliers) #observation to remove
}
indexes_to_remove <- unique(indexes_to_remove) #removes eventual duplicates
length(indexes_to_remove) #see how many observations you would remove
#which value to choose for z_score?? are we really sure to remove anything since we don't have a medical background?? let's talk about it...

#to fix the heights the smartest choice is to replace in the vectors all the values with the mode, because it is unlikely that someone significantly dropped or grew
known_columns <- c(13,63,65) #remove the index of the height and circ. vita
col <- 61

#function mode doesn't esist in R
moda <- function(x) {
  unique_x <- unique(x)                 
  freq_table <- table(x)                 
  moda_value <- unique_x[which.max(freq_table)] 
  return(moda_value)
}

for(i in 1:nrow(df)){
  modav <- moda(df[[col]][[i]])
  for(j in 1:length(df[[col]][[i]])){
    df[[col]][[i]][[j]] <- modav
  }
} 

#from visual inspection, circ. vita needs to be fixed in just one obs
df$Circonferenza_vita[[744]][3] <- 106 #found by previous inspection (run following code)

#saves <- df
#df <- saves

for(iter in 1:5){#by visual inspection, 5 seems to be enough
  cat('\n\n\n ITERATION NUMBER ',iter,'\n\n')
  if(iter==3)
    known_columns <- c(63,65) #PMAX fixes after 2 iterations
for(col in known_columns){
  cat('\n--------------------\n')
  cat('colonna ',colnames(df)[col])
  vars <- c() #create an auxiliary vector
  for(j in 1:nrow(df)){
    if(length(df[[col]][[j]])>1) #need at least 2 value to compute variance
      vars <- c(vars,var(df[[col]][[j]])) #compute variance of the i-th vector in j-th column
    else vars <- c(vars,0) #if the length is 1, then we have 1 value and so variance is zero
  }
  #once you have a vector with all variances, we compute z-scores to see if there are outliers (so vectors with really high variances)
  # Calculate mean and standard deviation of variances
  mean_var <- mean(vars)
  sd_var <- sd(vars)
  # Calculate Z-scores
  z_scores <- (vars - mean_var) / sd_var
  # Find observations with high Z-scores
  outliers <- which(abs(z_scores) > 5)
  if(length(outliers)>=1) #if there are no outliers nothing to print
    for(k in 1:length(outliers))
      cat('\n values of the outlier ',outliers[k],'with z-score ',z_scores[outliers[k]],': ',df[[col]][[outliers[k]]]) #print the vector of values of the outlier
  for(l in outliers){ #for each outlier we take the index of the maximum or minimum (depends which is the farthest one from the median) and we substitute
    #that value with the median (rounded) of that vector
    if(max(df[[col]][[l]] - median(df[[col]][[l]])) > abs(min(df[[col]][[l]] - median(df[[col]][[l]]))))
    idx <- which.max(df[[col]][[l]])
    else idx <- which.min(df[[col]][[l]])
    df[[col]][[l]][[idx]] <- round(median(df[[col]][[l]]))
  }
}
}
#here the problem is that you need to run this for loop multiple times to detect all the outliers, but it is not the same amount of times for all the columns,
#hence we need to do experiments (visual inspect if everything seems fine) to see how many times is needed for each column

#count the NA and for each cell it says how many NA there are in the vector
conta_na <- function(df) {
  #matrix initialization
  na_count_matrix <- matrix(0, nrow = nrow(df), ncol = ncol(df),
                            dimnames = list(rownames(df), colnames(df)))
  
  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      # Conta i NA nel vettore nella cella (i, j)
      na_count_matrix[i, j] <- sum(is.na(df[[j]][[i]]))
    }
  }
  
  return(na_count_matrix)
}

na_matrix <- conta_na(df)
sum(na_matrix > 0) #how many cells have a value larger than zero


#log transform numerical variables





