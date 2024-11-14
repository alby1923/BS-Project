# CANCELLARE DATI CON NA NELLE 5 TARGET

rows_to_remove = c()
names = c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
target_cols = which(colnames(df_wide) %in% names)
for(i in 1:nrow(df_wide)){
  flag <- 0
  for(col in target_cols){
    if(flag==1)
      break
    if(is.na(df_wide[i,col])){
      rows_to_remove <- c(rows_to_remove,i)
      flag <- 1
    }
  }
}

rows_to_remove <- unique(rows_to_remove)
length(rows_to_remove)
df_no_NA_target = df_wide[-rows_to_remove,]

df_dead <- df_wide[-rows_to_remove,]

length(unique(df_dead$CAI))

table(df_dead$CAI)

table(df_wide$CAI)
