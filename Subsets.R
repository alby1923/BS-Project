library(progress)
library(ggplot2)

df_wide <- readRDS('flattened_dataset.rds')

#DATASET WITH ALL RESPONSES NOT NA --------
rows_to_remove = c()
names = c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
target_cols = which(colnames(df_wide) %in% names)

pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
for(i in 1:nrow(df_wide)){
  flag <- 0
  pb$tick()
  for(col in target_cols){
    if(flag == 1)
      break
    if(is.na(df_wide[i, col])){
      rows_to_remove <- c(rows_to_remove, i)
      flag <- 1
    }
  }
}
length(rows_to_remove)
df <- df_wide[-rows_to_remove,]

#PLOT TO SEE HOW MANY NA -----------
#create dataset with how many NA for each column
na_counts <- data.frame(
  Variable = c("Glucosio", "Circonferenza_vita", "Colesterolo_Hdl", "PMAX", "Trigliceridi"),
  NA_Count = c(
    sum(is.na(df_wide$Glucosio)),
    sum(is.na(df_wide$Circonferenza_vita)),
    sum(is.na(df_wide$Colesterolo_Hdl)),
    sum(is.na(df_wide$PMAX)),
    sum(is.na(df_wide$Trigliceridi))
  )
)
#barplot
ggplot(na_counts, aes(x = Variable, y = NA_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Variable", y = "NA count", title = "Count number of NA for each column") +
  ylim(0, 105000) +
  geom_hline(yintercept = nrow(df_wide), linetype = "dashed", color = "red") + 
  geom_hline(yintercept = nrow(df), linetype = "dashed", color = "red") + 
  theme_minimal() 

#DATASET WITHOUT ANY NA IN A ROW ----------
pb <- progress_bar$new(total = nrow(df)) #just to visualize how much time remaining
rows_to_remove <- c()
for(i in 1:nrow(df)){
  flag <- 0
  pb$tick()
  for(col in 1:ncol(df)){
    if(flag == 1)
      break
    if(is.na(df[i, col])){
      rows_to_remove <- c(rows_to_remove, i)
      flag <- 1
    }
  }
}
length(rows_to_remove)

df_vogliopiangere <- df[-rows_to_remove,]
