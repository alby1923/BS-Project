library(dplyr)
library(ggplot2)
library(pROC)
library(progress)
library(patchwork)

#CREATE DATASET WITH PERCENTAGES ------------
#load data
df_stan_last <- readRDS('risk_scores_real.rds')
#build dataset with percentages of being in ranges for metabolic syndrome in that target variable
#Initialize dataset
df_classifier <- data.frame(
  Circonferenza_vita_percentage = numeric(nrow(df_stan_last)),
  Trigliceridi_percentage = numeric(nrow(df_stan_last)),
  Colesterolo_Hdl_percentage = numeric(nrow(df_stan_last)),
  Glucosio_percentage = numeric(nrow(df_stan_last)),
  PMAX_percentage = numeric(nrow(df_stan_last)),
  Metabolic_syndrome_probability = numeric(nrow(df_stan_last)),
  Risk_score = df_stan_last$risk_score
)

thresholds_m <- log(c(102,150,40,100,130))
thresholds_f <- log(c(88,150,50,100,130))

target_variables <- c('Circonferenza_vita','Trigliceridi','Colesterolo_Hdl','Glucosio','PMAX')
j=0
#compute probability y_test_point > threshold
for(target in target_variables){
  #iterate over columns of probabilities
  j <- j+1
  #load fit
  load(paste0('final_model_',target,'.Rdata'))
  #extract results
  output <- extract(fit)
  y_test <- output$y_test
  y_test_point <- colMeans(y_test, na.rm = TRUE)
  sigma <- mean(output$sigma_e)
  #for colestherol we need to check if value is below threshold, for the others if it is above
  if(target != 'Colesterolo_Hdl'){
    for (i in 1:ncol(y_test)) {
      #since distribution is a gaussian, we are provided of the point estimate of the mean and of the error
        if (df_stan_last$SESSO[i] == 'M') {
          df_classifier[i,j] <- pnorm(thresholds_m[j],mean = y_test_point[i],sd = sigma,lower.tail = FALSE)
        } 
        else {
          df_classifier[i,j] <- pnorm(thresholds_f[j],mean = y_test_point[i],sd = sigma,lower.tail = FALSE)
        }
    }
  }
  else {
    for (i in 1:ncol(y_test)) {
      if (df_stan_last$SESSO[i] == 'M') {
        df_classifier[i,j] <- pnorm(thresholds_m[j],mean = y_test_point[i],sd = sigma,lower.tail = TRUE)
      } 
      else {
        df_classifier[i,j] <- pnorm(thresholds_f[j],mean = y_test_point[i],sd = sigma,lower.tail = TRUE)
      }
    }
  }
}

#save dataset
saveRDS(df_classifier,'probabilities.rds')
#COMPUTE PROBABILITY OF METABOLIC SYNDROME -------------
df_classifier <- readRDS('probabilities.rds')
percentages <- df_classifier[,1:5]

calculate_combinations_with_sum <- function(probabilities) {
  
  #generate all binary combinations, set 1 if you want to take probability that observation i has that target in ranges of metabolic syndrome, 0 otherwise
  binary_combinations <- expand.grid(rep(list(c(0, 1)), 5))
  
  #compute total probability of that combination, so if set to one will take its probability while if set to zero will take 1-probability
  #the ^combination because if value set to zero than probabilities computes as well the product but elevates it to 0 so to obtain 1 and not
  #involve it in the final computation
  results <- apply(binary_combinations, 1, function(combination) {
    prod(probabilities ^ combination * (1 - probabilities) ^ (1 - combination))
  })
  
  #add total probability to each combination and counts how many 1 (so how many target in ranges) occurred in combination
  binary_combinations$num_ones <- rowSums(binary_combinations)
  binary_combinations$probability <- results
  
  #sum probabilities with at least 3 ones, hence to have metabolic syndrome
  sum_probabilities_at_least_3 <- sum(binary_combinations$probability[binary_combinations$num_ones >= 3])
  
  return(sum_probabilities_at_least_3)
}

#progress bar to see time remaining
pb <- txtProgressBar(min = 0, max = nrow(df_classifier), style = 3)
for(i in 1:nrow(df_classifier)){
  df_classifier$Metabolic_syndrome_probability[i] <- calculate_combinations_with_sum(percentages[i,])
  #update progress bar
  setTxtProgressBar(pb, i)
}
close(pb)

saveRDS(df_classifier,'completed_classifier.rds')
#ROC CURVE -------------
df_classifier <- readRDS('completed_classifier.rds')
#definition of real class (1 = metabolic syndrome, 0 = no metabolic syndrome)
df_classifier <- df_classifier %>%
  mutate(binary_target = ifelse(Risk_score > 2, 1, 0))

#compute ROC curve
roc_curve <- roc(
  response = df_classifier$binary_target,
  predictor = df_classifier$Metabolic_syndrome_probability,
  levels = c(0, 1)  #(1 = metabolic syndrome, 0 = no metabolic syndrome)
)

#extract ROC curve metrics
auc_value <- auc(roc_curve)  #compute AUC
roc_data <- data.frame(
  TPR = rev(roc_curve$sensitivities),  # True Positive Rate
  FPR = rev(1 - roc_curve$specificities)  # False Positive Rate
)

#plot
ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.2) +  # Linea ROC
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Linea casuale
  labs(
    title = paste("ROC Curve ( AUC =", round(auc_value, 3), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


#CLASSIFIER PERFORMANCE -------------
#compute percentages of classes in the test set
risk_score_summary <- df_classifier %>%
  summarise(
    perc_risk_above_2 = mean(Risk_score > 2) * 100,
    perc_risk_below_eq_2 = mean(Risk_score <= 2) * 100
  )
print(risk_score_summary)

#create column and YES if patient has metabolic syndrome, NO otherwise and set position just for the plot
df_classifier <- df_classifier %>%
  mutate(
    color_group = ifelse(Risk_score <= 2, "NO", "YES")
  )

#find the threshold which sets to zero the misclassifications of people having metabolic syndrome classified as not
df_classifier_sorted <- df_classifier %>%
  arrange(Metabolic_syndrome_probability)

#so we need to find the threshold which classifies correctly all the yes
optimal_threshold <- min(df_classifier_sorted$Metabolic_syndrome_probability[df_classifier_sorted$color_group == "YES"])
print(optimal_threshold)

#max 5% of YES misclassified
#optimal_threshold <- quantile(df_classifier_sorted$Metabolic_syndrome_probability[df_classifier_sorted$color_group == "YES"], 0.05)
#print(optimal_threshold)

#creates confusion matrix
df_classifier <- df_classifier %>%
  mutate(
    predicted_class = ifelse(Metabolic_syndrome_probability >= optimal_threshold, "YES", "NO")
  )
confusion_matrix <- table(
  Actual = df_classifier$color_group,
  Predicted = df_classifier$predicted_class
)
print(confusion_matrix)

#metrics of performance
true_positive <- confusion_matrix["YES", "YES"]
true_negative <- confusion_matrix["NO", "NO"]
false_positive <- confusion_matrix["NO", "YES"]
false_negative <- confusion_matrix["YES", "NO"]

accuracy <- (true_positive + true_negative) / sum(confusion_matrix)
sensitivity <- true_positive / (true_positive + false_negative)  # Recall
specificity <- true_negative / (true_negative + false_positive)
precision <- true_positive / (true_positive + false_positive)
f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)

metrics <- data.frame(
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision = precision,
  F1_Score = f1_score
)
print(metrics)

#PLOT OF DISTRIBUTIONS --------------
##FIRST VERSION ----------
ggplot(df_classifier, aes(x = Metabolic_syndrome_probability, fill = color_group)) +
  geom_histogram(binwidth = 0.05, alpha = 0.7, position = "identity", color = "black") +
  facet_wrap(~ color_group, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("NO" = "green", "YES" = "red")) +
  labs(
    title = "Distribution of Metabolic Syndrome Probabilities",
    x = "Metabolic Syndrome Probability",
    y = "Count",
    fill = "Metabolic Syndrome"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold")
  )
##SECOND VERSION -------------
plot_yes <- ggplot(df_classifier %>% filter(color_group == "YES"), 
                   aes(x = Metabolic_syndrome_probability)) +
  geom_histogram(bins = 50, fill = "red", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of Probabilities for YES (Metabolic Syndrome)",
    x = "Metabolic Syndrome Probability",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
plot_no <- ggplot(df_classifier %>% filter(color_group == "NO"), 
                  aes(x = Metabolic_syndrome_probability)) +
  geom_histogram(bins = 50, fill = "green", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of Probabilities for NO",
    x = "Metabolic Syndrome Probability",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

#combine two plots
plot_yes / plot_no

