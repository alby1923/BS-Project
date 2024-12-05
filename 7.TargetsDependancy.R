library(progress)
library(ggplot2)
library(rstan)
library(dplyr)

df_stan <- readRDS('ultimate_dataset.rds')

dir.create('STAN_...') #add your target instead of ...
folder_name <- paste0("STAN_...") #add your target instead of ...

#STAN MODEL -------
data_stan <- function(df,variables_to_keep,response){
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
  
  data_model = list(N = N, P = dim(covariates)[2], T = T, 
                    subj = subj, y = response, X = covariates)
  
  return(data_model)
}
bayes_R2 <- function(posterior_samples) {
  y_pred <- posterior_samples$y_hat
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

#Circonferenza_vita, Glucosio, PMAX trained one each

target <- as.vector(df_stan$...) #add your target instead of ...
data_for_model <- data_stan(df_stan,chosen_columns,target)

#change values for simulation
fit = stan(file = 'model_base.stan', 
           data = data_for_model, 
           chains = 2, 
           iter = 1000, 
           warmup = 500, 
           cores = 2,
           thin = 1,
           #control = list(adapt_delta = 0.95), may increase running time but better mix of chains. minimum value = 0.8 (default)
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

file_path <- file.path(folder_name, "Residuals_....Rdata") #add your target instead of ...but leave last . before Rdata
save(res, file = file_path)

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

#DEPENDANCY GRAPH --------------
library(BDgraph)

target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
res_list <- list()

for (target in target_variables){
  load(paste0("Residuals_", target, ".Rdata"))
  res_list[[target]] <- res
}

res_data <- do.call(rbind, res_list)

res_data <- t(res_data)

bdgraph = bdgraph.mpl(res_data,iter = 20000)
summary(bdgraph)

# Load necessary libraries
library(igraph)

# Define the p_links matrix (probability of connections)
p_links <- matrix(c(
  0, 0.95, 0, 1, 1,
  0, 0.00, 0, 0, 0,
  0, 0.00, 0, 1, 0, 
  0, 0.00, 0, 0, 1, 
  0, 0.00, 0, 0, 0  
), nrow = 5, byrow = TRUE)


# Create an igraph object
graph <- graph_from_adjacency_matrix(p_links, mode = "directed", diag = FALSE)

# Set vertex labels to the target variables
V(graph)$label <- target_variables

# Plot the graph
plot(graph, 
     vertex.size = 30,           # Adjust the vertex size
     vertex.color = "lightblue", # Color of the vertices
     vertex.label.color = "black",  # Color of the labels
     vertex.label.cex = 1.2,     # Font size of the labels
     edge.arrow.size = 0.5,      # Arrow size
     edge.width = 2,             # Edge width
     main = "Dependency Graph of Target Variables") # Title of the plot

#REAL DATA GRAPH ------------------
real_data <- df_stan[,target_variables]

bdgraph = bdgraph.mpl(real_data,iter = 20000)
summary(bdgraph)

# Load necessary libraries
library(igraph)

# Define the p_links matrix (probability of connections)
p_links <- matrix(c(0, 0, 0, 0, 0,
                    1, 0, 0, 0, 0,
                    1, 1, 0, 0, 0,
                    1, 1, 1, 0, 0,
                    1, 1, 1, 1, 0), 
                  nrow = 5, 
                  byrow = TRUE)

# Define the names of your target variables
target_variables <- c("Colesterolo_Hdl", "Circonferenza_vita", "Glucosio", "PMAX", "Trigliceridi")

# Create an igraph object
graph <- graph_from_adjacency_matrix(p_links, mode = "directed", diag = FALSE)

# Set vertex labels to the target variables
V(graph)$label <- target_variables

# Plot the graph
plot(graph, 
     vertex.size = 30,           # Adjust the vertex size
     vertex.color = "lightblue", # Color of the vertices
     vertex.label.color = "black",  # Color of the labels
     vertex.label.cex = 1.2,     # Font size of the labels
     edge.arrow.size = 0.5,      # Arrow size
     edge.width = 2,             # Edge width
     main = "Dependency Graph of Target Variables") # Title of the plot

