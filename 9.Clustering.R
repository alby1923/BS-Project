library(rstan)
library(dplyr)
library(ggplot2)

df_original <- readRDS('datasets//ultimate_dataset.rds')

# Select a subset of the data if desired
df_sub <- df_original %>% slice(1:20000)  # first 5000 observations

# Check column names
colnames(df_sub)
# Preview the relevant columns
head(df_sub %>% select(SESSO, AB0, Rh))

# Summarize the data
summary(df_sub %>% select(SESSO, AB0, Rh))

# Identify target variables
target_variables <- c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")

# we have to revert to factors(Sex and Blood Type) to their intial categorical form

# ==========================
# Check for scaling attributes
attributes(df_sub$SESSO)
attributes(df_sub$AB0)
attributes(df_sub$Rh)

# Reverse scaling
df_sub <- df_sub %>%
  mutate(
    SESSO_original = SESSO * attr(SESSO, "scaled:scale") + attr(SESSO, "scaled:center"),
    AB0_original = AB0 * attr(AB0, "scaled:scale") + attr(AB0, "scaled:center"),
    Rh_original = Rh * attr(Rh, "scaled:scale") + attr(Rh, "scaled:center")
  )

# Round to integers
df_sub <- df_sub %>%
  mutate(
    SESSO_recoded = round(SESSO_original),
    AB0_recoded = round(AB0_original),
    Rh_recoded = round(Rh_original)
  )

# Frequency tables
table(df_sub$SESSO_recoded)
table(df_sub$AB0_recoded)
table(df_sub$Rh_recoded)

# Create the group_var
df_sub <- df_sub %>%
  mutate(group_var = interaction(SESSO_recoded, AB0_recoded, Rh_recoded, drop = TRUE)) %>%
  mutate(group_var = as.integer(group_var))  # Convert to integers

# Check the number of unique groups
length(unique(df_sub$group_var))  # Should be <= 16

# Frequency table of group_var
table(df_sub$group_var)

# Contingency table for all combinations
table(df_sub$SESSO_recoded, df_sub$AB0_recoded, df_sub$Rh_recoded)

# ==========================

df_sub <- df_sub %>%
  select(
    group_var, 
    all_of(target_variables),  # Retain target variables
    CAI, Date
  )

# Validate the cleaned dataset
print(colnames(df_sub))  # Check remaining columns
print(dim(df_sub))       # Check dimensions

Y <- df_sub %>% select(all_of(target_variables)) %>% as.matrix()  # Target variables matrix

stan_data <- list(
  N = nrow(Y),                      # Number of observations
  D = ncol(Y),                      # Number of target variables
  G = length(unique(df_sub$group_var)),  # Number of unique groups
  K = 20,                           # Truncation level for DP
  Y = Y,                            # Target variables matrix
  group = df_sub$group_var,         # Group variable
  alpha = 1.0                       # DP concentration parameter
)

stan_model_path <- "./BS-Project/model_cluster.stan"  # Path to the Stan model file

fit <- stan(
  file = stan_model_path,  # Stan model
  data = stan_data,        # Prepared data
  iter = 10,             # Total iterations
  warmup = 1,           # Warmup iterations
  chains = 4,              # Number of chains
  cores = 4,               # Number of cores
  seed = 123               # Seed for reproducibility
)

# ==========================
# Analyze the results

print(fit, pars = c("v", "mu_raw", "gamma_raw", "sigma"))

# Extract posterior samples
results <- extract(fit)

calc_pi <- function(v_samps) {
  M <- nrow(v_samps)
  K <- ncol(v_samps)
  pi_out <- matrix(0, M, K)
  for (m in 1:M) {
    stick_elems <- numeric(K)
    stick_elems[1] <- v_samps[m, 1]
    for (kk in 2:K) {
      stick_elems[kk] <- v_samps[m, kk] * prod(1 - v_samps[m, 1:(kk-1)])
    }
    pi_out[m, ] <- stick_elems
  }
  pi_out
}

pi_samples <- calc_pi(results$v)
pi_mean <- colMeans(pi_samples)

# Barplot of cluster weights
barplot(pi_mean, main = "Average Cluster Weights", xlab = "Cluster", ylab = "Weight")

gamma_samples <- results$gamma_raw  # Dimensions: D x G x iterations
gamma_means <- apply(gamma_samples[1, , ], 2, mean)  # Mean for the first variable

barplot(gamma_means, main = "Group-Level Effects (Variable 1)", xlab = "Group", ylab = "Effect")

# Plot trace for the first few stick-breaking weights (v[1:5])
traceplot(fit, pars = c("v[1]", "v[2]", "v[3]", "v[4]", "v[5]"))
# Trace plots for group-level effects (gamma_raw) for the first group
traceplot(fit, pars = c("gamma_raw[1,1]", "gamma_raw[1,2]", "gamma_raw[1,3]"))
# Trace plots for global parameters like sigma
traceplot(fit, pars = "sigma")



