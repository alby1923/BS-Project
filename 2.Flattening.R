# FLATTEN THE DATA

# Libraries
library(dplyr)
library(tidyr)
library(purrr)

# Split of the dataset in time-variant variables and time invariant
df_Arrigoni <- readRDS("cleaned_dataset.rds")
df_timeinv <- df_Arrigoni[, c("CAI", "SESSO", "DATA_NASCITA", "Rh", "AB0")]
df_bloodvalues <- df_Arrigoni %>% select( -SESSO, -DATA_NASCITA, -Rh, -AB0,
                                          -Eta_glucosio, -Eta_colesterolo, 
                                          -Eta_trigliceridi, -Eta_circonferenza, 
                                          -Eta_PMAX, -Eta_min)
df_bloodvalues <- as.data.frame(df_bloodvalues)

# convert all in list to make code easier
for(i in 1:nrow(df_bloodvalues)){
  for(j in 2:ncol(df_bloodvalues)){
    if(!is.list(df_bloodvalues[i,j])){
      df_bloodvalues[i,j] <- list(df_bloodvalues[i,j])
    }
  }
}

# convert all factor in numeric
for (i in 1:nrow(df_bloodvalues)) {
  df_bloodvalues$Alcool[[i]] <- as.numeric(as.character(df_bloodvalues$Alcool[[i]]))
  df_bloodvalues$Fumo[[i]] <- as.numeric(as.character(df_bloodvalues$Fumo[[i]]))
  df_bloodvalues$Attivita_fisica[[i]] <- as.numeric(as.character(df_bloodvalues$Attivita_fisica[[i]]))
}

# Date e valori corrispondenti
colonne_date <- grep("^Data_", names(df_bloodvalues), value = TRUE)
#colonne_date
colonne_valori <- sub("^Data_", "", colonne_date)
#colonne_valori


# Pairs to flatten
variable_pairs <- lapply(seq_along(colonne_date), function(i) {
  list(date = colonne_date[i], observation = colonne_valori[i])
})
#variable_pairs

# Function to determine the type of elements within a list column, checking for specific types
get_list_type <- function(column) {
  if (all(sapply(column, function(x) is.double(unlist(x)))))
    return('double')
  else if (all(sapply(column, function(x) is.numeric(unlist(x))))) 
    return('numeric')
  else if (all(sapply(column, function(x) is.character(unlist(x)) || is.factor(unlist(x))))) 
    return('categorical')
  else {
    return('mixed')
  }
}

# Determine the type of each observation column
variable_types <- data.frame(
  Variable = colonne_valori,
  Type = sapply(colonne_valori, function(var) {
    if (var %in% names(df_bloodvalues)) {
      get_list_type(df_bloodvalues[[var]])  # Use the function to get the type of the list column
    } else {
      NA  # In case the column is not present, return NA
    }
  }),
  Aggregation = 'mean'  # Default aggregation method, adjust as needed
)

# Adjust Aggregation based on Type
variable_types$Aggregation <- ifelse(
  variable_types$Type == 'double' | variable_types$Type == 'numeric', 'mean', 
  ifelse(variable_types$Type == 'categorical', 'mode', 'N/A')
)
variable_types


# Define variable types and aggregation methods
variable_types <- data.frame(
  Variable = colonne_valori,
  Type = c(rep('numeric', 35)),
  Aggregation = c(rep('mean', 35))  # Adjust as needed
)
variable_types


# Initialize a list to store the long-format data frames
df_long_list <- list()

# Process each variable pair
for (pair in variable_pairs) {
  date_col <- pair$date
  obs_col <- pair$observation
  variable_name <- obs_col
  
  df_temp <- df_bloodvalues %>%
    select(CAI, all_of(date_col), all_of(obs_col)) %>%
    unnest(cols = c(all_of(date_col), all_of(obs_col))) %>%
    rename(Date = !!date_col, Value = !!obs_col) %>%
    mutate(Variable = variable_name)
  
  df_long_list[[variable_name]] <- df_temp
}

# Combine all long-format data frames
df_long <- bind_rows(df_long_list)

# Aggregate the data
df_long_agg <- df_long %>%
  group_by(CAI, Date, Variable) %>%
  summarize(
    Value = {
      var_type <- variable_types$Type[variable_types$Variable == first(Variable)]
      aggregation <- variable_types$Aggregation[variable_types$Variable == first(Variable)]
      values <- Value
      
      if (var_type == 'numeric') {
        values <- as.numeric(values)
        if (aggregation == 'mean') {
          mean(values, na.rm = TRUE)
        } else if (aggregation == 'sum') {
          sum(values, na.rm = TRUE)
        } else {
          mean(values, na.rm = TRUE)  # Default to mean
        }
      } else {
        # For categorical variables
        if (aggregation == 'mode') {
          # Return the mode
          ux <- unique(values)
          ux[which.max(tabulate(match(values, ux)))]
        } else if (aggregation == 'concatenate') {
          paste(unique(values), collapse = ", ")
        } else {
          paste(unique(values), collapse = ", ")
        }
      }
    },
    .groups = 'drop'
  )

# Pivot the data to wide format
df_wide <- df_long_agg %>%
  pivot_wider(names_from = Variable, values_from = Value)

# Convert Date column to Date type
df_wide$Date <- as.Date(df_wide$Date)

saveRDS(df_wide,'flattened_dataset.rds')
