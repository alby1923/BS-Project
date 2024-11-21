library(progress)
library(ggplot2)
library(dplyr)

df_wide <- readRDS('flattened_dataset.rds')

#FILLER OF NA IN DATASET ----------------

library(lubridate) # For handling date differences
library(data.table)  # Added for efficient data manipulation

# Ensure the Date column is in Date format
df_wide <- df_wide %>%
  mutate(Date = as.Date(Date))

# Function to fill missing values based on the nearest available value within a specified delay
fill_missing_values <- function(data, cols_to_fill, id_col, date_col, delay_days) {
  # Create a progress bar
  pb <- progress_bar$new(
    format = "  Filling [:bar] :percent eta: :eta",
    total = length(cols_to_fill), clear = FALSE, width=60
  )
  
  # Iterate over each column to fill
  for (col in cols_to_fill) {
    pb$tick()
    # For each patient
    data <- data %>%
      group_by(!!sym(id_col)) %>%
      arrange(!!sym(date_col)) %>%
      mutate(
        # Create a dataframe of non-NA values
        non_na_dates = list(Date[!is.na(.data[[col]])]),
        non_na_values = list(.data[[col]][!is.na(.data[[col]])])
      ) %>%
      rowwise() %>%
      mutate(
        filled_value = ifelse(
          is.na(.data[[col]]),
          {
            current_date <- !!sym(date_col)
            # Get non-NA dates and values for this patient
            nn_dates <- non_na_dates
            nn_values <- non_na_values
            if (length(nn_dates) == 0) {
              NA
            } else {
              # Calculate the difference in days
              diffs <- abs(as.numeric(difftime(current_date, nn_dates, units = "days")))
              # Find the minimum difference within delay_days
              if (min(diffs) <= delay_days) {
                # Get the index of the closest date
                closest_idx <- which.min(diffs)
                nn_values[closest_idx]
              } else {
                NA
              }
            }
          },
          .data[[col]]
        )
      ) %>%
      ungroup() %>%
      mutate(
        !!sym(col) := filled_value
      ) %>%
      select(-non_na_dates, -non_na_values, -filled_value)
  }
  
  return(data)
}

delay <- 90
columns_to_fill <- c(
  "Alanina_aminotransferasi_alt",
  "Albumina",
  "Altezza",
  "Colesterolo_Hdl",
  "Colesterolo_totale",
  "Creatinina",
  "Distribuzione_di_volume",
  "Ematocrito_hct",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_hb",
  "Emoglobina_massa_media_mch",
  "Eosinofili_perc",
  "Eritrociti_rbc",
  "Ferritina",
  "Ferro_totale",
  "Glucosio",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Monociti_perc",
  "PMAX",
  "Peso",
  "Piastrine",
  "Polso",
  "Proteine_totali",
  "S_alfa_1_globuline",
  "S_alfa_2_globuline",
  "S_beta_1_globuline",
  "S_beta_2_globuline",
  "S_gamma_globuline",
  "Trigliceridi",
  "Volume_medio",
  #  "Alcool",
  #  "Attivita_fisica",
  "Circonferenza_vita"
  #  "Fumo"
)


df_filled <- fill_missing_values(
  data = df_wide,
  cols_to_fill = columns_to_fill,
  id_col = "CAI",
  date_col = "Date",
  delay_days = delay
)

# calculate rows with at least `threshold` missing values
count_rows_with_na <- function(data, columns, threshold = 1) {
  if (!is.data.frame(data)) {
    stop("input not data frame")
  }
  if (!all(columns %in% colnames(data))) {
    stop("columns does not exist")
  }
  
  data %>%
    dplyr::select(all_of(columns)) %>%  # Ensure using dplyr's select
    rowwise() %>%
    mutate(na_count = sum(is.na(c_across(everything())))) %>%  # Count NA in each row
    ungroup() %>%
    summarise(rows_with_na = sum(na_count >= threshold)) %>%  # Count rows meeting the threshold
    pull(rows_with_na)
}

# Calculate for the initial dataset
rows_with_na_initial <- count_rows_with_na(df_wide, columns_to_fill)
rows_with_na_final <- count_rows_with_na(df_filled, columns_to_fill)
rows_with_na_initial
rows_with_na_final
#DATASET WITH NO NA IN GIVEN COLUMNS --------

subset_dataset <- function(df_wide,names){
rows_to_remove = c()
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
return(df)
}

names = c("Colesterolo_Hdl","Circonferenza_vita","Glucosio","PMAX","Trigliceridi")
df_responses_filled <- subset_dataset(df_filled,names)

#PLOT TO SEE HOW MANY NA -----------

# Define the function
plot_na_counts <- function(df, variables, ylim_max) {
  # Create dataset with NA counts
  na_counts <- data.frame(
    Variable = variables,
    NA_Count = sapply(variables, function(var) sum(is.na(df[[var]])))
  )
  
  # Create barplot
  plot <- ggplot(na_counts, aes(x = Variable, y = NA_Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Variable", y = "NA count", title = "Count of NA for each column") +
    ylim(0, ylim_max) +
    theme_minimal()
  
  
  # Print the plot
  print(plot)
}

# Funzione per calcolare e plottare le percentuali di NA
plot_na_percentages <- function(dataset) {
  # Calcola la percentuale di NA per ogni colonna
  na_percentages <- sapply(dataset, function(col) mean(is.na(col)) * 100)
  
  # Crea un dataframe per ggplot
  na_data <- data.frame(
    Column = names(na_percentages),
    NA_Percentage = na_percentages
  )
  
  # Ordina per percentuale decrescente
  na_data <- na_data[order(-na_data$NA_Percentage), ]
  
  # Crea il grafico
  ggplot(na_data, aes(x = reorder(Column, -NA_Percentage), y = NA_Percentage)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Percentuale di Valori Mancanti per Colonna",
      x = "Colonne",
      y = "Percentuale di NA (%)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


plot_na_counts(
  df = df_wide, 
  variables = names, 
  ylim_max = 105000
)
plot_na_percentages(df_responses_filled)
#DATASET WITHOUT ANY NA IN A ROW ----------
columns_to_keep <- c(
  "Alanina_aminotransferasi_alt",
  #"Albumina",
  "Altezza",
  "Colesterolo_Hdl",
  "Colesterolo_totale",
  "Creatinina",
  "Distribuzione_di_volume",
  "Ematocrito_hct",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_hb",
  "Emoglobina_massa_media_mch",
  "Eosinofili_perc",
  "Eritrociti_rbc",
  "Ferritina",
  #"Ferro_totale",
  "Glucosio",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Monociti_perc",
  "PMAX",
  "Peso",
  "Piastrine",
  "Polso",
  "Proteine_totali",
  #"S_alfa_1_globuline",
  #"S_alfa_2_globuline",
  #"S_beta_1_globuline",
  #"S_beta_2_globuline",
  #"S_gamma_globuline",
  "Trigliceridi",
  "Volume_medio",
  #  "Alcool",
  #  "Attivita_fisica",
  "Circonferenza_vita"
  #  "Fumo"
)

df_opt <- subset_dataset(df_responses_filled,columns_to_keep)
df_opt <- df_opt[,columns_to_keep]

saveRDS(df_opt,'filled_dataset.rds')
