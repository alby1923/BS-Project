library(progress)
library(ggplot2)
library(dplyr)

df_wide <- readRDS('../datasets/flattened_dataset.rds')

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
df_responses <- subset_dataset(df_wide,names)

#FILLER OF NA IN DATASET ----------------

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

delay <- 180
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
    stop("The input 'data' must be a data frame or tibble.")
  }
  if (!all(columns %in% colnames(data))) {
    stop("Some columns in 'columns' are not present in the data frame.")
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

saveRDS(df_filled, file = "df_filled.rds")
