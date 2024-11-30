library(progress)
library(ggplot2)
library(dplyr)

df_wide <- readRDS('flattened_dataset.rds')

#ADD FACTORS ---------
df_to_add <- readRDS('cleaned_dataset.rds')
#to modify in numerical
#CHANGE FACTORS TO NUMERICAL VARIABLES
# Define mappings for each factor column
map_rh <- function(x) {
  as.numeric(factor(x, levels = c("NEG", "POS"), labels = c(1, 2)))
}

map_AB0 <- function(x) {
  as.numeric(factor(x, levels = c("0", "A", 'AB', 'B'), labels = c(1, 2, 3, 4)))
}

map_SESSO <- function(x) {
  as.numeric(factor(x, levels = c("1", "2"), labels = c(1, 2)))
}

# Apply mappings to the respective columns
df_to_add$Rh <- lapply(df_to_add$Rh, map_rh)
df_to_add$AB0 <- lapply(df_to_add$AB0, map_AB0)
df_to_add$SESSO <- lapply(df_to_add$SESSO, map_SESSO)

# Convert lists to vectors if necessary (if the structure requires)
df_to_add$Rh <- sapply(df_to_add$Rh, unlist)
df_to_add$AB0 <- sapply(df_to_add$AB0 , unlist)
df_to_add$SESSO <- sapply(df_to_add$SESSO, unlist)
  
# Merge based on the patient ID column
df_to_add_subset <- df_to_add[, c('CAI','Rh','AB0','SESSO'), drop = FALSE]
df_wide <- merge(df_wide, df_to_add_subset, by = 'CAI', all.x = TRUE)

#STANDARDIZE COVARIATES -----------
#needed to fix a little bug after flattening dataset

adjust_round <- function(vec) {
  for (i in seq_along(vec)) {
    if (!is.na(vec[i]) && vec[i] %% 1 != 0) {  #check value if not integer
      # find preceding line not na
      prev_value <- NA
      for (j in seq(i - 1, 1, -1)) {
        if (!is.na(vec[j])) {
          prev_value <- vec[j]
          break
        }
      }
      
      #round based on condition
      if (!is.na(prev_value)) {
        if (prev_value < vec[i]) {
          vec[i] <- ceiling(vec[i])  # ceiling
        } else {
          vec[i] <- floor(vec[i])   # floor
        }
      }
    }
  }
  return(vec)
}

df_wide$Alcool <- adjust_round(df_wide$Alcool)
df_wide$Fumo <- adjust_round(df_wide$Fumo)
df_wide$Attivita_fisica <- adjust_round(df_wide$Attivita_fisica)


for(j in 1:ncol(df_wide)){
  if(!(j %in% c(1,2,6,18,22,32,36))){
    df_wide[,j] <- scale(df_wide[,j],center = TRUE, scale = TRUE)
  }
}
#FILLER OF NA IN DATASET ----------------

library(lubridate) # For handling date differences
library(data.table)  # Added for efficient data manipulation

# Ensure the Date column is in Date format
df_wide <- df_wide %>%
  mutate(Date = as.Date(Date))
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
              # Calculate the difference in days with previous and future dates
              diffs_previous <- abs(as.numeric(difftime(current_date, nn_dates, units = "days")))
              diffs_future <- abs(as.numeric(difftime(nn_dates, current_date, units = "days")))
              
              # Combine both differences (previous and future)
              diffs_combined <- c(diffs_previous, diffs_future)
              nn_combined_values <- c(nn_values, nn_values)
              
              # Find the minimum difference within delay_days
              if (min(diffs_combined) <= delay_days) {
                # Get the index of the closest date (either previous or future)
                closest_idx <- which.min(diffs_combined)
                nn_combined_values[closest_idx]
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

delay <- 15
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
    "Alcool",
   "Attivita_fisica",
  "Circonferenza_vita",
    "Fumo"
)

df_filled <- fill_missing_values(
  data = df_wide,
  cols_to_fill = columns_to_fill,
  id_col = "CAI",
  date_col = "Date",
  delay_days = delay
)
#DATASET WITH NO NA IN GIVEN COLUMNS --------

subset_dataset <- function(df_wide,names){
rows_to_remove = c()
target_cols = which(colnames(df_wide) %in% names)
pb <- progress_bar$new(total = nrow(df_wide)) #just to visualize how much time remaining
for(i in 1:nrow(df_wide)){
  flag <- 0
  pb$tick()
  for(col in target_cols){
    if(col == 'Date'){
      next
    }
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
  
  # Barplot
  plot <- ggplot(na_counts, aes(x = reorder(Variable, NA_Count), y = NA_Count)) +
    geom_bar(stat = "identity", aes(fill = NA_Count), width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = NA_Count), vjust = -0.5, size = 4, color = "black") +
    scale_fill_gradient(low = "skyblue", high = "darkblue") +
    labs(
      x = "Variable",
      y = "Count of NA"
    ) +
    ylim(0, ylim_max) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(),  # Rimuove le linee di griglia principali
      panel.grid.minor = element_blank(),  # Rimuove le linee di griglia minori
      panel.background = element_blank(),  # Imposta lo sfondo del pannello su bianco
      plot.background = element_blank()    # Imposta lo sfondo del grafico su bianco
    )
  
  ggsave("NAresponses.pdf", plot, width = 10, height = 8, dpi = 300)
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
  
  na_data$Highlight <- ifelse(na_data$NA_Percentage > 5, "Remove", "Keep")
  
  # Barplot
  plot <- ggplot(na_data, aes(x = reorder(Column, -NA_Percentage), y = NA_Percentage)) +
    geom_bar(
      stat = "identity",
      aes(fill = Highlight),
      width = .85,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("Remove" = "red", "Keep" = "steelblue")) +
    geom_text(
      aes(label = paste0(round(NA_Percentage, 1), "%")),
      hjust = -0.1,
      size = 3.5,
      color = "black"
    ) +
    coord_flip() +
    labs(
      x = "Columns",
      y = "Percentage of NA (%)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(size = 12, color = "darkgray"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(),  # Rimuovi le linee di griglia orizzontali
      panel.grid.minor = element_blank(),  # Rimuovi le linee di griglia minori
      panel.background = element_blank(),  # Imposta lo sfondo del pannello su bianco
      plot.background = element_blank()    # Imposta lo sfondo generale su bianco
    ) +
    ylim(0, max(na_data$NA_Percentage) + 5)  
  
  print(plot)
  
  ggsave("NApercentage.pdf", plot, width = 10, height = 8, dpi = 300)
}

plot_na_counts(
  df = df_wide, 
  variables = names, 
  ylim_max = 105000
)
plot_na_percentages(df_responses_filled)
#DATASET WITHOUT ANY NA IN A ROW ----------
columns_to_keep <- c(
  "CAI",
  "Date",
  #"Alanina_aminotransferasi_alt",
  #"Albumina",
  "Altezza",
  "Colesterolo_Hdl",
  "Colesterolo_totale",
  #"Creatinina",
  "Distribuzione_di_volume",
  "Ematocrito_hct",
  "Emoglobina_conc_media_mchc",
  "Emoglobina_hb",
  "Emoglobina_massa_media_mch",
  "Eosinofili_perc",
  "Eritrociti_rbc",
  #"Ferritina",
  #"Ferro_totale",
  "Glucosio",
  "Leucociti_wbc",
  "Linfociti_perc",
  "Monociti_perc",
  "PMAX",
  "Peso",
  "Piastrine",
  "Polso",
  #"Proteine_totali",
  #"S_alfa_1_globuline",
  #"S_alfa_2_globuline",
  #"S_beta_1_globuline",
  #"S_beta_2_globuline",
  #"S_gamma_globuline",
  "Trigliceridi",
  "Volume_medio",
  "Alcool",
  "Attivita_fisica",
  "Circonferenza_vita",
  "Fumo",
  'Rh',
  'SESSO',
  'AB0'
)

df_opt1 <- subset_dataset(df_responses_filled,columns_to_keep)
df_opt1 <- df_opt1[,columns_to_keep]

#DELTA DATE COLUMN -------------
#check all dates to be after a certain date
reference_date <- as.Date("2019-01-01")

# counts how many dates are before the reference date
dates_before <- sum(df_opt1$Date < reference_date, na.rm = TRUE)
dates_before
dates <- df_opt1$Date[df_opt1$Date < reference_date]
dates

#create delta date column for timestamps
calculate_delta_date <- function(df, patient_col, date_col) {
  
  df$delta_date <- unlist(by(
    df, 
    df[[patient_col]], 
    function(sub_df) {
      
      first_date <- sub_df[[date_col]][1]
     
      delta <- as.numeric(difftime(sub_df[[date_col]], first_date, units = "days"))
      return(delta)
    }
  ))
  
  return(df)
}

df_opt1 <- df_opt1[order(df_opt1$CAI), ]
df_opt1 <- calculate_delta_date(df_opt1, "CAI", "Date")
save <- df_opt1

#add age column
ages <- readRDS('final_dataset.rds')
ages <- ages %>%
  distinct(CAI, .keep_all = TRUE)
df_opt1 <- df_opt1 %>%
  left_join(ages %>% select(CAI, DATA_NASCITA), by = "CAI")
df_opt1 <- df_opt1 %>%
  mutate(eta = as.numeric(difftime(Date, DATA_NASCITA, units = "days")) / 365) %>%
  select(-DATA_NASCITA)

df_opt1$eta_std <- scale(df_opt1$eta, center = TRUE, scale = TRUE)

saveRDS(df_opt1,'filled_dataset.rds')
