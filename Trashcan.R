#SUBSET FILE ---------
##count na after filling ----------
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