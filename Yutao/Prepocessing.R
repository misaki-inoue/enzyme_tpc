# Load necessary libraries
library(dplyr)
library(stringr)
library(MASS)
library(ggplot2)
setwd("~/Documents/FinalProject")

# Read the datasets
thermRespData <- read.csv("data/ThermRespData.csv")
fungalData <- read.csv("data/Fungaldata.csv")

# Inspect the structure of the datasets
str(thermRespData)


# Check column names to confirm the ID column and the names of other columns
print(colnames(thermRespData))
#create subset
all_subset_therm <- split(thermRespData, thermRespData$ID)
#select possible useful columns therm
# Load the dplyr package
library(dplyr)

# Ensure you are using the correct column names
selected_subset_therm <- thermRespData %>%
  dplyr::select(trait_value = OriginalTraitValue , OriginalTraitUnit, temp = ConTemp, Unit = ConTempUnit, ID)



fungalData <- fungalData %>%
  mutate(
    ID = as.numeric(str_extract(OriginalID, "[0-9]+")) + 904  # Extract, convert, and add in one step
  )

# Validate the updated IDs
unique(fungalData$ID)
# Adding a check for NA in IDs after conversion
sum(is.na(fungalData$ID))
# Filtering out rows with NA IDs if necessary
fungalData <- fungalData %>% filter(!is.na(ID))

# Alternatively, handle or flag entries without a numeric part
fungalData <- fungalData %>%
  mutate(
    ID = if_else(is.na(ID), -1, ID)  # Assign a special value to IDs that couldn't be processed
  )

# Select specific columns and rename them as needed
selected_subset_fungal <- fungalData %>%
  dplyr::select(
    trait_value = OriginalTraitValue, 
    OriginalTraitUnit, 
    temp = Interactor1Temp, 
    Unit = Interactor1TempUnit, 
    ID
  )
# Check the transformation
print(head(fungalData$ID))

# Checking the updated IDs to confirm the change
print(head(selected_subset_fungal$ID))
#merge together
# Full join to include all records from both subsets and match where possible
selected_subset <- full_join(selected_subset_fungal, selected_subset_therm, by = "ID")
# Combining rows from both subsets
selected_subset <- bind_rows(selected_subset_fungal, selected_subset_therm)
# Check the structure of the combined dataset
str(selected_subset)

# View the first few rows to confirm
head(selected_subset)

##########Prepocessing
# Filter out non-positive trait measurements
cleaned_data <- selected_subset %>%
  filter(trait_value > 0)

# Viewing the first few rows of the cleaned dataset to ensure correct filtering
head(cleaned_data)

# Checking the minimum value of 'TraitValue' to ensure no non-positive values remain
min(cleaned_data$trait_value)

# Na check
cleaned_data <- selected_subset %>%
  filter(!is.na(trait_value) & trait_value > 0)

# Split subset 
cleaned_data_subset <- split(cleaned_data, cleaned_data$ID)


# Check if 'Temp' exists in all subsets
all(sapply(cleaned_data_subset, function(data) "temp" %in% names(data)))

# Initialize a list to hold data subsets with non-negative slopes
relevant_data_subsets <- list()

# Fit a linear model to each subset and collect those with non-negative slopes
for (subset_name in names(cleaned_data_subset)) {
  data <- cleaned_data_subset[[subset_name]]
  if ("temp" %in% names(data)) {
    model <- lm(trait_value ~ temp, data = data)
    coef_summary <- summary(model)$coefficients
    # Check if Temp's coefficient is non-negative and collect the data subset
    if (!is.null(coef_summary) && "temp" %in% rownames(coef_summary) && coef_summary["temp", "Estimate"] >= 0) {
      relevant_data_subsets[[subset_name]] <- data
    }
  }
}
# Check each subset for NA or infinite values
check_data <- function(df) {
  sum(is.na(df$trait_value)) + sum(is.na(df$temp)) +
    sum(is.infinite(df$trait_value)) + sum(is.infinite(df$temp))
}

na_inf_counts <- lapply(relevant_data_subsets, check_data)
print(na_inf_counts)  # This will show you if any subset contains NA or infinite values

# Cleaning data by removing rows with NA or infinite values
clean_data <- function(df) {
  df <- df[!is.na(df$trait_value) & !is.na(df$temp) &
             !is.infinite(df$trait_value) & !is.infinite(df$temp), ]
  return(df)
}

cleaned_subsets <- lapply(relevant_data_subsets, function(df) {
  if("ID" %in% names(df)) {
    names(df)[names(df)=="ID"] <- "id"
  }
  return(df)
})

# Assume subset "84" starts at the 85th position
final_subsets <- cleaned_subsets[503:length(cleaned_subsets)]








