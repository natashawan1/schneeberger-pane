## Title: Correlation Map
## Last update: 11/18/2025
## Written by: Natasha Wan
## Lab: Schneeberger-Pané
## ---------------------------------------------------------------------------

# Some set up...

getwd() # working directory
dir()   # files in  working directory
rm(list = ls()) # clean up environment
setwd("~/Desktop/Marc Lab Work/Female:Male Social Project/Social Homeostasis Data/Raw Data") # set working directory

library(dplyr)
library(stringr)
library(purrr)
library(readxl)
library(scales)
library(patchwork)
library(ggplot2)


# Note: maybe break up functions

# Define input and output files
input_files_14d <- c("14d_control_isolated_Fem.csv", "14d_control_overcrowded_Fem.csv",
                 "14d_control_isolated_Male.csv", "14d_control_overcrowded_Male.csv")
output_files_14d <- c("14d_control_isolated_Fem_processed.csv", 
                  "14d_control_overcrowded_Fem_processed.csv",
                  "14d_control_isolated_Male_processed.csv", "14d_control_overcrowded_Male_processed.csv")


input_files_48hr <- c("48hr_control_isolated_Fem.csv", "48hr_control_overcrowded_Fem.csv",
                    "48hr_control_isolated_Male.csv", "48hr_control_overcrowded_Male.csv")
output_files_48hr <- c("48hr_control_isolated_Fem_processed.csv", "48hr_control_overcrowded_Fem_processed.csv",
                       "48hr_control_isolated_Male_processed.csv", "48hr_control_overcrowded_Male_processed.csv")

processed_files <- function(input_files, output_files) {
  
  if (length(input_files) != length(output_files)) {
    stop("Input and output files must have the same length.")
  }
  
  for (i in seq_along(input_files)) {
    # Read the data
    data <- read.csv(input_files[i])
    
    # Keep only rows where hemisphere != 0
    data <- data[data$hemisphere != 0, ]
    
    # Remove non-brain region rows
    values_to_remove <- c("No label", "universe")
    data <- data[!(data$name %in% values_to_remove), ]
    
    # Create a new data frame to store results
    df <- data.frame(
      name = data$name,
      p_values = data$p_value,
      q_value = data$q_value
    )

    # Determine which columns exist for condition means and control means
    condition_cols <- c("mean_isolated", "mean_overcrowded", "mean_isolated_female", "mean_overcrowded_female", "mean_overcrowed_female")
    control_cols <- c("mean_control", "mean_control_female")

    existing_condition <- intersect(names(data), condition_cols)
    existing_control <- intersect(names(data), control_cols)

    if (length(existing_condition) == 0 | length(existing_control) == 0) {
      stop("No matching condition or control columns found in ", input_files[i])
    }

    # Take the first existing column if multiple match (you can change logic if needed)
    condition_mean <- data[[existing_condition[1]]]
    control_mean <- data[[existing_control[1]]]

    # Compute log2 fold change
    df$log2FoldChange <- log2(condition_mean / control_mean)

    # Compute direction: positive if condition > control, negative otherwise
    df$direction <- ifelse(condition_mean > control_mean, "up", "down")

    # Order the df by log2FoldChange descending 
    df <- df[order(df$log2FoldChange, decreasing = TRUE), ]
   
    # Save the processed data
    write.csv(df, output_files[i], row.names = FALSE)

    message("Processed: ", input_files[i], " → ", output_files[i]) }}

processed_files(input_files_48hr, output_files_48hr)
processed_files(input_files_14d, output_files_14d)


# Function to read each file and extract metadata
read_with_metadata <- function(filename) {
  data <- read.csv(filename)
  
  # Extract Condition (text between "control_" and last "_F"/"_M"/"_Fem"/"_Male")
  condition <- str_extract(filename, "(?<=control_).*?(?=_[FM]|_Fem|_Male)")
  
  # Extract Gender
  gender <- case_when(
    str_detect(filename, "_F") | str_detect(filename, "Fem") ~ "F",
    str_detect(filename, "_M") | str_detect(filename, "Male") ~ "M",
    TRUE ~ NA_character_
  )
  
  # Add new columns
  data <- data %>%
    mutate(
      Condition = condition,
      Gender = gender
    )
  
  return(data)
}

# Combine all files into one big data frame
combined_data_48hr <- map_dfr(output_files_48hr, read_with_metadata)
combined_data_14d <- map_dfr(output_files_14d, read_with_metadata)


# Adding abbreviations to the brain regions
ref <- read_excel("Allen Brain Brain Regions (names and abbriviations).xlsx", col_names = FALSE)

# Assign the second row as headers
colnames(ref) <- as.character(ref[2, ])
ref <- ref[-c(1,2), ]  # remove first three rows
colnames(ref)

ref <- ref %>%
  select(
    name = `full structure name`,
    abbreviation,
    level = `depth in tree`
  )

biased_circuit <- c("Ventral tegmental area", "Nucleus accumbens",
                    # "Cortical amygdalar area", "Lateral amygdalar nucleus", "Basolateral amygdalar nucleus", "Basomedial amygdalar nucleus", "Posterior amygdalar nucleus", "Striatum-like amygdalar nuclei",
                    # "Somatosensory areas",
                    "Hypothalamus",
                    "Field CA1", "Field CA2", "Field CA3")


filtered_data_48hr <- combined_data_48hr %>%
  group_by(name) %>%
  filter(name %in% biased_circuit |
           str_detect(name, "^Dentate gyrus"))

filtered_data_14d <- combined_data_14d %>%
  group_by(name)  %>%
  filter(name %in% biased_circuit |
           str_detect(name, "^Dentate gyrus"))

common_nuclei <- intersect(filtered_data_48hr$name, filtered_data_14d$name)
filtered_data_48hr <- filtered_data_48hr %>% filter(name %in% common_nuclei)
filtered_data_14d <- filtered_data_14d %>% filter(name %in% common_nuclei)

abbrevlevel_data_48hr <- filtered_data_48hr %>%
  left_join(ref, by = "name")  # Adds abbreviation and level

abbrevlevel_data_14d <- filtered_data_14d %>%
  left_join(ref, by = "name")  # Adds abbreviation and level

all48hr <-abbrevlevel_data_48hr
all14d <- abbrevlevel_data_14d

combined <- full_join(all48hr,all14d)

## Graphing

# set color scale mean the same thing in every heatmap
lims <- max(
  abs(c(all48hr$log2FoldChange, all14d$log2FoldChange)),
  na.rm = TRUE
)

graph <- function(data, title) {
  
    data <- data %>%
      mutate(
        category = paste(Gender, Condition),
        category = factor(
          category,
          levels = c("F isolated", "F overcrowded", "M isolated", "M overcrowded")
        ),
        significance_symbol = ifelse(p_values < 0.05, "*", ""),
        nuclei = paste0(abbreviation, " — ", name),
        nuclei = factor(nuclei, levels = unique(nuclei))
      )
    
  ggplot(data, aes(x = category, y = nuclei, fill = log2FoldChange)) +
    geom_tile(color = "white", width = 1, height = 1) +
    geom_text(aes(label = significance_symbol),
              size = 5, color = "white", vjust = 0.8) +
    scale_fill_gradient2(
      low = "#4575B5",
      mid = "white",
      high = "#D73127",
      midpoint = 0,
      limits = c(-lims,lims),
      name = "log2 fold change"
    ) +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 75, hjust = 1),
      legend.position = "right"
    ) +
    labs(title = title)
}


p1 <- graph(data = all48hr, title = "48 hour")
p2 <- graph(data = all14d,  title = "14 days")

combined_plot <- p1 + p2 +
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "right")
print(combined_plot)

