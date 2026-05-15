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

## READ ME: This code first 1) processes all of the files from ClearMap to isolate only half of the brain, 2) remove brain regions that have no labels,
## extract condition and control means, and 3) compute and find the directionality of the log2foldchange.



# Define input and output files

input_files_fem_48hr <- c("48hr_control_isolated_Fem.csv", "48hr_control_overcrowded_Fem.csv")
output_files_fem_48hr <- c("48hr_control_isolated_Fem_processed.csv", "48hr_control_overcrowded_Fem_processed.csv")

input_files_fem_14d <- c("14d_control_isolated_Fem.csv", "14d_control_overcrowded_Fem.csv")
output_files_fem_14d <- c("14d_control_isolated_Fem_processed.csv", "14d_control_overcrowded_Fem_processed.csv")

input_files_male_48hr <- c("48hr_control_isolated_Male.csv", "48hr_control_overcrowded_Male.csv")
output_files_male_48hr <- c("48hr_control_isolated_Male_processed.csv", "48hr_control_overcrowded_Male_processed.csv")

input_files_male_14d <- c("14d_control_isolated_Male.csv", "14d_control_overcrowded_Male.csv")
output_files_male_14d <- c("14d_control_isolated_Male_processed.csv", "14d_control_overcrowded_Male_processed.csv")


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
    # df$direction <- ifelse(df$log2FoldChange > 0, "up", "down")
    
    # Save the processed data
    write.csv(df, output_files[i], row.names = FALSE)
    
    message("Processed: ", input_files[i], " → ", output_files[i]) }}

processed_files(input_files_fem_48hr, output_files_fem_48hr)
processed_files(input_files_fem_14d, output_files_fem_14d)
processed_files(input_files_male_48hr, output_files_male_48hr)
processed_files(input_files_male_14d, output_files_male_14d)


## READ ME: Next we are working to 1) extract the metadata for sex and condition and 2) add in the abbreviations from the Allen Brain Atlas,
## depth in its "tree" (larger numbers are more specific).

# Function to read each file and extract metadata
read_with_metadata <- function(filename) {
  data <- read.csv(filename)
  
  # Extract Condition (text between "control_" and last "_F"/"_M"/"_Fem"/"_Male")
  condition <- str_extract(filename, "(?<=control_).*?(?=_[FM]|_Fem|_Male)")
  
  # # Extract Time
  # time <- case_when(
  #   str_detect(filename, "48hr") ~ "48hr",
  #   str_detect(filename, "14d")  ~ "14d",
  #   TRUE ~ NA_character_
  # )
  
  # Add new columns
  data <- data %>%
    mutate(Condition = condition)
  
  return(data)
}

# Combine all files into one big data frame
combined_data_fem_48hr <- map_dfr(output_files_fem_48hr, read_with_metadata)
combined_data_fem_14d <- map_dfr(output_files_fem_14d, read_with_metadata)
combined_data_male_48hr <- map_dfr(output_files_male_48hr, read_with_metadata)
combined_data_male_14d <- map_dfr(output_files_male_14d, read_with_metadata)


# Adding abbriviations to the brain regions
ref <- read_excel("Allen Brain Brain Regions (names and abbriviations).xlsx", col_names = FALSE)

# Assign the second row as headers
colnames(ref) <- as.character(ref[2, ])
ref <- ref[-c(1,2), ]  # remove first two rows

ref <- ref %>%
  rename ( name = 'full structure name',
           abbreviation = abbreviation,
           level = 'depth in tree'
  ) %>%
  select(name, abbreviation, level)  


## Next, the data sets are combined and filtered for brain regions that appear in all data sets.
## There is an unused code for filtering out brain regions of a certain p-value. Currently this is also filtering for layer 6.


filter_consistent <- function(df, other_df, p_cutoff = 0.05) {
  df %>%
    group_by(name) %>%
    filter(n() > 1) %>%
    filter(all(direction == "up") | all(direction == "down")) %>%
    filter(sum(p_values < p_cutoff) >= 0) %>%   # at least x significant p-values
    filter(name %in% unique(other_df$name)) %>%
    ungroup()
}

filtered_data_48hr_f  <- filter_consistent(combined_data_fem_48hr,  combined_data_fem_14d)
filtered_data_14d_f <- filter_consistent(combined_data_fem_14d, combined_data_fem_48hr)

common_nuclei_f <- intersect(filtered_data_48hr_f$name,
                             filtered_data_14d_f$name)

filtered_data_48hr_f  <- filtered_data_48hr_f  %>% filter(name %in% common_nuclei_f)
filtered_data_14d_f <- filtered_data_14d_f %>% filter(name %in% common_nuclei_f)

# Join iso + over by nucleus to check directions
dircheck_f <- filtered_data_48hr_f %>%
  select(name, direction_48h = direction) %>%
  distinct() %>%
  inner_join(
    filtered_data_14d_f %>%
      select(name, direction_14d = direction) %>%
      distinct(),
    by = "name"
  ) %>%
  # keep only direction *mismatch* nuclei
  filter(direction_48h != direction_14d)

# valid mismatching nuclei
valid_nuclei_f <- dircheck_f$name

# keep only those in both datasets
filtered_data_48hr_f  <- filtered_data_48hr_f  %>% filter(name %in% valid_nuclei_f)
filtered_data_14d_f <- filtered_data_14d_f %>% filter(name %in% valid_nuclei_f)

filtered_data_48hr_m  <- filter_consistent(combined_data_male_48hr,  combined_data_male_14d)
filtered_data_14d_m <- filter_consistent(combined_data_male_14d, combined_data_male_48hr)

common_nuclei_m <- intersect(filtered_data_48hr_m$name,
                             filtered_data_14d_m$name)

filtered_data_48hr_m  <- filtered_data_48hr_m  %>% filter(name %in% common_nuclei_m)
filtered_data_14d_m <- filtered_data_14d_m %>% filter(name %in% common_nuclei_m)

dircheck_m <- combined_data_male_48hr %>%
  select(name, direction_48h = direction) %>%
  distinct() %>%
  inner_join(
    filtered_data_14d_m %>%
      select(name, direction_14d = direction) %>%
      distinct(),
    by = "name"
  ) %>%
  filter(direction_48h != direction_14d)

valid_nuclei_m <- dircheck_m$name

filtered_data_48hr_m  <- filtered_data_48hr_m  %>% filter(name %in% valid_nuclei_m)
filtered_data_14d_m <- filtered_data_14d_m %>% filter(name %in% valid_nuclei_m)

add_abbrev <- function(df) {
  df %>% left_join(ref, by = "name") %>%
    filter(level %in% c(6))
  
}

abbrevlevel_data_48h_f  <- add_abbrev(filtered_data_48hr_f)
abbrevlevel_data_14d_f <- add_abbrev(filtered_data_14d_f)
abbrevlevel_data_48h_m  <- add_abbrev(filtered_data_48hr_m)
abbrevlevel_data_14d_m <- add_abbrev(filtered_data_14d_m)

## READ ME: Now, the heatmap's color is being set up to find the max and min log2FoldChange values. The "lims" were hard coded for
## Natasha's thesis so mulitple heatmaps could have the same color gradient and be compared visually to each other.

all48h_f <-abbrevlevel_data_48h_f
all14d_f <- abbrevlevel_data_14d_f

all48h_m <-abbrevlevel_data_48h_m
all14d_m <- abbrevlevel_data_14d_m

combined_f <- full_join(all48h_f,all14d_f)
combined_m <- full_join(all48h_m,all14d_m)
combined <- full_join (combined_f,combined_m)



# Compute gradient limits
max_val <- max(combined$log2FoldChange, na.rm = TRUE)
min_val <- min(combined$log2FoldChange, na.rm = TRUE)
lims <- c(min_val, max_val)
lims <- c(-5.330374, 4.169925) # hard coded for thesis

print(lims)  # check --> -3.341037  4.169925

## READ ME: Next came graphing where labels were first combined so the main observations were
## the temporal changes unique to either sex

graph <- function(data, nuclei, brain_change, title) {
  
  # Define condition category (combining Gender and Time)
  data <- data %>%
    mutate(
      category = paste(Condition),
      category = factor(
        category,
        levels = c("isolated","overcrowded")
      ),
      significance_symbol = ifelse(p_values < 0.05, "*", ""),
      nuclei = abbreviation
    )
  
  # Create the heatmap
  p <- ggplot(data, aes(x = category, y = nuclei, fill = log2FoldChange)) +
    geom_tile(color = "black", width = 0.6, height = 0.6) +
    geom_text(aes(label = significance_symbol),
              size = 5, color = "white", vjust = 0.8) +
    scale_fill_gradientn(
      colors = c("#4575B5", "white", "#D73127"),
      values = rescale(c(min_val, 0, max_val)),
      limits = c(min_val, max_val),
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



p1 <- graph (data = all48h_f, nuclei = all48h_f$abbreviation, brain_change = all48h_f$direction, title = "Females 48hr")
p2 <- graph (data = all14d_f, nuclei = all14d_f$abbreviation, brain_change = all14d_f$direction, title = "Females 14d")
p3 <- graph (data = all48h_m, nuclei = all48h_m$abbreviation, brain_change = all48h_m$direction, title = "Males 48hr")
p4 <- graph (data = all14d_m, nuclei = all14d_m$abbreviation, brain_change = all14d_m$direction, title = "Males 14d")

combined_plot_1 <- p1 + p2 + plot_layout(ncol = 2) + labs(fill = "brain_change")
combined_plot_2 <- p3 + p4 + plot_layout(ncol = 2) + labs(fill = "brain_change")
print(combined_plot_1)
print(combined_plot_2)

## READ ME: These heatmaps were combined to visualize the regions all together.

# Overlap and unique (Looking at Time and Gender)
# all48h_f, all14d_f, all48h_m, all14d_m already exist

# Female regions
female_regions <- union(all48h_f$name, all14d_f$name)

# Male regions
male_regions <- union(all48h_m$name, all14d_m$name)

# Non-overlapping sets
female_only_regions <- setdiff(female_regions, male_regions)
male_only_regions   <- setdiff(male_regions,  female_regions)

h48_f <- all48h_f  %>% mutate(Sex="Female", Time="48h")
h14_f <- all14d_f %>% mutate(Sex="Female", Time="14d")

h48_m <- all48h_m  %>% mutate(Sex="Male", Time="48h")
h14_m <- all14d_m %>% mutate(Sex="Male", Time="14d")

combined_all <- bind_rows(
  h48_f, h14_f,
  h48_m, h14_m
)

sig_counts <- combined_all %>%
  group_by(name, abbreviation, Sex, Time) %>%
  summarise(n_sig = sum(p_values < 0.05, na.rm = TRUE), .groups = "drop")

non_overlap <- combined_all %>%
  filter(
    name %in% female_only_regions |
      name %in% male_only_regions
  ) %>%
  mutate(
    SexCategory = case_when(
      name %in% female_only_regions ~ "Female Only",
      name %in% male_only_regions ~ "Male Only"
    )
  )

non_overlap2 <- non_overlap %>%
  left_join(sig_counts, by = c("name", "abbreviation", "Sex", "Time")) %>%
  mutate(
    sig_label = case_when(
      n_sig == 0 ~ "",
      n_sig == 1 ~ "*",
      n_sig == 2 ~ "**",
      n_sig >= 3 ~ "***"
    )
  )

plot_nonoverlap_heatmap <- function(data, title) {
  
  data <- data %>%
    mutate(
      Time = factor(Time, levels = c("48h", "14d"))
    )
  
  ggplot(data, aes(x = Time, y = abbreviation, fill = log2FoldChange)) +
    geom_tile(color = "black", width = 0.7, height = 0.7) +
    geom_text(aes(label = sig_label),
              color = "white", size = 5, vjust = 0.7) +
    scale_fill_gradientn(
      colors = c("#4575B5", "white", "#D73127"),
      values = rescale(c(min_val, 0, max_val)),
      limits = c(min_val, max_val),
      name = "log2 fold change"
    ) +
    facet_grid(SexCategory ~ Sex, scales="free_y", space="free_y") +
    theme_minimal(base_size=12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust=1, size=10),
      strip.text = element_text(size=11, face="bold"),
      legend.position = "none"
    ) +
    labs(title = title)
}


p <- ggplot(data, aes(x = category, y = nuclei, fill = log2FoldChange)) +
  geom_tile(color = "black", width = 0.6, height = 0.6) +
  geom_text(aes(label = significance_symbol),
            size = 5, color = "white", vjust = 0.8) +
  scale_fill_gradientn(
    colors = c("#4575B5", "white", "#D73127"),
    values = rescale(c(min_val, 0, max_val)),
    limits = c(min_val, max_val),
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


non_overlap_plot <- plot_nonoverlap_heatmap(
  non_overlap2,
  title = "Brain Regions Unique to One Sex (Non-Overlapping)"
)

print(non_overlap_plot)

