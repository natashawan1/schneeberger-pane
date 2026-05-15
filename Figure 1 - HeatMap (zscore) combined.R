## ---------------------------------------------------------------------------

## Heat Maps - Male/Female 48hr iso/over
## Last update: 3-24-2025
## Written by: Natasha Wan
## Adapted from: sanbomics_scripts
## Lab: Schneeberger-Pané

## ---------------------------------------------------------------------------

# Some set up...

rm(list = ls()) # clean up environment

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(tibble)
library(grid)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(stringr)


# READ ME: This script takes in data which has been normalized (z-scores) to allow for
# comparisons between experiments taken at different time points. 

load_and_clean <- function(file_path) {
  df <- read.csv(file_path, row.names = 1)
  
  df_clean <- df %>%
    tidyr::drop_na() %>%
    dplyr::select(OM_mean, OF_mean, IM_mean, IF_mean) %>%
    rownames_to_column("name") %>%
    # Remove everything after comma or the word "layer" to collapse sub-regions
    mutate(name = str_remove(name, ",.*$"),
           name = str_remove(name, " layer.*$"),
           name = str_trim(name)) %>%
    group_by(name) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames("name")
  
  # Scale values for Z-score comparison
  df_scaled <- as.data.frame(lapply(df_clean, function(x) as.numeric(scale(x))))
  rownames(df_scaled) <- rownames(df_clean)
  return(df_scaled)
}

### Part 1: MAPPING TO ATLAS & HIEARCHY -------------------------------------------------------------------------

process_with_atlas <- function(df, AllenBrainFile, levels_to_keep, region_list) {
  
  ref <- read_excel(AllenBrainFile, col_names = FALSE)
  colnames(ref) <- as.character(ref[2, ])
  ref <- ref[-c(1,2), ]
  
  # YOUR BRAIN ORDER & HIERARCHY LOGIC
  brain_order <- c("Cerebral cortex", "Hippocampus", "Striatum", "Pallidum", 
                   "Thalamus", "Hypothalamus", "Interbrain", "Midbrain", 
                   "Pons", "Medulla", "Hindbrain", "Cerebellum", "Other")
  
  ref <- ref %>%
    select(name = `full structure name`, abbreviation, level = `depth in tree`,
           `structure ID`, structure_id_path) %>%
    mutate(level = as.numeric(level), name = str_trim(name)) %>%
    mutate(parent_name = case_when(
      str_detect(structure_id_path, "/997/8/567/688/695/1089/") ~ "Hippocampus",
      str_detect(structure_id_path, "/997/8/567/688/") ~ "Cerebral cortex",
      str_detect(structure_id_path, "/997/8/567/623/477/") ~ "Striatum",
      str_detect(structure_id_path, "/997/8/567/623/803/") ~ "Pallidum",
      str_detect(structure_id_path, "/997/8/343/1129/549/") ~ "Thalamus", # Specific first
      str_detect(structure_id_path, "/997/8/343/1129/1097/") ~ "Hypothalamus",
      str_detect(structure_id_path, "/997/8/343/1129/") ~ "Interbrain",   # General second
      str_detect(structure_id_path, "/997/8/343/313/") ~ "Midbrain",
      str_detect(structure_id_path, "/997/8/343/1065/771/") ~ "Pons",
      str_detect(structure_id_path, "/997/8/343/1065/354/") ~ "Medulla",
      str_detect(structure_id_path, "/997/8/343/1065/") ~ "Hindbrain",
      str_detect(structure_id_path, "/997/8/512/") ~ "Cerebellum",
      TRUE ~ "Other"
    )) %>%
    filter(level %in% levels_to_keep) %>%
    mutate(parent_name = factor(parent_name, levels = brain_order))
  
  # Filter data for ONLY the regions in our target list (the common regions)
  final_names <- ref %>%
    filter(name %in% region_list) %>%
    arrange(parent_name, structure_id_path) %>%
    pull(name)
  
  mat_heatmap <- df[final_names, ]
  
  # Map Abbreviations for display
  name_to_abbrev <- ref %>% select(name, abbreviation)
  mat_display <- mat_heatmap %>%
    as.data.frame() %>%
    rownames_to_column("name") %>%
    left_join(name_to_abbrev, by = "name") %>%
    mutate(name = abbreviation) %>%
    select(-abbreviation) %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  return(list(ref = ref, mat = mat_display))
}


### Setting up for graphing ------------------------------------------------------------------------------

# A. Load and clean both
df_48h <- load_and_clean("Condition_at_48hr_Males_Fem_HM.csv")
df_14d <- load_and_clean("Condition_at_14d_Males_Fem_HM.csv")

# B. Find Top 100 for each to find the overlap
# (We use a temporary run to get the names)
get_names <- function(df, top_n) {
  df %>% rownames_to_column("name") %>%
    mutate(avg = rowMeans(select(., -name))) %>%
    slice_max(order_by = avg, n = top_n) %>% pull(name)
}

top_48h <- get_names(df_48h, 200)
top_14d <- get_names(df_14d, 200)
common_regions <- intersect(top_48h, top_14d)

# C. Process both using the SAME common regions
atlas_file <- "Allen Brain Brain Regions (names and abbriviations).xlsx"
res_48h <- process_with_atlas(df_48h, atlas_file, 5:8, common_regions)
res_14d <- process_with_atlas(df_14d, atlas_file, 5:8, common_regions)

settingupgraph <- function(heatmapvalues, referencenames) {
  ## Adding the scales for each value above for z-scores
  col_zscores <- colorRamp2(quantile(heatmapvalues, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), c("#4575B5","#BAD9E9","#FBF7c2","#FBB272","#D73127"))
  
  ## Identifying parent regions for each depth = 6 brain region (map abbrev → full name → parent name)
  
  row_info <- heatmapvalues %>%
    as.data.frame() %>%
    tibble::rownames_to_column("abbreviation") %>%
    left_join(referencenames %>% select(name, abbreviation, parent_name), by = "abbreviation")
  
  ## Visually group each parent group
  
  # Modify parent_vec to create labels only for the first row in each parent group
  parent_vec <- factor(row_info$parent_name)
  
  # Get unique parent names in the same order as they appear
  parent_groups <- unique(parent_vec)
  
  # Assign colors to each parent group
  num_parents <- length(parent_groups)
  palette_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(num_parents)
  
  # Name the colors by parent group
  parent_colors <- setNames(palette_colors, parent_groups)
  
  # Create a vector for annotations that leaves empty labels for repeated parent regions
  parent_labels <- sapply(1:length(parent_vec), function(i) {
    if (i == 1 || parent_vec[i] != parent_vec[i-1]) {
      return(as.character(parent_vec[i]))
    } else {
      return("")  # Empty label for repeated parent regions
    }
  })
  return(list(
    parent_vec = parent_vec,
    parent_colors = parent_colors,
    parent_labels = parent_labels,
    col_zscores = col_zscores
  ))
}


# Generate the master parameters based on the 48h matrix (same regions as 14d)
graphing_master <- settingupgraph(res_48h$mat, res_48h$ref)

# Ensure we have our unified color scale
all_vals <- c(res_48h$mat, res_14d$mat)
unified_col <- colorRamp2(
  quantile(all_vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 
  c("#4575B5","#BAD9E9","white","#FBB272","#D73127")
)


# --- A. Create the Left Side (Annotations) ---
# This uses your exact width/styling logic
left_ha <- rowAnnotation(
  ParentRegion = graphing_master$parent_vec,
  col = list(ParentRegion = graphing_master$parent_colors),
  show_legend = TRUE,
  width = unit(0.5, "cm")
)

left_labels <- rowAnnotation(
  ParentLabel = anno_text(
    graphing_master$parent_labels,
    gp = gpar(fontsize = 8),
    just = "left",
    location = 0,
    width = unit(2.5, "cm")
  ),
  show_annotation_name = FALSE
)

# --- B. Create the Heatmaps ---
h_48h <- Heatmap(res_48h$mat, 
                 cluster_rows = F, cluster_columns = F, 
                 name = "48h Z-Score",
                 column_title = "48 Hours",
                 col = unified_col,
                 show_row_names = FALSE) # Hide middle names

h_14d <- Heatmap(res_14d$mat, 
                 cluster_rows = F, cluster_columns = F, 
                 name = "14d Z-Score",
                 column_title = "14 Days",
                 col = unified_col,
                 show_row_names = TRUE) # Show names on the far right

# --- C. Combine Everything ---
final_h <- left_ha + left_labels + h_48h + h_14d

# --- D. Draw ---
draw(final_h, merge_legend = TRUE, column_title = "Stress Signature Overlap: 48h vs 14d")

