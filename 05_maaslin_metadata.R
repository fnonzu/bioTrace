install.packages("rio")


# INPUT
# 1. Abundance Table (Taxa Table):
#Format: A tab-delimited (TSV) file.
# Rows: Each row represents a microbial feature (e.g., OTU, ASV, species, genus).
# Columns: Each column represents a sample
# 
# 2. Metadata Table:
#   Format: Also a tab-delimited (TSV) file.
# Rows: Each row corresponds to a sample.
# Columns: Each column represents a metadata variable related to the samples.
# Variables: These can include categorical variables (e.g., diagnosis, dysbios

######################## ABUNDANCE TABLE ######################

# Aggregate at the Genus level
phy_genus <- tax_glom(phy_regions_filtered_subset[["V5V7"]], taxrank = "Genus")

# Extract the abundance table from the phyloseq object
abundance_table <- as.data.frame(otu_table(phy_genus))

# Check orientation: MaAsLin 3 expects features as rows and samples as columns.
if (!taxa_are_rows(phy_genus)) {
  abundance_table <- t(abundance_table)
}

# Extract taxonomy information from the phyloseq object
tax_info <- as.data.frame(tax_table(phy_genus))

# Replace the row names (currently sequence/OTU names) with the Genus names
rownames(abundance_table) <- tax_info$Genus


######################## METDATA TABLE ######################

library(rio) 
library(dplyr)
# reading data from all sheets 
trace_metals_data <- import_list("~/p952/20240918_Habkern_compiled_soil_data_v2.xlsx") 

trace_metals <- trace_metals_data[["Trace metals"]]

colnames(trace_metals) <- trace_metals[1, ]
trace_metals <- trace_metals[-2, ]
trace_metals <- trace_metals[8:25, ]
names(trace_metals)[2] <- "Comments"
names(trace_metals)[3] <- "Soil weight"

trace_metals <- trace_metals %>%
  rename(
    Pb = '208  Pb  [ He ]',
    Cu = '63  Cu  [ He ]',
    Sb_121 = '121  Sb  [ He ]',
    Sb_123 = '123  Sb  [ He ]'
  )
trace_metals <- as.data.frame(trace_metals)

#average Sb
trace_metals <- trace_metals %>%
  mutate(
    Sb_121 = as.numeric(as.character(Sb_121)),
    Sb_123 = as.numeric(as.character(Sb_123))
  ) %>%
  dplyr::mutate(Sb = (Sb_121 + Sb_123)/2) %>%
  dplyr::select(-Sb_121, -Sb_123)


# Create a small data frame with row names as the sample ID
trace_metals_metadata <- trace_metals[, c("Sample", "Pb", "Cu", "Sb")]
rownames(trace_metals_metadata) <- trace_metals_metadata$Sample
#trace_metals_metadata$Sample <- NULL
rownames(trace_metals_metadata) <- NULL

############## NORMILIZE TRACE METALS - USED

## transform metadata for z-scores to fit it in and also make it numeric
trace_metals_metadata_scaled <- trace_metals_metadata %>%
  mutate(across(c(Pb, Cu, Sb), as.numeric)) %>% 
  mutate(across(c(Pb, Cu, Sb), ~ as.numeric(scale(.x))))

##############################

# Convert sample_data to a data frame again if needed
metadata_df <- as.data.frame(sample_data(phy_regions_filtered_subset[["V5V7"]]))

# 0–300 = Low
# 301–800 = Medium
# 801+ = High
metadata_df$XRF_cat <- cut(as.numeric(metadata_df$XRF),
                           breaks = c(0, 300, 1000, Inf),
                           labels = c("Low", "Medium", "High"),
                           right = FALSE)


# MERGE PHYLO AND METAL

trace_metals_metadata$Sample_orig <- gsub("\\.", "_", trace_metals_metadata$Sample)
metadata_plain <- as.data.frame(unclass(metadata_df))
trace_metals_plain <- as.data.frame(unclass(trace_metals_metadata))

rownames(metadata_df) <- metadata_df$Sample


merged_metadata <- merge(
  x = metadata_plain,
  y = trace_metals_plain,
  by = "Sample_orig",
  all.x = TRUE
)
rownames(merged_metadata) <- merged_metadata$Sample.y


# STANDIRTIZE MERGE PHYLO AND METAL

trace_metals_metadata_scaled$Sample_orig <- gsub("\\.", "_", trace_metals_metadata_scaled$Sample)
metadata_plain <- as.data.frame(unclass(metadata_df))
trace_metals_plain <- as.data.frame(unclass(trace_metals_metadata_scaled))

rownames(metadata_df) <- metadata_df$Sample


merged_metadata_scaled <- merge(
  x = metadata_plain,
  y = trace_metals_plain,
  by = "Sample_orig",
  all.x = TRUE
)
rownames(merged_metadata_scaled) <- merged_metadata_scaled$Sample.y

#### MEtaDATA ANALYSIS
merged_metadata$Pb <- as.numeric(as.character(merged_metadata$Pb))
merged_metadata$Cu <- as.numeric(as.character(merged_metadata$Cu))
merged_metadata$Sb <- as.numeric(as.character(merged_metadata$Sb))
merged_metadata$XRF <- as.numeric(as.character(merged_metadata$XRF))
# Repeat for Sb and XRF
# Compute Pearson correlation matrix (handling missing values)
cor_matrix <- cor(merged_metadata[, c("Pb", "Cu", "Sb", "XRF")], 
                  use = "complete.obs", 
                  method = "pearson")

# For Spearman (non-parametric):
cor_matrix_spearman <- cor(merged_metadata[, c("Pb", "Cu", "Sb", "XRF")], 
                           use = "complete.obs", 
                           method = "spearman")

#### MEtaDATA ANALYSIS SCALED NO XRF
merged_metadata_scaled$Pb <- as.numeric(as.character(merged_metadata_scaled$Pb))
merged_metadata_scaled$Cu <- as.numeric(as.character(merged_metadata_scaled$Cu))
merged_metadata_scaled$Sb <- as.numeric(as.character(merged_metadata_scaled$Sb))
merged_metadata_scaled$XRF <- as.numeric(as.character(merged_metadata_scaled$XRF))
# Repeat for Sb and XRF
# Compute Pearson correlation matrix (handling missing values)
cor_matrix <- cor(merged_metadata_scaled[, c("Pb", "Cu", "Sb")], 
                  use = "complete.obs", 
                  method = "pearson")

# For Spearman (non-parametric):
cor_matrix_spearman <- cor(merged_metadata_scaled[, c("Pb", "Cu", "Sb")], 
                           use = "complete.obs", 
                           method = "spearman")
# Print the matrix
print(cor_matrix)
# Install packages if needed
# install.packages("corrplot")
# install.packages("GGally")

library(corrplot)
corrplot(cor_matrix, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         addCoef.col = "black") # Adds correlation coefficients



#########
####### RUN MAASLIN ##########
# 1. Prepare metadata: Set sample names as row names
rownames(merged_metadata) <- merged_metadata$Sample_orig
merged_metadata$Sample_orig <- NULL  # remove duplicate column

# 2. Ensure XRF_cat is a factor with the proper order (adjust levels if needed)
merged_metadata$XRF_cat <- factor(merged_metadata$XRF_cat, levels = c("Low", "Medium", "High"))

# 3. Ensure that Pb, Cu, and Sb are numeric 
merged_metadata$Pb <- as.numeric(merged_metadata$Pb)
merged_metadata$Cu <- as.numeric(merged_metadata$Cu)
merged_metadata$Sb <- as.numeric(merged_metadata$Sb)




#### RUN MaAsLin2 with a new data

######################## ABUNDANCE TABLE (ADJUSTED FOR MAASLIN2) ######################

# Aggregate at the Genus level (ACTUALLY ON SPECIES?)
#phy_genus <- tax_glom(phy_regions_filtered_subset[["V2V3"]], taxrank = "Genus")
phy_genus <- tax_glom(phy_final[["V2V3"]], taxrank = "Genus")

# Extract and transpose abundance table for MaAsLin2 (samples as rows, features as columns)
abundance_table <- as.data.frame(otu_table(phy_genus))
if (taxa_are_rows(phy_genus)) {
  abundance_table <- as.data.frame(t(abundance_table))
}

# Extract taxonomy information
tax_info <- as.data.frame(tax_table(phy_genus))

# Use Genus names for features
colnames(abundance_table) <- tax_info$Genus

######################## METADATA TABLE (UNCHANGED) ######################

# [Keep your existing metadata processing code here]
# ... (content from your metadata preparation section) ...

######################## MAASLIN2 ANALYSIS ######################


# Ensure metadata row names match abundance table row names (critical!)
rownames(merged_metadata) <- rownames(abundance_table)

# Run MaAsLin2 with combined continuous + categorical variables
fit_out <- Maaslin2(
  input_data = abundance_table,
  input_metadata = merged_metadata,
  output = "maaslin2_output_PbCu_orig_V5V7",
  fixed_effects = c("Pb", "Cu"), 
  normalization = "TSS",
  transform = "LOG",
  #reference = "XRF_cat, Low",
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)


########### MAASLIN2 ON THE MOST ABUNDUNT TAXA ####################
# Convert raw counts to relative abundance (proportions)
abundance_rel <- apply(abundance_table, 2, function(x) x / sum(x))
n_top <- 50  # Adjust based on your data
mean_abundance_sorted <- sort(colMeans(abundance_rel), decreasing = TRUE)
top_taxa <- names(mean_abundance_sorted[1:n_top])
abundance_filtered <- abundance_table[, top_taxa, drop = FALSE]
rownames(merged_metadata) <- rownames(abundance_filtered)

fit_out <- Maaslin2(
  input_data = abundance_filtered,  # Use filtered abundance
  input_metadata = merged_metadata,
  output = "maaslin2_output_top_taxa_V2V3",
  fixed_effects = c("Pb", "Cu"),
  normalization = "TSS",  # Optional if already normalized
  transform = "LOG",
  plot_heatmap = TRUE
)



##################### MAASLIN2 LOOP WITH TOP TAXA FILTERING ###############
n_top <- 50  # Set the number of top taxa to keep (adjust as needed)

for (region in v_regions) {
  if (!region %in% names(phy_final)) {
    warning(paste("Region", region, "not found in phy_final. Skipping."))
    next
  }
  
  tryCatch({
    # Extract phyloseq object for current region
    phy_region <- phy_final[[region]]
    
    # Agglomerate taxa at Genus level
    phy_genus <- tax_glom(phy_region, taxrank = "Genus")
    
    # Extract and transpose abundance table
    abundance_table <- as.data.frame(otu_table(phy_genus))
    if (taxa_are_rows(phy_genus)) {
      abundance_table <- as.data.frame(t(abundance_table))
    }
    
    # Extract taxonomy and assign column names
    tax_info <- as.data.frame(tax_table(phy_genus))
    colnames(abundance_table) <- tax_info$Genus
    
    # Top taxa selection 
    # Calculate relative abundance
    abundance_rel <- t(apply(abundance_table, 1, function(x) x/sum(x)))
    
    # Select top taxa based on mean relative abundance
    mean_abundance <- colMeans(abundance_rel, na.rm = TRUE)
    mean_abundance_sorted <- sort(mean_abundance, decreasing = TRUE)
    top_taxa <- names(mean_abundance_sorted[1:n_top])
    
    # Filter abundance table to top taxa
    abundance_filtered <- abundance_table[, top_taxa, drop = FALSE]
    
    # Get sample identifiers
    phy_samples_orig <- phy_genus@sam_data$Sample_orig
    phy_sample_ids <- rownames(abundance_filtered)
    
    # Create ordering frame for metadata alignment
    phy_samples_df <- data.frame(
      Sample_orig = phy_samples_orig,
      index = seq_along(phy_samples_orig)
    ) 
    
    # Merge with metadata while preserving phy sample order
    metadata_subset <- merge(phy_samples_df, merged_metadata, 
                             by = "Sample_orig", all.x = TRUE)
    
    # Check for missing metadata
    if (any(is.na(metadata_subset$Pb))) {
      missing_samples <- metadata_subset$Sample_orig[is.na(metadata_subset$Pb)]
      warning(paste("Region", region, "has samples with missing metadata:", 
                    paste(missing_samples, collapse = ", "), ". Skipping."))
      next
    }
    
    # Restore original order and clean up
    metadata_subset <- metadata_subset[order(metadata_subset$index), ]
    metadata_subset$index <- NULL
    
    # Set rownames to match phyloseq sample IDs
    rownames(metadata_subset) <- phy_sample_ids
    
    # Final dimension check (updated to use filtered abundance)
    if (!all(rownames(abundance_filtered) == rownames(metadata_subset))) {
      stop(paste("Sample mismatch in region", region, 
                 "- Filtered abundance and metadata tables don't align"))
    }
    
    # Create unique output directory
    output_dir <- paste0("maaslin2_output_abundant_", region)
    
    # Run MaAsLin2 analysis with filtered abundance
    fit_out <- Maaslin2(
      input_data = abundance_filtered,  # Use filtered data
      input_metadata = metadata_subset,
      output = output_dir,
      fixed_effects = c("Pb", "Cu"),
      normalization = "TSS",
      transform = "LOG",
      #reference = "XRF_cat, Low",
      plot_heatmap = TRUE,
      plot_scatter = TRUE
    )
    
  }, error = function(e) {
    warning(paste("Error processing region", region, ":", e$message))
  })
}

#################### MAASLIN2 LOOP ###############


for (region in v_regions) {
  if (!region %in% names(phy_final)) {
    warning(paste("Region", region, "not found in phy_final. Skipping."))
    next
  }
  
  tryCatch({
    # Extract phyloseq object for current region
    phy_region <- phy_final[[region]]
    
    # Agglomerate taxa at Genus level
    phy_genus <- tax_glom(phy_region, taxrank = "Genus")
    
    # Extract and transpose abundance table
    abundance_table <- as.data.frame(otu_table(phy_genus))
    if (taxa_are_rows(phy_genus)) {
      abundance_table <- as.data.frame(t(abundance_table))
    }
    
    # Extract taxonomy and assign column names
    tax_info <- as.data.frame(tax_table(phy_genus))
    colnames(abundance_table) <- tax_info$Genus
    
    # Get sample identifiers
    phy_samples_orig <- phy_genus@sam_data$Sample_orig
    phy_sample_ids <- rownames(abundance_table)
    
    # Create ordering frame for metadata alignment
    phy_samples_df <- data.frame(
      Sample_orig = phy_samples_orig,
      index = seq_along(phy_samples_orig)
    )
    
    # Merge with metadata while preserving phy sample order
    metadata_subset <- merge(phy_samples_df, merged_metadata, 
                             by = "Sample_orig", all.x = TRUE)
    
    # Check for missing metadata
    if (any(is.na(metadata_subset$XRF_cat))) {
      missing_samples <- metadata_subset$Sample_orig[is.na(metadata_subset$XRF_cat)]
      warning(paste("Region", region, "has samples with missing metadata:", 
                    paste(missing_samples, collapse = ", "), ". Skipping."))
      next
    }
    
    # Restore original order and clean up
    metadata_subset <- metadata_subset[order(metadata_subset$index), ]
    metadata_subset$index <- NULL
    
    # Set rownames to match phyloseq sample IDs
    rownames(metadata_subset) <- phy_sample_ids
    
    # Final dimension check
    if (!all(rownames(abundance_table) == rownames(metadata_subset))) {
      stop(paste("Sample mismatch in region", region, 
                 "- Abundance and metadata tables don't align"))
    }
    
    # Create unique output directory
    output_dir <- paste0("maaslin2_output_", region)
    
    # Run MaAsLin2 analysis
    fit_out <- Maaslin2(
      input_data = abundance_table,
      input_metadata = metadata_subset,
      output = output_dir,
      fixed_effects = c("XRF_cat", "Cu"),
      normalization = "TSS",
      transform = "LOG",
      reference = "XRF_cat, Low",
      plot_heatmap = TRUE,
      plot_scatter = TRUE
    )
    
  }, error = function(e) {
    warning(paste("Error processing region", region, ":", e$message))
  })
}



####### SIMPLE ANALYSIS

region_paths <- c(
  V1V2 = "maaslin2_output_V1V2/significant_results.tsv",
  V2V3 = "maaslin2_output_V1V2/significant_results.tsv",
  V3V4 = "maaslin2_output_V3V4/significant_results.tsv",
  V4V5 = "maaslin2_output_V1V2/significant_results.tsv",
  V5V7 = "maaslin2_output_V5V7/significant_results.tsv",
  V7V9 = "maaslin2_output_V1V2/significant_results.tsv"
)
library(readr)
region_dfs <- lapply(region_paths, read_tsv)

# Step 1: Extract "feature" column from each df and find the intersection
common_features <- Reduce(intersect, lapply(region_dfs, function(df) df$feature))

# Step 2: Filter each sub-data frame to retain only common "feature" values
region_dfs_filtered <- lapply(region_dfs, function(df) df %>% filter(feature %in% common_features))

# Process each region to get the minimal q value per feature
processed_regions <- map(names(region_dfs_filtered), ~{
  region_data <- region_dfs_filtered[[.x]]
  region_data %>%
    group_by(feature) %>%
    summarise(min_qval = min(qval), .groups = 'drop') %>%
    mutate(region = .x)
}) %>% bind_rows()

# Create a wide dataframe with regions as columns
qval_wide <- processed_regions %>%
  pivot_wider(names_from = region, values_from = min_qval)

# Identify features present in at least two regions
overlap_features <- qval_wide %>%
  filter(rowSums(!is.na(select(., -feature))) >= 2) %>%
  pull(feature)


######## WORK #####
# Function to process data for specific metadata values
process_metadata <- function(metadata_values, title_suffix) {
  # Process each region for specified metadata values
  processed <- map(names(region_dfs_filtered), ~{
    region_data <- region_dfs_filtered[[.x]] %>%
      filter(value %in% metadata_values) %>%
      group_by(feature) %>%
      summarise(min_qval = min(qval), .groups = 'drop') %>%
      mutate(region = .x)
  }) %>% bind_rows()
  
  # Create wide format and find overlapping features
  qval_wide <- processed %>%
    pivot_wider(names_from = region, values_from = min_qval)
  
  overlap_features <- qval_wide %>%
    filter(rowSums(!is.na(select(., -feature))) >= 2) %>%
    pull(feature)
  
  # Prepare plot data
  plot_data <- qval_wide %>%
    filter(feature %in% overlap_features) %>%
    pivot_longer(-feature, names_to = "region", values_to = "min_qval") %>%
    filter(!is.na(min_qval))
  
  # Generate heatmap
  ggplot(plot_data, aes(x = region, y = feature, fill = -log10(min_qval))) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "white", 
      high = "red",
      name = "-log10(q)",
      limits = c(0, NA)  # Ensure same scale for both plots
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste("Significant Features for", title_suffix),
      x = "Region", 
      y = "Feature"
    )
}

# Create plots for each contamination type
xrf_plot <- process_metadata(
  metadata_values = c("High", "Medium"),  # Add other XRF categories if present
  title_suffix = "XRF Contamination"
)

cu_plot <- process_metadata(
  metadata_values = "Cu",
  title_suffix = "Cu Contamination"
)

# Display plots side by side
library(patchwork)
xrf_plot + cu_plot + plot_layout(ncol = 2)




# Modified processing function for coefficients
process_metadata_coef <- function(metadata_values, title_suffix) {
  processed <- map(names(region_dfs_filtered), ~{
    region_data <- region_dfs_filtered[[.x]] %>%
      filter(value %in% metadata_values) %>%
      group_by(feature) %>%
      slice_min(qval, n = 1) %>%  # Get most significant association per feature
      summarise(coef = first(coef), .groups = 'drop') %>%
      mutate(region = .x)
  }) %>% bind_rows()
  
  qval_wide <- processed %>%
    pivot_wider(names_from = region, values_from = coef)
  
  overlap_features <- qval_wide %>%
    filter(rowSums(!is.na(select(., -feature))) >= 2) %>%
    pull(feature)
  
  plot_data <- qval_wide %>%
    filter(feature %in% overlap_features) %>%
    pivot_longer(-feature, names_to = "region", values_to = "coef") %>%
    filter(!is.na(coef))
  
  ggplot(plot_data, aes(x = region, y = feature, fill = coef)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "Coefficient",
      limits = c(-max(abs(plot_data$coef)), max(abs(plot_data$coef)))  # Symmetric scale
    ) +
    theme_minimal(14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "italic"),  # Emphasize taxonomic names
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(
      title = paste("Feature Associations for", title_suffix),
      x = NULL, 
      y = NULL
    ) +
    geom_text(aes(label = round(coef, 2)), size = 3.5)  # Add coefficient values
}

# Generate plots
xrf_coef_plot <- process_metadata_coef(
  metadata_values = c("High", "Medium"),
  title_suffix = "XRF Contamination"
)

cu_coef_plot <- process_metadata_coef(
  metadata_values = "Cu",
  title_suffix = "Cu Contamination"
)

# Arrange plots
(xrf_coef_plot | cu_coef_plot) + plot_annotation(title = "Microbial Associations with Contamination Types")

# 
# ################## merging results from different regions:###############
# region_paths <- c(
#   V1V2 = "maaslin2_outputV1V2/significant_results.tsv",
#   V3V4 = "maaslin2_outputV3V4/significant_results.tsv",
#   V5V7 = "maaslin2_outputV5V7/significant_results.tsv"
# )
# 
# all_results <- lapply(region_paths, function(path) {
#   df <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
#   df$Region <- gsub("maaslin2_|_output.*", "", basename(path))  # Add region label
#   return(df)
# })
# 
# # Combine into one data frame
# combined_results <- do.call(rbind, all_results)
# combined_results$Region <- substr(rownames(combined_results), start = 1, stop = 4)
# 
# overlapping_taxa <- combined_results %>%
#   group_by(feature) %>%
#   filter(n_distinct(Region) >= 2) %>%  # Taxa significant in ≥2 regions
#   arrange(feature, Region)
# 
# # Venn diagram (for 2-3 regions)
# library(VennDiagram)
# venn.diagram(
#   x = list(
#     V1V2 = unique(filter(combined_results, Region == "V1V2")$feature),
#     V3V4 = unique(filter(combined_results, Region == "V3V4")$feature),
#     V5V7 = unique(filter(combined_results, Region == "V5V7")$feature)
#   ),
#   filename = "taxa_overlap.png",
#   fill = c("blue", "green", "red")
# )
# 
# 
# ##### HEATMAP COMBINED #######
# library(tidyverse)
# library(pheatmap)
# 
# # 1. Read and combine results from all regions
# region_paths <- c(
#   V1V2 = "maaslin2_output_V1V2/significant_results.tsv",
#   V3V4 = "maaslin2_output_V3V4/significant_results.tsv",
#   V5V7 = "maaslin2_output_V5V7/significant_results.tsv"
# )
# 
# combined_results <- map_dfr(region_paths, ~{
#   read.delim(.x) %>%
#     mutate(Region = str_extract(.x, "V[0-9]V[0-9]"))
# }, .id = "Region")
# 
# # 2. Filter for overlapping taxa (significant in ≥2 regions)
# overlap_taxa <- combined_results %>%
#   group_by(feature, metadata) %>%
#   filter(n_distinct(Region) >= 2) %>%
#   ungroup()
# 
# # 3. Split data into XRF categories and metals
# xrf_data <- overlap_taxa %>%
#   filter(metadata == "XRF_cat") %>%
#   mutate(metadata = paste(metadata, value, sep = "_"))
# 
# metal_data <- overlap_taxa %>%
#   filter(metadata %in% c("Sb", "Cu", "Pb"))
# 
# # 4. Create heatmap matrices
# create_matrix <- function(data) {
#   data %>%
#     group_by(feature, metadata) %>%
#     summarise(AvgCoef = mean(coef)) %>%  # Average coefficients across regions
#     pivot_wider(names_from = metadata, values_from = AvgCoef) %>%
#     column_to_rownames("feature") %>%
#     as.matrix()
# }
# 
# xrf_matrix <- create_matrix(xrf_data)
# metal_matrix <- create_matrix(metal_data)
# 
# # 5. Combined heatmap with annotation
# # Modified matrix creation with alignment
# create_aligned_matrix <- function(data, all_taxa) {
#   data %>%
#     group_by(feature, metadata) %>%
#     summarise(AvgCoef = mean(coef), .groups = 'drop') %>%
#     complete(feature = all_taxa, metadata, fill = list(AvgCoef = NA)) %>%
#     pivot_wider(names_from = metadata, values_from = AvgCoef) %>%
#     column_to_rownames("feature") %>%
#     as.matrix()
# }
# 
# # Get unique taxa from ALL significant features
# all_taxa <- unique(overlap_taxa$feature)
# 
# # Create aligned matrices
# xrf_matrix <- create_aligned_matrix(xrf_data, all_taxa)
# metal_matrix <- create_aligned_matrix(metal_data, all_taxa)
# 
# # Now combine them
# combined_matrix <- cbind(xrf_matrix, metal_matrix)
# 
# # Remove taxa with all NA values (optional)
# combined_matrix <- combined_matrix[rowSums(!is.na(combined_matrix)) > 0, ]
# annotation_col <- data.frame(
#   Category = c(rep("XRF", ncol(xrf_matrix)), 
#                rep("Metals", ncol(metal_matrix)))
# )
# rownames(annotation_col) <- colnames(combined_matrix)
# 
# # Impute NA values for clustering only
# clust_matrix <- combined_matrix
# clust_matrix[is.na(clust_matrix)] <- 0  # Replace NAs with 0 for clustering
# 
# # Compute the hierarchical clustering (hclust object)
# row_hclust <- hclust(dist(clust_matrix), method = "complete")
# 
# # Now generate the heatmap using the original matrix (with NA values)
# pheatmap(
#   combined_matrix,
#   na_col = "gray90",  # Color for NA values
#   color = colorRampPalette(c("blue", "white", "red"))(50),
#   annotation_col = annotation_col,
#   main = "Cross-Region Associations",
#   cluster_rows = row_hclust,  # Pass the hclust object directly
#   cluster_cols = FALSE,
#   filename = "combined_heatmap.png"
# )
# 


################# UPSET####

# Define output structure
results_dir <- "combined_analysis"
dir.create(results_dir, showWarnings = FALSE)

# 1. Enhanced Results Merging with Error Handling -------------------------------
region_paths <- c(
  V1V2 = "maaslin2_output_V1V2/significant_results.tsv",
  V2V3 = "maaslin2_output_V1V2/significant_results.tsv",
  V3V4 = "maaslin2_output_V3V4/significant_results.tsv",
  V4V5 = "maaslin2_output_V1V2/significant_results.tsv",
  V5V7 = "maaslin2_output_V5V7/significant_results.tsv",
  V7V9 = "maaslin2_output_V1V2/significant_results.tsv"
)

# Safe reading with error handling
all_results <- lapply(names(region_paths), function(region) {
  path <- region_paths[[region]]
  tryCatch({
    df <- read.delim(path, sep = "\t", check.names = FALSE) %>% 
      mutate(Region = region,
             Association = paste0(metadata, ":", value),
             Direction = ifelse(coef > 0, "Positive", "Negative"))
    return(df)
  }, error = function(e) {
    warning(paste("Failed to read", path, ":", e$message))
    return(NULL)
  })
}) %>% bind_rows()

# 2. Consensus Analysis -------------------------------------------------------
# Calculate consensus metrics
consensus_results <- all_results %>%
  group_by(feature, Association) %>%
  summarise(
    Region_Count = n_distinct(Region),
    Regions = paste(unique(Region), collapse = ";"),
    Avg_Coef = mean(coef),
    Coef_Consistency = sd(coef)/abs(mean(coef)),  # Lower = more consistent
    Min_qval = min(qval),
    .groups = "drop"
  ) %>%
  arrange(-Region_Count, Coef_Consistency, Min_qval)

# Save consensus table
write_csv(consensus_results, file.path(results_dir, "consensus_associations.csv"))

# 3. Enhanced Visualization ----------------------------------------------------
library(UpSetR)

upset_data <- all_results %>%
  distinct(feature, Region) %>%
  mutate(Value = 1) %>%
  pivot_wider(names_from = Region, values_from = Value, values_fill = 0) %>%
  as.data.frame()

pdf(file.path(results_dir, "taxa_overlap_upset.pdf"), width = 10, height = 6)
upset(upset_data, 
      nsets = length(region_paths),
      mainbar.y.label = "Significant Taxa",
      sets.x.label = "Taxa per Region")
dev.off()

# Taxa unique to V5V7
consensus_results %>% 
  filter(str_detect(Regions, "V5V7") & !str_detect(Regions, "[;]"))

taxa_presence <- all_results %>%
  group_by(feature) %>%
  summarise(
    Total_Regions = n_distinct(Region),
    Regions = paste(unique(Region), collapse = ";")
  ) %>%
  filter(Total_Regions >= 5 & Total_Regions <= 6) %>%
  arrange(-Total_Regions, feature)

# 4. Statistical Analysis of Overlaps ------------------------------------------
# Fisher's exact test for association overlap
create_contingency <- function(taxa) {
  matrix(c(
    sum(taxa %in% all_results$feature),  # Significant in all regions
    length(unique(all_results$feature)) - sum(taxa %in% all_results$feature),
    sum(!taxa %in% all_results$feature),  # Not significant
    nrow(upset_data) - sum(!taxa %in% all_results$feature)
  ), nrow = 2)
}

# Test for each association type
association_types <- unique(all_results$Association)
enrichment_results <- lapply(association_types, function(a) {
  taxa <- unique(all_results$feature[all_results$Association == a])
  mat <- create_contingency(taxa)
  ft <- fisher.test(mat)
  data.frame(
    Association = a,
    Odds_Ratio = ft$estimate,
    p_value = ft$p.value,
    Significant_Taxa = length(taxa)
  )
}) %>% bind_rows() %>% arrange(p_value)

# Save enrichment results
write_csv(enrichment_results, file.path(results_dir, "association_enrichment.csv"))


############# Heatmap ###########

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Safe reading with error handling
combined_results <- map_dfr(names(region_paths), function(region) {
  path <- region_paths[[region]]
  tryCatch({
    read.delim(path) %>%
      mutate(Region = region)  # Use actual region name from vector
  }, error = function(e) {
    warning(paste("Failed to read", path, ":", e$message))
    return(NULL)
  })
})

# 2. Filter for overlapping taxa ----------------------------------------------
overlap_taxa <- combined_results %>%
  group_by(feature, metadata) %>%
  filter(n_distinct(Region) >= 2) %>%
  ungroup()

# 3. Enhanced Matrix Creation -------------------------------------------------
create_aligned_matrix <- function(data, all_taxa, min_regions = 2) {
  data %>%
    group_by(feature, metadata) %>%
    filter(n() >= min_regions) %>%  # Require presence in ≥2 regions
    summarise(AvgCoef = mean(coef, na.rm = TRUE), .groups = 'drop') %>%
    complete(feature = all_taxa, metadata, fill = list(AvgCoef = NA)) %>%
    pivot_wider(names_from = metadata, values_from = AvgCoef) %>%
    column_to_rownames("feature") %>%
    as.matrix()
}

# Get unique taxa from overlapping features
all_taxa <- unique(overlap_taxa$feature)
# 2. Filter for overlapping taxa ----------------------------------------------
overlap_taxa <- combined_results %>%
  group_by(feature, metadata) %>%
  filter(n_distinct(Region) >= 5) %>%  # Changed from 2 to 5
  ungroup()

# 3. Enhanced Matrix Creation -------------------------------------------------
create_aligned_matrix <- function(data, all_taxa, min_regions = 5) {  # Changed min_regions to 5
  data %>%
    group_by(feature, metadata) %>%
    filter(n() >= min_regions) %>%  # Updated comment: Require presence in ≥5 regions
    summarise(AvgCoef = mean(coef, na.rm = TRUE), .groups = 'drop') %>%
    complete(feature = all_taxa, metadata, fill = list(AvgCoef = NA)) %>%
    pivot_wider(names_from = metadata, values_from = AvgCoef) %>%
    column_to_rownames("feature") %>%
    as.matrix()
}
# 4. Create and combine matrices with error handling --------------------------
tryCatch({
  xrf_matrix <- combined_results %>%
    filter(metadata == "XRF_cat") %>%
    mutate(metadata = paste(metadata, value, sep = "_")) %>%
    create_aligned_matrix(all_taxa)
  
  metal_matrix <- combined_results %>%
    filter(metadata %in% c("Sb", "Cu", "Pb")) %>%
    create_aligned_matrix(all_taxa)
  
  combined_matrix <- cbind(xrf_matrix, metal_matrix)
}, error = function(e) {
  stop(paste("Matrix creation failed:", e$message))
})

# 5. Enhanced Clustering with NA handling -------------------------------------
safe_cluster <- function(mat) {
  # Convert NA to 0 for clustering only
  clust_mat <- mat
  clust_mat[is.na(clust_mat)] <- 0
  
  # Remove rows with zero variance
  row_vars <- apply(clust_mat, 1, var, na.rm = TRUE)
  clust_mat <- clust_mat[row_vars > 0, ]
  
  # Calculate distance with error handling
  tryCatch({
    hclust(dist(clust_mat), method = "complete")
  }, error = function(e) {
    message("Using alternative clustering due to error: ", e$message)
    hclust(dist(clust_mat, method = "manhattan"), method = "ward.D2")
  })
}

row_hclust <- safe_cluster(combined_matrix)

# 6. Robust Heatmap Visualization --------------------------------------------
heatmap_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(50)

pheatmap(
  combined_matrix,
  na_col = "grey90",
  color = rev(heatmap_colors),  # Blue=negative, Red=positive
  annotation_col = data.frame(
    Category = factor(rep(c("XRF", "Metals"), 
                          c(ncol(xrf_matrix), ncol(metal_matrix)))),
    row.names = colnames(combined_matrix)
  ),
  main = "Cross-Region Associations (≥2 regions)",
  cluster_rows = row_hclust,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 9,
  angle_col = 45,
  gaps_col = ncol(xrf_matrix),
  filename = file.path(results_dir, "combined_heatmap.png"),
  width = 12,
  height = 14
)

# 6. Robust Heatmap Visualization --------------------------------------------
pheatmap(
  combined_matrix,
  na_col = "grey90",
  color = rev(heatmap_colors),
  annotation_col = data.frame(
    Category = factor(rep(c("XRF", "Metals"), 
                          c(ncol(xrf_matrix), ncol(metal_matrix)))),
    row.names = colnames(combined_matrix)
  ),
  main = "Cross-Region Associations (≥5 regions)",  # Updated title
  cluster_rows = row_hclust,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 9,
  angle_col = 45,
  gaps_col = ncol(xrf_matrix),
  filename = file.path(results_dir, "combined_heatmap2.png"),
  width = 12,
  height = 14
)

############################ OTHER METADATA ############################
## PH ##
pH <- trace_metals_data[["pH"]]
names(pH)[3] <- "pH in CaCl2"
pH[["pH in CaCl2"]] <- as.numeric(pH[["pH in CaCl2"]])

average_pH <- mean(pH[1:3, "pH in CaCl2"], na.rm = TRUE)
pH[1:3, "pH in CaCl2"] <- average_pH
pH <- pH[-c(2, 3), ]
average_pH <- mean(pH[6:8, "pH in CaCl2"], na.rm = TRUE)
pH[6:8, "pH in CaCl2"] <- average_pH
pH <- pH[-c(7, 8), ]
average_pH <- mean(pH[12:14, "pH in CaCl2"], na.rm = TRUE)
pH[12:14, "pH in CaCl2"] <- average_pH
pH <- pH[-c(13, 14), ]
rownames(pH) <- NULL

## ORGANIC MATTER 
org_mat <- trace_metals_data[["Organic Matter %"]]

## GRAIN SIZE
grain <- trace_metals_data[["grain size"]]
colnames(grain)[1] <- "Sample Name"
grain <- grain[-1, ]
average_grain <- grain %>%
  filter(grepl("- Average$", `Sample Name`))  # Keep only averages

average_grain <- average_grain %>%
  mutate(
    `clay fraction (%)` = as.numeric(as.character(`clay fraction (%)`)),
    `silt fraction (%)` = as.numeric(as.character(`silt fraction (%)`)),
    `sand fraction (%)` = as.numeric(as.character(`sand fraction (%)`))
  )
average_grain <- average_grain %>% 
  mutate(sample_group = sub("^([0-9]+).*", "\\1", `Sample Name`))
average_grain <- average_grain %>%
  group_by(sample_group) %>%
  summarise(
    clay_fraction = mean(`clay fraction (%)`, na.rm = TRUE),
    silt_fraction = mean(`silt fraction (%)`, na.rm = TRUE),
    sand_fraction = mean(`sand fraction (%)`, na.rm = TRUE)
  )
average_grain$sample_group <- as.numeric(as.character(average_grain$sample_group))
rownames(average_grain) <- average_grain$sample_group # WEIRD

#PLOTTING
library(ggplot2)
library(tidyr)
library(gridExtra)

# Add sample numbers from row names to each data frame
pH$Sample <- as.numeric(rownames(pH))
org_mat$Sample <- as.numeric(rownames(org_mat))

# Standardize grain size data to sum to 100% per row
average_grain[, 2:4] <- average_grain[, 2:4] / rowSums(average_grain[, 2:4]) * 100

# Prepare grain data for stacked bars
grain_long <- average_grain %>% 
  mutate(Sample = as.numeric(average_grain$sample_group)) %>% 
  pivot_longer(cols = 2:4, names_to = "Grain_Size", values_to = "Percentage")

# Open a PNG device to save the plot
png("plots/metadata_plot_3inrow.png", width = 1500, height = 700)
# Create pH plot
p_pH <- ggplot(pH, aes(x = Sample, y = pH[,3])) +
  geom_line(color = "blue", linewidth = 0.8) +
  geom_point(color = "blue", size = 2) +
  labs(title = "pH in CaCl2", x = "Sample", y = "pH value") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal()

# Create organic matter plot
p_org <- ggplot(org_mat, aes(x = Sample, y = org_mat[,3])) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_point(color = "darkgreen", size = 2) +
  labs(title = "Organic matter", x = "Sample", y = "Organic matter percentage %") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal()

# Create stacked grain size plot
p_grain <- ggplot(grain_long, aes(x = factor(Sample), y = Percentage, fill = Grain_Size)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Grain composition", x = "Sample", y = "Percentage (%)") +
  scale_x_discrete(breaks = 1:12) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine all plots in a grid
grid.arrange(p_pH, p_org, p_grain, 
             ncol = 3,
             widths = c(1, 1, 1.2))  # Adjust relative widths
dev.off()
# Add sample numbers from row names to each data frame
pH$Sample <- as.numeric(rownames(pH))
org_mat$Sample <- as.numeric(rownames(org_mat))

# Standardize grain size data to sum to 100% per row
average_grain[, 2:4] <- average_grain[, 2:4] / rowSums(average_grain[, 2:4]) * 100

# Prepare grain data for stacked bars
grain_long <- average_grain %>% 
  mutate(Sample = as.numeric(average_grain$sample_group)) %>% 
  pivot_longer(cols = 2:4, names_to = "Grain_Size", values_to = "Percentage")

# Open a PNG device to save the plot
png("plots/metadata_plot.png", width = 1500, height = 900)  # Increased height

# Create pH plot
p_pH <- ggplot(pH, aes(x = Sample, y = pH[,3])) +
  geom_line(color = "blue", linewidth = 0.8) +
  geom_point(color = "blue", size = 2) +
  labs(title = "pH in CaCl₂", x = "Sample", y = "pH value") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 5, "pt"))

# Create organic matter plot
p_org <- ggplot(org_mat, aes(x = Sample, y = org_mat[,3])) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_point(color = "darkgreen", size = 2) +
  labs(title = "Organic matter", x = "Sample", y = "Organic matter (%)") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 5, "pt"))

# Create stacked grain size plot
p_grain <- ggplot(grain_long, aes(x = factor(Sample), y = Percentage, fill = Grain_Size)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Grain composition", x = "Sample", y = "Percentage (%)") +
  scale_x_discrete(breaks = 1:12) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.margin = margin(5, 15, 5, 15, "pt"))

# Combine plots in new layout
grid.arrange(
  arrangeGrob(p_pH, p_org, ncol = 1, heights = c(1, 1)),  # Left column
  p_grain,  # Right column
  ncol = 2, 
  widths = c(1, 1.5)  # Right plot is 1.5x wider than left column
)

dev.off()
