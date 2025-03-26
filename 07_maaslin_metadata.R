library(rio)
library(corrplot)
library(Maaslin2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
######################## ABUNDANCE TABLE ######################

# Aggregate at the Genus level
phy_genus <- tax_glom(phy_final[["V5V7"]], taxrank = "Genus")

# Extract the abundance table from the phyloseq object
abundance_table <- as.data.frame(otu_table(phy_genus))
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

corrplot(cor_matrix, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         addCoef.col = "black")



#########
####### RUN MAASLIN ##########
# 1. Prepare metadata: Set sample names as row names
rownames(merged_metadata) <- merged_metadata$Sample_orig
merged_metadata$Sample_orig <- NULL  # remove duplicate column


# 3. Ensure that Pb, Cu, and Sb are numeric 
merged_metadata$Pb <- as.numeric(merged_metadata$Pb)
merged_metadata$Cu <- as.numeric(merged_metadata$Cu)
merged_metadata$Sb <- as.numeric(merged_metadata$Sb)


######################## ABUNDANCE TABLE (ADJUSTED FOR MAASLIN2) ######################

# Aggregate at the Genus level
#phy_genus <- tax_glom(phy_regions_filtered_subset[["V2V3"]], taxrank = "Genus")
phy_genus <- tax_glom(phy_final[["V5V7"]], taxrank = "Genus")

# Extract and transpose abundance table for MaAsLin2 (samples as rows, features as columns)
abundance_table <- as.data.frame(otu_table(phy_genus))
if (taxa_are_rows(phy_genus)) {
  abundance_table <- as.data.frame(t(abundance_table))
}

abundance_table <- abundance_table[rownames(abundance_table) != "7_3_V5V7", ]
# Extract taxonomy information
tax_info <- as.data.frame(tax_table(phy_genus))

# Use Genus names for features
colnames(abundance_table) <- tax_info$Genus


######################## MAASLIN2 ANALYSIS ######################

merged_metadata <- merged_metadata[rownames(merged_metadata) != "7_3_V5V7", ]

# Ensure metadata row names match abundance table row names
rownames(merged_metadata) <- rownames(abundance_table)

# Run MaAsLin2 with combined continuous + categorical variables
fit_out <- Maaslin2(
  input_data = abundance_table,
  input_metadata = merged_metadata,
  output = "maaslin2_output_PbCu_orig_V5V7",
  fixed_effects = c("Pb", "Cu"), 
  normalization = "TSS",
  transform = "LOG",
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)


########### MAASLIN2 ON THE MOST ABUNDUNT TAXA ####################
# Convert raw counts to relative abundance (proportions)
abundance_rel <- apply(abundance_table, 2, function(x) x / sum(x))
n_top <- 50  # Adjust 
mean_abundance_sorted <- sort(colMeans(abundance_rel), decreasing = TRUE)
top_taxa <- names(mean_abundance_sorted[1:n_top])
abundance_filtered <- abundance_table[, top_taxa, drop = FALSE]
rownames(merged_metadata) <- rownames(abundance_filtered)

# Ensure metadata row names match abundance table row names
rownames(merged_metadata) <- rownames(abundance_filtered)

fit_out <- Maaslin2(
  input_data = abundance_filtered,  # Use filtered abundance
  input_metadata = merged_metadata,
  output = "maaslin2_output_top_taxa_V5V7",
  fixed_effects = c("Pb", "Cu"),
  normalization = "TSS",
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




################# Additional plots for V regions comparison ####

# Define output structure
results_dir <- "combined_analysis"
dir.create(results_dir, showWarnings = FALSE)

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

# Consensus Analysis -------------------------------------------------------
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

# Enhanced Visualization ----------------------------------------------------
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

# Statistical Analysis of Overlaps ------------------------------------------
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

#  Heatmap Visualization --------------------------------------------
heatmap_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(50)
pheatmap(
  combined_matrix,
  na_col = "grey90",
  color = rev(heatmap_colors),
  annotation_col = data.frame(
    Category = factor(rep(c("XRF", "Metals"), 
                          c(ncol(xrf_matrix), ncol(metal_matrix)))),
    row.names = colnames(combined_matrix)
  ),
  main = "Cross-region associations (≥5 regions)",
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


