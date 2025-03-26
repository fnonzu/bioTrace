# ENVIRONMENTAL DATA ANALYSIS WITH CCA -----------------------------------------
# This script performs Canonical Correspondence Analysis (CCA) to explore 
# relationships between microbial communities and environmental variables (Pb, Cu)

# LIBRARIES --------------------------------------------------------------------
library(phyloseq)
library(vegan)
library(tidyverse)

# DATA PREPARATION -------------------------------------------------------------
# Original metadata processing
# Convert sample names to match phyloseq format (replace "." with "_")
merged_metadata <- merged_metadata %>%
  dplyr::mutate(Sample_orig = str_replace(Sample.y, "\\.", "_")) 
rownames(merged_metadata) <- NULL
merged_metadata <- column_to_rownames(merged_metadata, "Sample_orig")

# Extract environmental variables and check sample overlap
# env_data <- merged_metadata %>% dplyr::select(Pb, Cu)
# common_samples <- intersect(rownames(otu_table), rownames(env_data))

# Handle case with no overlapping samples
if(length(common_samples) == 0){
  warning("No overlapping samples for region: ", region)
  next  # Skip rest of loop iteration
}

# Scaled metadata processing (duplicate section for normalized data)
# TODO: Consider creating function for duplicate processing steps
merged_metadata_scaled <- merged_metadata_scaled %>%
  dplyr::mutate(Sample_orig = str_replace(Sample.y, "\\.", "_")) 
rownames(merged_metadata_scaled) <- NULL
merged_metadata_scaled <- column_to_rownames(merged_metadata_scaled, "Sample_orig")

env_data <- merged_metadata_scaled %>% dplyr::select(Pb, Cu)
common_samples <- intersect(rownames(otu_table), rownames(env_data))
if(length(common_samples) == 0){
  warning("No overlapping samples for region: ", region)
  next
}
# MAIN ANALYSIS PIPELINE -------------------------------------------------------
# Configure output directory
plots_dir <- path.expand("~/plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Process each target region
for (region in c('V1V2','V2V3','V4V5','V5V7')) {
  # DATA PREPROCESSING ---------------------------------------------------------
  # Extract and prepare phyloseq object
  phy_subset <- phy_final[[region]] %>%
    tax_glom(taxrank = "Genus")  # Aggregate at genus level
  
  # TAXONOMY PROCESSING --------------------------------------------------------
  # Create meaningful taxon labels
  tax_table <- as.data.frame(tax_table(phy_subset), stringsAsFactors = FALSE)
  tax_table$Genus <- ifelse(is.na(tax_table$Genus),
                            paste0("Unk_", tax_table$Family),  # Use Family if Genus missing
                            tax_table$Genus)
  taxa_names(phy_subset) <- make.unique(tax_table$Genus)
  
  # SAMPLE NAME PROCESSING -----------------------------------------------------
  sample_names(phy_subset) <- sub("(_V.*)$", "", sample_names(phy_subset))  # Remove region suffix
  
  # OTU TABLE PROCESSING -------------------------------------------------------
  otu_table_raw <- as.data.frame(otu_table(phy_subset))
  
  # Filter and clean OTU data
  otu_table_df <- otu_table_raw %>%
    rownames_to_column("Sample_orig") %>%
    filter(Sample_orig %in% rownames(env_data) & !Sample_orig %in% c("7_3")) %>%
    column_to_rownames("Sample_orig") %>%
    .[rowSums(.) > 0, ]  # Remove zero-count samples
  
  # Skip region if no valid samples
  if(nrow(otu_table_df) == 0){
    warning("No valid samples for region: ", region)
    next
  }
  
  # COMMUNITY ANALYSIS ---------------------------------------------------------
  # Hellinger transformation for compositionality
  otu_hellinger <- decostand(otu_table_df, method = "hellinger")
  
  # Environmental data subsetting
  env_subset <- env_data[rownames(otu_hellinger), , drop = FALSE]
  
  # VARIABLE SELECTION ---------------------------------------------------------
  # Initial CCA for VIF calculation
  cca_initial <- cca(otu_hellinger ~ ., env_subset)
  vif_scores <- vif.cca(cca_initial)
  
  # Filter variables with VIF < 10
  env_filtered <- env_subset[, vif_scores < 10, drop = FALSE]
  
  # Skip region if no variables remain
  if(ncol(env_filtered) == 0){
    warning("No valid environmental variables for region: ", region)
    next
  }
  
  # CCA EXECUTION --------------------------------------------------------------
  cca_result <- cca(otu_hellinger ~ ., env_filtered)
  
  # Variance calculation
  variance <- round(cca_result$CCA$eig / sum(cca_result$CCA$eig) * 100, 1)
  x_label <- paste0("CCA1 (", variance[1], "%)")
  y_label <- paste0("CCA2 (", variance[2], "%)")
  
  # STATISTICAL TESTING --------------------------------------------------------
  cca_test <- anova(cca_result, permutations = 999)
  
  # VISUALIZATION --------------------------------------------------------------
  file_path <- file.path(plots_dir, paste0("CCA_plot_PbCu_norm_", region, ".png"))
  png(file_path, width = 1500, height = 1500, res = 300)
  
  # Base plot setup
  plot(cca_result, type = "n", main = paste("CCA for", region),
       xlab = x_label, ylab = y_label, scaling = 2)
  
  # Plot elements
  mtext("normalized Pb concentration", side = 1, line = 4, adj = 1, col = "darkgreen", cex = 0.5)
  points(cca_result, display = "sites", pch = 16, col = "blue", cex = 0.8)
  text(cca_result, display = "sites", col = "blue", cex = 0.6, pos = 3)
  text(cca_result, display = "bp", col = "red", cex = 0.8, font = 2)
  
  # Taxa labeling (commented filter - consider activating)
  # taxa_dist <- sqrt(taxa_scores[,1]^2 + taxa_scores[,2]^2)
  # top_taxa <- names(taxa_dist)[taxa_dist > quantile(taxa_dist, 0.95)]
  #orditorp(cca_result, display = "species", col = "darkgray", 
  #        cex = 0.7, air = 1.2, priority = taxa_dist)
  
  # Environmental contours
  ordisurf(cca_result, env_filtered$Pb, add = TRUE, 
           col = "darkgreen", labcex = 0.8)
  
  # Legend and annotations
  legend("topleft",
         legend = c(
           paste("Explained variance (CCA1 + CCA2):",
                 round(sum(cca_result$CCA$eig[1:2]/cca_result$tot.chi*100), 1), "%"),
           paste("Model p-value:", cca_test$`Pr(>F)`[1])
         ),
         bty = "n", cex = 0.8)
  
  dev.off()
}


summary(cca_result)
anova(cca_result)