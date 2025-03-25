# Load libraries
library(phyloseq)
library(vegan)
library(tidyverse)

# Fix sample names in metadata to match phyloseq (replace "." with "_")
merged_metadata <- merged_metadata %>%
  dplyr::mutate(Sample_orig = str_replace(Sample.y, "\\.", "_")) #%>% # Converts "3.1" to "3_1"
rownames(merged_metadata) <- NULL
merged_metadata <- column_to_rownames(merged_metadata, "Sample_orig")

# Extract environmental variables (XRF, Cu, Pb, Sb) and sample IDs
#env_data <- merged_metadata %>% dplyr::select(XRF, Cu, Pb, Sb)
env_data <- merged_metadata %>% dplyr::select(XRF, Cu)
#env_data <- merged_metadata %>% dplyr::select(XRF)
common_samples <- intersect(rownames(otu_table), rownames(env_data))
if(length(common_samples) == 0){
  warning("No overlapping samples for region: ", region)
  next  # This skips the rest of the loop iteration
}



#### SCALED METADATA NO XRF METADATA IS GENERATED IN SCRIPT 05
# Fix sample names in metadata to match phyloseq (replace "." with "_")
merged_metadata_scaled <- merged_metadata_scaled %>%
  dplyr::mutate(Sample_orig = str_replace(Sample.y, "\\.", "_")) #%>% # Converts "3.1" to "3_1"
rownames(merged_metadata_scaled) <- NULL
merged_metadata_scaled <- column_to_rownames(merged_metadata_scaled, "Sample_orig")

env_data <- merged_metadata_scaled %>% dplyr::select(Pb, Cu)
common_samples <- intersect(rownames(otu_table), rownames(env_data))
if(length(common_samples) == 0){
  warning("No overlapping samples for region: ", region)
  next
}

####### WORK
## HERE change the variables for XRF of Pb to compare and check the data before to make sure it's the right metadata in env_data
# Create the plots directory if it doesn't exist
plots_dir <- path.expand("~/plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

for (region in c('V1V2','V2V3','V4V5','V5V7')) {
  # Extract the phyloseq subset for the region
  phy_subset <- phy_final[[region]]
  phy_subset <- tax_glom(phy_subset, taxrank = "Genus")
  # ------ FIX TAXON LABELS: Replace sequence IDs with genus names ------
  
  # Extract taxonomy table and replace taxa names with Genus (or Family if NA)
  tax_table <- as.data.frame(tax_table(phy_subset), stringsAsFactors = FALSE)
  tax_table$Genus <- ifelse(is.na(tax_table$Genus), 
                            paste0("Unk_", tax_table$Family),  # Fallback to Family + "Unk"
                            tax_table$Genus)
  taxa_names(phy_subset) <- make.unique(tax_table$Genus)  # Ensure unique names
  
  # Adjust sample names in the phyloseq object to remove the region suffix
  sample_names(phy_subset) <- sub("(_V.*)$", "", sample_names(phy_subset))
  
  # Extract the OTU table as a data frame
  otu_table_raw <- as.data.frame(otu_table(phy_subset))
  
  # Filter and format OTU table to include only samples that overlap with env_data
  otu_table_df <- otu_table_raw %>%
    rownames_to_column("Sample_orig") %>%
    filter(Sample_orig %in% rownames(env_data) & !Sample_orig %in% c("7_3")) %>%
    column_to_rownames("Sample_orig")
  
  
  # Remove samples with zero counts to ensure all row sums > 0
  otu_table_df <- otu_table_df[rowSums(otu_table_df) > 0, ]
  
  if(nrow(otu_table_df) == 0){
    warning("No samples with non-zero counts for region: ", region)
    next
  }
  
  # Check for overlapping samples
  if(nrow(otu_table_df) == 0){
    warning("No overlapping samples for region: ", region)
    next  # Skip to the next iteration if none found
  }
  
  # Perform Hellinger transformation on the OTU table
  otu_hellinger <- decostand(otu_table_df, method = "hellinger")
  
  # Extract the corresponding environmental data subset
  env_subset <- env_data[rownames(otu_hellinger), , drop = FALSE]
  
  # Calculate VIF scores using an initial CCA model
  cca_initial <- cca(otu_hellinger ~ ., env_subset)
  vif_scores <- vif.cca(cca_initial)
  
  # Filter environmental variables based on VIF threshold (< 10)
  env_filtered <- env_subset[, vif_scores < 10, drop = FALSE]
  
  if(ncol(env_filtered) == 0){
    warning("No environmental variables left after VIF filtering for region: ", region)
    next  # Skip to next iteration if none left
  }
  # Run CCA
  cca_result <- cca(otu_hellinger ~ ., env_filtered)
  variance <- round(cca_result$CCA$eig / sum(cca_result$CCA$eig) * 100, 1)
  x_label <- paste0("CCA1 (", variance[1], "%)")
  y_label <- paste0("CCA2 (", variance[2], "%)")
  
  # Permutation test
  cca_test <- anova(cca_result, permutations = 999)
  
  # Save plot to file
  file_path <- file.path(plots_dir, paste0("CCA_plot_XRFCu_orig_", region, ".png"))
  png(file_path)
  
  # --- Enhanced Plot ---
  plot(cca_result, type = "n", main = paste("CCA for", region), 
       xlab = x_label, ylab = y_label, scaling = 2)
  
  # Add samples (points and labels)
  points(cca_result, display = "sites", pch = 16, col = "blue", cex = 0.8)
  text(cca_result, display = "sites", labels = rownames(otu_hellinger), 
       col = "blue", cex = 0.6, pos = 3)
  
  # Add environmental vectors
  text(cca_result, display = "bp", col = "red", cex = 0.8, font = 2)
  
  # --- Add taxa labels CLEANLY ---
  taxa_scores <- scores(cca_result, display = "species")
  
  # Filter taxa: Only label those with strong contributions (top 10% from origin)
  taxa_dist <- sqrt(taxa_scores[,1]^2 + taxa_scores[,2]^2)
  top_taxa <- names(taxa_dist)[taxa_dist > quantile(taxa_dist, 0.90)]  # Top 10%
  
  # Use orditorp() for non-overlapping labels (vegan)
  orditorp(cca_result, display = "species", labels = top_taxa, 
           col = "darkgray", cex = 0.7, air = 1.2, priority = taxa_dist)
  
  # Add XRF contour
  ordisurf(cca_result, env_filtered$XRF, add = TRUE, col = "darkgreen", labcex = 0.8)

  legend("topleft",
         legend = c(
           paste("Explained variance (CCA1 + CCA2):",
                 round(sum(cca_result$CCA$eig[1:2]/cca_result$tot.chi*100), 1), "%"),
           paste("Model p-value:", cca_test$`Pr(>F)`[1])
         ),
         bty = "n", cex = 0.8)
  
  # Close the plotting device
  dev.off()
}
summary(cca_result)
anova(cca_result)