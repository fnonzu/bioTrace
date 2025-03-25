# Function to inspect the structure of each region
explore_region_structure <- function(phy_regions) {
  for (region in names(phy_regions)) {
    cat("\n--- Region:", region, "---\n")
    phy_obj <- phy_regions[[region]]
    cat("Number of Samples:", nsamples(phy_obj), "\n")
    cat("Number of Taxa:", ntaxa(phy_obj), "\n")
  }
}

# Split V regions
V_regions <- c("V1V2", "V2V3", "V3V4", "V4V5", "V5V7", "V7V9")

phy_regions <- list()

for (region in V_regions) {
  # Select samples matching the current V region
  region_samples <- sample_data_df %>%
    filter(grepl(region, Sample)) %>%
    pull(Sample)
  
  # Subset the phyloseq object to include only samples from the current V region
  phy_region <- subset_samples(phy_obj_filtered, Sample %in% region_samples)
  
  # Remove taxa with zero counts
  phy_region <- prune_taxa(taxa_sums(phy_region) > 0, phy_region)
  
  # Add the phyloseq object to the list
  phy_regions[[region]] <- phy_region
}

# Explore the structure of each region
explore_region_structure(phy_regions)


### CLEANED


phy_regions_filtered <- list()

for (region in V_regions) {
  # Select samples matching the current V region from the filtered sample data
  region_samples <- sample_data_df_filtered %>%
    filter(grepl(region, Sample)) %>%
    pull(Sample)
  
  # Subset the filtered phyloseq object to include only samples from the current V region
  phy_region <- subset_samples(phy_obj_filtered, Sample %in% region_samples)
  
  # Remove taxa with zero counts
  phy_region <- prune_taxa(taxa_sums(phy_region) > 0, phy_region)
  
  # Add the phyloseq object to the list
  phy_regions_filtered[[region]] <- phy_region
}

# Explore the structure of each region
explore_region_structure(phy_regions_filtered)


