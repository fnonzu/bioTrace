library(phyloseq)
library(ggplot2)
library(dplyr)

# Input/Output Documentation
# IN_DATA:  Phyloseq object path (.RDS)
# OUT_DATA: phy_obj (original data), 
#           phy_obj_filtered (filtered shooting range samples),
#           phy_regions (isolated V regions)
#           phy_regions_filtered (only pollution V regions)

# ---------------------------
# Function to exlore phyloseq object
# ---------------------------

explore_phyloseq_structure <- function(phy_obj) {
  message("\n=== Phyloseq Object Structure ===\n")
  print(phy_obj)
  
  message("\n=== Sample Data Preview ===\n")
  print(head(sample_data(phy_obj)))
}

explore_region_structure <- function(phy_regions) {
  for (region in names(phy_regions)) {
    cat("\n--- Region:", region, "---\n")
    phy_obj <- phy_regions[[region]]
    cat("Number of Samples:", nsamples(phy_obj), "\n")
    cat("Number of Taxa:", ntaxa(phy_obj), "\n")
  }
}
# ---------------------------
# Data Import and Preparation
# ---------------------------

# Import phyloseq object
phy_obj <- readRDS("~/p952/16S_dada2_pipeline/results/16S_phyloseq_object.RDS")

# Correct naming error in sample descriptions
sample_data(phy_obj)$description <- gsub(
  "Hardkern", 
  "Habkern", 
  sample_data(phy_obj)$description
)

# ---------------------------
# Data Filtering
# ---------------------------

# Filter samples to include only pollution gradient samples
phy_obj_filtered <- subset_samples(
  phy_obj,
  test_name == "pollution_gradient"
)

# ---------------------------
# Data Exploration
# ---------------------------

explore_phyloseq_structure(phy_obj)

# ---------------------------
# V Region Extraction
# ---------------------------

v_regions <- c("V1V2", "V2V3", "V3V4", "V4V5", "V5V7", "V7V9")

phy_regions <- list()

for (region in v_regions) {
  # Identify samples for current region
  region_samples <- grep(
    pattern = region,
    x = sample_names(phy_obj),
    value = TRUE
  )
  
  # Create subset for current region
  phy_region <- subset_samples(
    phy_obj,
    sample_names(phy_obj) %in% region_samples
  )
  
  # Remove empty taxa
  phy_region <- prune_taxa(taxa_sums(phy_region) > 0, phy_region)
  
  # Store processed object
  phy_regions[[region]] <- phy_region
}

phy_regions_filtered <- list()

for (region in v_regions) {
  # Identify samples for current region
  region_samples <- grep(
    pattern = region,
    x = sample_names(phy_obj_filtered),
    value = TRUE
  )
  
  # Create subset for current region
  phy_region <- subset_samples(
    phy_obj_filtered,
    sample_names(phy_obj_filtered) %in% region_samples
  )
  
  # Remove empty taxa
  phy_region <- prune_taxa(taxa_sums(phy_region) > 0, phy_region)
  
  # Store processed object
  phy_regions_filtered[[region]] <- phy_region
}

explore_region_structure(phy_regions_filtered)


# ---------------------------
# V Region Subsetting (samples 3-8)
# ---------------------------


phy_regions_filtered_subset <- lapply(phy_regions_filtered, function(ps_obj) {
  sample_data(ps_obj)$sample_num <- as.numeric(sub("_.*", "", sample_data(ps_obj)$description))
  subset_samples(ps_obj, sample_num %in% 3:8)
})
