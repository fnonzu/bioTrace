# Load Required Packages ------------------------------------------------------
library(phyloseq)
library(tidyverse)
library(RColorBrewer)

# Constants and Configuration -------------------------------------------------
ABUNDANCE_THRESHOLD <- 0.01
V_REGIONS <- c("V1V2", "V2V3", "V3V4", "V4V5", "V5V7", "V7V9")

# Expected Reference Data -----------------------------------------------------
EXPECTED_GENERA <- tibble::tribble(
  ~Genus, ~Expected_Relative_Abundance,
  "Pseudomonas",             4.2,
  "Escherichia/Shigella",    10.1,
  "Salmonella",              10.4,
  "Lactobacillus",           18.4,
  "Enterococcus",            9.9,
  "Staphylococcus",          15.5,
  "Listeria",                14.1,
  "Bacillus",                17.4,
  "Saccharomyces",           0,
  "Cryptococcus",            0
) %>%
  mutate(MeanAbundance = Expected_Relative_Abundance / sum(Expected_Relative_Abundance))

GRAM_STAIN_INFO <- tibble::tribble(
  ~Genus, ~Gram_Stain,
  "Pseudomonas",          "Gram-negative",
  "Escherichia/Shigella", "Gram-negative",
  "Salmonella",           "Gram-negative",
  "Lactobacillus",        "Gram-positive",
  "Enterococcus",         "Gram-positive",
  "Staphylococcus",       "Gram-positive",
  "Listeria",             "Gram-positive",
  "Bacillus",             "Gram-positive",
  "Saccharomyces",        "Fungi",
  "Cryptococcus",         "Fungi"
)

# Data Processing Functions ---------------------------------------------------
process_region_data <- function(phy_region, threshold = ABUNDANCE_THRESHOLD) {
  # Transform to relative abundance
  phy_rel <- transform_sample_counts(phy_region, function(x) x / sum(x))
  
  # Aggregate at genus level
  phy_glom <- tax_glom(phy_rel, taxrank = "Genus", NArm = FALSE)
  
  # Melt and clean data
  phy_melt <- psmelt(phy_glom) %>%
    mutate(Genus = coalesce(Genus, "Unassigned"))
  
  # Identify low abundance taxa
  low_abundance_taxa <- phy_melt %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    filter(MeanAbundance < threshold) %>%
    pull(Genus)
  
  # Apply your requested ifelse modification
  phy_melt <- phy_melt %>%
    mutate(Genus = ifelse(Genus %in% low_abundance_taxa, "Other", Genus)) %>%
    group_by(Sample, Genus) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop")
  
  return(phy_melt)
}

process_low_abundance <- function(abundance_df, threshold) {
  abundance_df %>%
    mutate(Genus = if_else(MeanAbundance < threshold, "Other", Genus)) %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = sum(MeanAbundance), .groups = "drop")
}

# Main Processing Pipeline -----------------------------------------------------
# Extract Zymo Community samples
zymo_comm_regions <- map(V_REGIONS, ~ {
  region_samples <- phy_regions[[.x]] %>%
    subset_samples(grepl("Zymo_Comm", Sample)) %>%
    prune_samples(sample_sums(.) > 0, .)
  
  if (nsamples(region_samples) > 0) region_samples else NULL
}) %>% set_names(V_REGIONS)

# Process all regions
observed_abundance <- compact(zymo_comm_regions) %>%
  imap_dfr(~ {
    process_region_data(.x) %>%
      mutate(V_region = .y)
  })

# Prepare combined dataset
combined_abundance <- bind_rows(
  observed_abundance,
  EXPECTED_GENERA %>%
    select(Genus, MeanAbundance) %>%
    mutate(V_region = "Theoretical")
) %>%
  mutate(V_region = factor(V_region, levels = c(V_REGIONS, "Theoretical")))

# Visualization Functions ------------------------------------------------------
create_taxonomic_plot <- function(data) {
  custom_pal <- colorRampPalette(brewer.pal(8, "Set3"))(n_distinct(data$Genus))
  
  ggplot(data, aes(x = V_region, y = MeanAbundance, fill = Genus)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = custom_pal) +
    labs(title = "Taxonomic Composition by V Region",
         x = "Amplification Region", 
         y = "Relative Abundance") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5))
}

create_gram_stain_plot <- function(data) {
  data %>%
    left_join(GRAM_STAIN_INFO, by = "Genus") %>%
    mutate(Gram_Stain = coalesce(Gram_Stain, "Other")) %>%
    group_by(V_region, Gram_Stain) %>%
    summarise(TotalAbundance = sum(MeanAbundance), .groups = "drop") %>%
    ggplot(aes(x = V_region, y = TotalAbundance, fill = Gram_Stain)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Gram Stain Composition by V Region",
         x = "Amplification Region",
         y = "Total Relative Abundance") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
}

# Generate Plots --------------------------------------------------------------
taxonomic_plot <- create_taxonomic_plot(combined_abundance)
gram_plot <- create_gram_stain_plot(combined_abundance)

print(taxonomic_plot)
print(gram_plot)