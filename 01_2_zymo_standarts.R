library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

##################
# this script 
##################

# Create a list to hold Zymo_Comm samples for each region
zymo_comm_regions <- list()

for (region in v_regions) {
  # Get the phyloseq object for the region
  phy_region <- phy_regions[[region]]

  # Subset to Zymo_Comm samples within the region
  zymo_comm_samples <- subset_samples(phy_region, grepl("Zymo_Comm", Sample))

  # Check if there are Zymo_Comm samples in the region
  if (nsamples(zymo_comm_samples) > 0) {
    # Normalize the data
    zymo_comm_rel <- transform_sample_counts(zymo_comm_samples, function(x) x / sum(x))

    # Add to the list
    zymo_comm_regions[[region]] <- zymo_comm_rel
  }
}

# Expected abundances for the ZymoBIOMICS Microbial Community Standard at genus level
expected_abundance <- data.frame(
  Genus = c(
    "Pseudomonas",
    "Escherichia/Shigella",
    "Salmonella",
    "Lactobacillus",
    "Enterococcus",
    "Staphylococcus",
    "Listeria",
    "Bacillus",
    "Saccharomyces",
    "Cryptococcus"
  ),
  Expected_Relative_Abundance = c(4.2, 10.1, 10.4, 18.4, 9.9, 15.5, 14.1, 17.4, 0, 0)
)

# Convert percentages to proportions
expected_abundance$Expected_Relative_Abundance <- expected_abundance$Expected_Relative_Abundance / 100

# # Create a list to store observed abundances for each region
# observed_abundances <- list()
# 
# for (region in names(zymo_comm_regions)) {
#   # Get the normalized phyloseq object for the region
#   zymo_comm_rel <- zymo_comm_regions[[region]]
# 
#   # Handle missing genus assignments if necessary
#   tax_data <- as.data.frame(tax_table(zymo_comm_rel))
#   tax_data$Genus <- as.character(tax_data$Genus)
#   tax_data$Genus[is.na(tax_data$Genus) | tax_data$Genus == ""] <- "Unclassified_Genus"
#   tax_table(zymo_comm_rel) <- as.matrix(tax_data)
# 
#   # Agglomerate taxa at the genus level
#   zymo_comm_genus <- tax_glom(zymo_comm_rel, taxrank = "Genus")
# 
#   # Melt the phyloseq object to long format
#   zymo_comm_melt <- psmelt(zymo_comm_genus)
# 
#   # Filter to genera of interest
#   zymo_comm_melt <- zymo_comm_melt %>%
#     filter(Genus %in% expected_abundance$Genus)
# 
#   # Calculate mean observed relative abundance for each genus
#   observed_comm_abundance <- zymo_comm_melt %>%
#     group_by(Genus) %>%
#     summarize(Observed_Relative_Abundance = mean(Abundance))
# 
#   # Merge observed and expected data
#   comparison_comm_df <- merge(observed_comm_abundance, expected_abundance, by = "Genus", all = TRUE)
# 
#   # Replace NA with 0 for genera not detected
#   comparison_comm_df[is.na(comparison_comm_df)] <- 0
# 
#   # Store in the list
#   observed_abundances[[region]] <- comparison_comm_df
# }
# 
# # Create a list to store plots
# comparison_plots <- list()
# 
# ## Ploting the mean arelative abundance
# for (region in names(observed_abundances)) {
#   comparison_comm_df <- observed_abundances[[region]]
# 
#   # Prepare data for plotting
#   comparison_comm_long <- comparison_comm_df %>%
#     pivot_longer(cols = c("Observed_Relative_Abundance", "Expected_Relative_Abundance"),
#                  names_to = "Type",
#                  values_to = "Relative_Abundance")
# 
#   # Create the plot
#   plot <- ggplot(comparison_comm_long, aes(x = Genus, y = Relative_Abundance, fill = Type)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = paste("Observed vs. Expected Relative Abundance -", region),
#          y = "Relative Abundance",
#          x = "Genus") +
#     scale_fill_manual(values = c("Observed_Relative_Abundance" = "#1f77b4",
#                                  "Expected_Relative_Abundance" = "#ff7f0e"))
# 
#   # Save the plot in the list
#   comparison_plots[[region]] <- plot
# }
# 
# # Display the plots
# for (region in names(comparison_plots)) {
#   print(comparison_plots[[region]])
# }

# 
# 
# # ---- Taxonomic Composition ----
# 
# # Transform counts to relative abundances for V5V7
# phy_V1V2_zymo <- transform_sample_counts(phy_regions[["V1V2"]], function(x) x / sum(x))
# 
# # Plot taxonomic composition at the Phylum level
# tax_plot <- plot_bar(phy_V1V2_zymo, fill = "Genus") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Aggregated Taxonomic Composition for V1V2 Region", x = "Samples", y = "Relative Abundance")
# print(tax_plot)
# 
# 
# # Function to plot taxonomic composition at a specified taxonomic level
# plot_taxonomic_composition <- function(phy_obj, tax_rank, abundance_threshold = 0.01, title = NULL) {
#   # Transform counts to relative abundances
#   phy_rel <- transform_sample_counts(phy_obj, function(x) x / sum(x))
#   
#   # Agglomerate taxa at the specified taxonomic rank
#   phy_glom <- tax_glom(phy_rel, taxrank = tax_rank, NArm = FALSE)
#   
#   # Melt the phyloseq object to long format
#   phy_melt <- psmelt(phy_glom)
#   
#   # Replace NA in taxonomic rank with "Unassigned"
#   phy_melt[[tax_rank]][is.na(phy_melt[[tax_rank]])] <- "Unassigned"
#   
#   # Calculate mean relative abundance of each taxon across all samples
#   taxon_abundance <- phy_melt %>%
#     group_by(!!sym(tax_rank)) %>%
#     summarise(MeanAbundance = mean(Abundance)) %>%
#     arrange(desc(MeanAbundance))
#   
#   # Identify taxa below the abundance threshold
#   low_abundance_taxa <- taxon_abundance %>%
#     filter(MeanAbundance < abundance_threshold) %>%
#     pull(!!sym(tax_rank))
#   
#   # Replace low abundance taxa with "Other"
#   phy_melt[[tax_rank]] <- ifelse(
#     phy_melt[[tax_rank]] %in% low_abundance_taxa,
#     "Other",
#     phy_melt[[tax_rank]]
#   )
#   
#   # Recalculate mean relative abundance after grouping "Other"
#   phy_melt <- phy_melt %>%
#     group_by(Sample, !!sym(tax_rank)) %>%
#     summarise(Abundance = sum(Abundance)) %>%
#     ungroup()
#   
#   # Create the plot
#   tax_plot <- ggplot(phy_melt, aes(x = Sample, y = Abundance, fill = !!sym(tax_rank))) +
#     geom_bar(stat = "identity") +
#     scale_fill_viridis_d() +
#     theme_minimal() +
#     labs(
#       title = ifelse(is.null(title), paste("Taxonomic Composition at", tax_rank, "Level"), title),
#       x = "Samples",
#       y = "Relative Abundance",
#       fill = tax_rank
#     ) +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       legend.position = "right"
#     )
#   
#   return(tax_plot)
# }
# 
# 
# genus_plot <- plot_taxonomic_composition(
#   phy_obj = zymo_comm_regions[["V1V2"]],
#   tax_rank = "Genus",
#   abundance_threshold = 0.01,  # Adjust threshold as needed
#   title = "Taxonomic Composition at Genus Level for V1V2 Region"
# )
# print(genus_plot)
# 


# printing all plots for V regions and the standart

# Initialize an empty list to store data frames
phy_melt_list <- list()

abundance_threshold <- 0.01  # Adjust as needed

for (V_region in v_regions) {
  phy_obj <- zymo_comm_regions[[V_region]]
  
  # Transform counts to relative abundances
  phy_rel <- transform_sample_counts(phy_obj, function(x) x / sum(x))
  
  # Agglomerate taxa at the Genus level
  phy_glom <- tax_glom(phy_rel, taxrank = "Genus", NArm = FALSE)
  
  # Melt the phyloseq object to long format
  phy_melt <- psmelt(phy_glom)
  
  # Replace NA in Genus with "Unassigned"
  phy_melt$Genus[is.na(phy_melt$Genus)] <- "Unassigned"
  
  # Calculate mean relative abundance of each taxon across all samples in that V region
  taxon_abundance <- phy_melt %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    arrange(desc(MeanAbundance))
  
  # Identify taxa below the abundance threshold
  low_abundance_taxa <- taxon_abundance %>%
    filter(MeanAbundance < abundance_threshold) %>%
    pull(Genus)
  
  # Replace low abundance taxa with "Other"
  phy_melt$Genus <- ifelse(
    phy_melt$Genus %in% low_abundance_taxa,
    "Other",
    phy_melt$Genus
  )
  
  # Recalculate mean relative abundance after grouping "Other"
  phy_melt <- phy_melt %>%
    group_by(Sample, Genus) %>%
    summarise(Abundance = sum(Abundance)) %>%
    ungroup()
  
  # Calculate mean abundance per Genus for the V region
  mean_abundance <- phy_melt %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    ungroup()
  
  # Add V_region as a column
  mean_abundance$V_region <- V_region
  
  # Append to the list
  phy_melt_list[[V_region]] <- mean_abundance
}

# Combine all data frames into one
Observed_abundance <- bind_rows(phy_melt_list)

# Expected abundances for the ZymoBIOMICS Microbial Community Standard at genus level
expected_abundance <- data.frame(
  Genus = c(
    "Pseudomonas",
    "Escherichia/Shigella",
    "Salmonella",
    "Lactobacillus",
    "Enterococcus",
    "Staphylococcus",
    "Listeria",
    "Bacillus",
    "Saccharomyces",
    "Cryptococcus"
  ),
  Expected_Relative_Abundance = c(4.2, 10.1, 10.4, 18.4, 9.9, 15.5, 14.1, 17.4, 0, 0)
)

# Normalize expected abundances to sum to 1
expected_abundance$MeanAbundance <- expected_abundance$Expected_Relative_Abundance / sum(expected_abundance$Expected_Relative_Abundance)

# Create expected abundance data frame
expected_abundance_adjusted <- expected_abundance %>%
  select(Genus, MeanAbundance) %>%
  mutate(V_region = "Expected")

# Combine observed and expected abundances
Combined_abundance <- bind_rows(Observed_abundance, expected_abundance)

# Replace NA in a specific column
Combined_abundance <- Combined_abundance %>%
  mutate(V_region = ifelse(is.na(V_region), "Standart", V_region))

# Create the plot
ggplot(Combined_abundance, aes(x = Genus, y = MeanAbundance, fill = V_region)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(
    title = "Mean Taxonomic Composition for Zymo standarts",
    x = "Genus",
    y = "Mean Relative Abundance",
    fill = "V Region"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# 
# 
# primer_efficiency <- observed_abundances %>%
#   bind_rows(.id = "Region") %>%
#   group_by(Region) %>%
#   summarise(
#     Pearson = cor(Observed_Relative_Abundance, Expected_Relative_Abundance),
#     RMSE = sqrt(mean((Observed_Relative_Abundance - Expected_Relative_Abundance)^2))
#   )
# 
# ggplot(primer_efficiency, aes(x = Region, y = Pearson, fill = RMSE)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   ggtitle("Primer Region Performance Metrics")

# Ensure Combined_abundance is correctly structured
Combined_abundance <- bind_rows(Observed_abundance, expected_abundance)
Combined_abundance$V_region[is.na(Combined_abundance$V_region)] <- "Theoretical"
Combined_abundance$V_region <- as.factor(Combined_abundance$V_region)

library(RColorBrewer)
custom_palette <- colorRampPalette(brewer.pal(8, "Set3"))(15)

ggplot(Combined_abundance, aes(x = V_region, y = MeanAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_palette) +
  labs(title = "Taxonomic Composition",
       x = "V Region",
       y = "Relative Abundance",
       fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Define Gram staining for expected genera
gram_stain <- data.frame(
  Genus = c("Pseudomonas", "Escherichia/Shigella", "Salmonella", 
            "Lactobacillus", "Enterococcus", "Staphylococcus", 
            "Listeria", "Bacillus", "Saccharomyces", "Cryptococcus"),
  Gram_Stain = c(rep("Gram-negative", 3), rep("Gram-positive", 5), rep("Fungi", 2))
)

# Merge with data and summarize
Combined_abundance_gram <- Combined_abundance %>%
  left_join(gram_stain, by = "Genus") %>%
  mutate(Gram_Stain = replace_na(Gram_Stain, "Other"))

gram_summary <- Combined_abundance_gram %>%
  group_by(V_region, Gram_Stain) %>%
  summarise(TotalAbundance = sum(MeanAbundance), .groups = 'drop')

# Plot Gram staining composition
ggplot(gram_summary, aes(x = V_region, y = TotalAbundance, fill = Gram_Stain)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Gram Stain Composition by V Region",
       x = "V Region",
       y = "Total Relative Abundance") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#########################
#### ABSOLUTE ########### DONT WORK
#########################
zymo_comm_rel <- zymo_comm_samples  # Keep raw counts
observed_comm_abundance <- zymo_comm_melt %>%
  group_by(Genus) %>%
  summarize(Observed_Absolute_Abundance = mean(Abundance))  # Mean raw count
# List to store alpha diversity results
alpha_results <- list()

for (region in names(zymo_comm_regions)) {
  phy_region <- zymo_comm_regions[[region]]
  
  # Calculate alpha diversity
  alpha_div <- estimate_richness(phy_region, measures = "Observed")
  alpha_div$Region <- region
  
  alpha_results[[region]] <- alpha_div
}

# Combine results
alpha_df <- bind_rows(alpha_results)

# Plot alpha diversity
ggplot(alpha_df, aes(x = Region, y = Observed)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity (Observed Genera) by V Region",
       y = "Number of Observed Genera") +
  theme_minimal()
# Prepare absolute abundance data
absolute_abundances <- bind_rows(lapply(names(observed_abundances), function(region) {
  df <- observed_abundances[[region]]
  df$Region <- region
  df
}))

# Plot absolute counts
ggplot(absolute_abundances, aes(x = Genus, y = Observed_Absolute_Abundance, fill = Region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean Absolute Abundance by V Region",
       y = "Mean Absolute Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
