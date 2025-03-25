# SUBSET samples 3-8 
####

phy_regions_filtered_subset <- lapply(phy_regions_filtered, function(ps_obj) {
  sample_data(ps_obj)$sample_num <- as.numeric(sub("_.*", "", sample_data(ps_obj)$description))
  subset_samples(ps_obj, sample_num %in% 3:8)
})

############################
### RAREFACTION CURVE #####
#############################

# Rarefaction curve using ggrare()

library(phyloseq)
library(vegan)
library(ggplot2)


# USE phy_regions_filtered or phy_regions_filtered_subset depending on samples needed (1-9 or 3-8)
for (region in V_regions) {
  # 1. Create a copy of the phyloseq object for the current region
  ps_copy <- phy_regions_filtered_subset[[region]]
  
  # 2. Extract the OTU table as a matrix
  otu_mat <- as(otu_table(ps_copy), "matrix")
  
  # 3. Check if transposition is needed
  if (taxa_are_rows(ps_copy)) {
    otu_mat <- t(otu_mat)
  }
  
  
  # 4. Plot the rarefaction curve
  #  save plots
  file_name <- paste0("plots/rarefaction_curve_subset_", region, ".png")
  png(file_name, width = 800, height = 600)
  
  rarecurve(otu_mat, step = 100, sample = min(rowSums(otu_mat)), 
            col = rainbow(nrow(otu_mat)), label = FALSE, 
            main = paste("Rarefaction Curve", region), 
            xlab = "Number of Sequences", ylab = "Observed Species")
  
  # Extract legend labels from the current region's sample data
  legend_labels <- ps_copy@sam_data[["Sample"]]
  
  # Add a legend with smaller text in the bottom right
  legend("bottomright", 
         legend = legend_labels, 
         col = rainbow(length(legend_labels)), 
         lty = 1, 
         bty = "n", 
         cex = 0.8,   # Adjust text size
         pt.cex = 0.7 # Adjust point symbol size
  )
  
  # Close the device to save the plot file
  dev.off()
}



############################
### RAREFACTION CURVE REMOVE#####
#############################

# Rarefaction curve using ggrare()

library(phyloseq)
library(vegan)
library(ggplot2)

# Initialize list to store excluded samples
excluded_samples <- list()

# USE phy_regions_filtered or phy_regions_filtered_subset depending on samples needed (1-9 or 3-8)
for (region in V_regions) {
  # 1. Create a copy of the phyloseq object for the current region
  ps_copy <- phy_regions_filtered_subset[[region]]
  
  # 2. Identify and exclude samples with <1000 reads
  sample_read_sums <- sample_sums(ps_copy)
  samples_to_remove <- names(which(sample_read_sums < 2000))
  excluded_samples[[region]] <- samples_to_remove
  
  # Print excluded samples for current region
  if (length(samples_to_remove) > 0) {
    message("Region ", region, " excluded samples: ", paste(samples_to_remove, collapse = ", "))
  } else {
    message("Region ", region, " had no samples below 5000 reads")
  }
  
  # Filter the phyloseq object
  ps_filtered <- prune_samples(sample_read_sums >= 2000, ps_copy)
  
  # 3. Extract the OTU table as a matrix
  otu_mat <- as(otu_table(ps_filtered), "matrix")
  
  # 4. Check if transposition is needed
  if (taxa_are_rows(ps_filtered)) {
    otu_mat <- t(otu_mat)
  }
  
  # Skip region if no samples remain
  if (nrow(otu_mat) == 0) {
    message("Skipping region ", region, " - no samples remaining after filtering")
    next
  }
  
  # 5. Plot the rarefaction curve
  file_name <- paste0("plots/rarefaction_curve_subset_filtered_", region, ".png")
  png(file_name, width = 800, height = 600)
  
  rarecurve(otu_mat, step = 100, sample = min(rowSums(otu_mat)), 
            col = rainbow(nrow(otu_mat)), label = FALSE, 
            main = paste("Rarefaction Curve", region), 
            xlab = "Number of Sequences", ylab = "Observed Species")
  
  # Extract legend labels from the filtered samples
  legend_labels <- ps_filtered@sam_data[["Sample"]]
  
  # Add a legend with smaller text in the bottom right
  legend("bottomright", 
         legend = legend_labels, 
         col = rainbow(length(legend_labels)), 
         lty = 1, 
         bty = "n", 
         cex = 0.8,   # Adjust text size
         pt.cex = 0.7 # Adjust point symbol size
  )
  
  # Close the device to save the plot file
  dev.off()
}

# Print summary of excluded samples
cat("\nSummary of excluded samples:\n")
for (region in names(excluded_samples)) {
  cat(paste0(region, ": ", paste(excluded_samples[[region]], collapse = ", "), "\n"))
}



########## EXCLUDED SAMPLES PHY OBJ############

# Create new filtered phyloseq object list
phy_final <- list()

# Initialize list to track excluded samples
excluded_samples <- list()

for(region in names(phy_regions_filtered_subset)) {
  # Get region-specific phyloseq object
  ps_region <- phy_regions_filtered_subset[[region]]
  
  # Calculate sample read counts
  sample_reads <- sample_sums(ps_region)
  
  # Identify samples to keep (≥2000 reads)
  keep_samples <- names(which(sample_reads > 2000))
  
  # Track excluded samples
  excluded_samples[[region]] <- setdiff(sample_names(ps_region), keep_samples)
  
  # Filter the phyloseq object
  ps_filtered <- prune_samples(keep_samples, ps_region)
  
  # Remove any unused taxa if needed
  ps_filtered <- filter_taxa(ps_filtered, function(x) sum(x) > 0, TRUE)
  
  # Add to final object only if samples remain
  if(nsamples(ps_filtered) > 0) {
    phy_final[[region]] <- ps_filtered
  }
}

# Optional: Remove empty regions from phy_final
phy_final <- phy_final[!sapply(phy_final, is.null)]

# Print exclusion report
cat("Samples excluded per region:\n")
for(region in names(excluded_samples)) {
  if(length(excluded_samples[[region]]) > 0) {
    cat(paste0(
      region, ": ", 
      length(excluded_samples[[region]]), 
      " samples excluded - ",
      paste(excluded_samples[[region]], collapse = ", "),
      "\n"
    ))
  } else {
    cat(paste0(region, ": No samples excluded\n"))
  }
}

# Check sample counts before/after
original_counts <- sapply(phy_regions_filtered_subset, nsamples)
filtered_counts <- sapply(phy_final, nsamples)

data.frame(
  Region = names(original_counts),
  Original = original_counts,
  Filtered = ifelse(names(original_counts) %in% names(filtered_counts),
                    filtered_counts[match(names(original_counts), names(filtered_counts))],
                    0)
)




###EXAMPLE FOR ONE SAMPLE< THE FUNCITON IS BELOW
# Calculate alpha diversity indices
alpha_div <- estimate_richness(phy_regions_filtered_subset[["V3V4"]], measures = c("Observed", "Chao1", "Shannon"))
head(alpha_div)  # view the first few rows of the diversity estimates

plot_richness(phy_regions_filtered_subset[["V5V7"]],
              x = "sample_num", 
              measures = c("Shannon", "InvSimpson", "Chao1")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Sample Number") +
  scale_x_continuous(breaks = 3:8, labels = paste("Sample", 3:8))


# Calculate beta diversity using the Bray-Curtis dissimilarity metric
bray_dist <- distance(phy_regions_filtered_subset[["V3V4"]], method = "bray")

# Perform NMDS ordination based on the Bray-Curtis distance matrix
nmds_ord <- ordinate(phy_regions_filtered_subset[["V3V4"]], method = "NMDS", distance = bray_dist)

# Plot the NMDS ordination, coloring samples by SampleType
plot_ordination(phy_regions_filtered_subset[["V3V4"]], nmds_ord, color = "description") +
  theme_minimal()




##### FUCNTIONS FOR APLHA AND BETA DIVERSITY ####
library(phyloseq)
library(ggplot2)
library(gridExtra)  # For arranging multiple plots

### APLHA ###
plot_all_alpha_diversity <- function(phy_list) {
  plots <- list()
  
  for (region in names(phy_list)) {
    # Define the file name dynamically for each region
    file_name <- paste0("plots/alpha_diversity_final_", region, ".png")
    
    # Open a PNG device to save the plot
    png(file_name, width = 800, height = 600)
    
    p <- plot_richness(phy_list[[region]],
                       x = "sample_num", 
                       measures = c("Shannon", "InvSimpson", "Chao1")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      xlab("Sample Number") +
      scale_x_continuous(breaks = 3:8, labels = paste("Sample", 3:8)) +
      ggtitle(paste("Alpha Diversity -", region))
    
    plots[[region]] <- p
    print(p)  # Print each plot
    dev.off()
  }
  
  return(plots)
}




### ALPHA SHANNON FOR ONE PLOT
### ALPHA DIVERSITY - SHANNON COMPARISON ###
plot_combined_shannon <- function(phy_list) {
  # Create combined dataframe
  combined_df <- data.frame()
  
  # Collect data from all regions
  for (region in names(phy_list)) {
    ps <- phy_list[[region]]
    if (nsamples(ps) > 0) {  # Skip empty regions
      # Calculate Shannon diversity
      div_df <- estimate_richness(ps, measures = "Shannon")
      
      # Add metadata and region information
      div_df$Sample <- sample_names(ps)
      div_df$Region <- region
      
      combined_df <- rbind(combined_df, div_df)
    }
  }
  
  # Create unified plot with fixed y-axis
  p <- ggplot(combined_df, aes(x = Region, y = Shannon, color = Region)) +
    geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.3) +
    scale_y_continuous(limits = c(0, max(combined_df$Shannon) + 0.5)) +
    labs(title = "Shannon Diversity Across V Regions",
         x = "Hypervariable Region",
         y = "Shannon Diversity Index") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.grid.major.x = element_blank())
  
  # Save and return plot
  ggsave("plots/combined_shannon_diversity.png", p, width = 10, height = 6)
  return(p)
}

# Usage example:
plot_combined_shannon(phy_final)





### BETA ###
plot_all_beta_diversity <- function(phy_list) {
  plots <- list()
  
  for (region in names(phy_list)) {
    # Define the file name dynamically for each region
    file_name <- paste0("plots/beta_diversity_final_", region, ".png")
    
    # Open a PNG device to save the plot
    png(file_name, width = 800, height = 600)
    
    # Compute Bray-Curtis distance
    bray_dist <- distance(phy_list[[region]], method = "bray")
    
    # Perform NMDS ordination
    nmds_ord <- ordinate(phy_list[[region]], method = "NMDS", distance = bray_dist)
    
    # Plot NMDS ordination
    p <- plot_ordination(phy_list[[region]], nmds_ord, color = "description") +
      theme_minimal() +
      ggtitle(paste("Beta Diversity (NMDS) -", region))
    
    plots[[region]] <- p
    print(p)  # Print each plot
    dev.off()
  }
  
  return(plots)
}

### PLOT
# alpha_plots <- plot_all_alpha_diversity(phy_regions_filtered_subset)
# beta_plots <- plot_all_beta_diversity(phy_regions_filtered_subset)

alpha_plots <- plot_all_alpha_diversity(phy_final)
beta_plots <- plot_all_beta_diversity(phy_final)

##### ADDITIONAL #####
#TO deal with single plot and explore
# Exclude sample
phy_filtered <- subset_samples(phy_regions_filtered_subset[["V3V4"]], Sample != "3_2_V3V4")
bray_dist <- distance(phy_filtered, method = "bray")
nmds_ord <- ordinate(phy_filtered, method = "NMDS", distance = bray_dist)
plot_ordination(phy_filtered, nmds_ord, color = "description") +
  theme_minimal()

# Use few samples -- doesn't make any sense
phy_filtered <- subset_samples(phy_regions_filtered_subset[["V3V4"]], Sample %in% c("3_2_V3V4", "3_1_V3V4", "3_3_V3V4"))
bray_dist <- distance(phy_filtered, method = "bray")
nmds_ord <- ordinate(phy_filtered, method = "NMDS", distance = bray_dist)
plot_ordination(phy_filtered, nmds_ord, color = "description") +
  theme_minimal()



##### TAXONOMIC COMPOSITION 
# 
# # Agglomerate taxa at the Genus level
# phy_genus <- tax_glom(phy_regions_filtered_subset[["V5V7"]], taxrank = "Genus")
# # Transform counts to relative abundance
# phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
# 
# # Summarize relative abundances per genus
# genus_abundance <- psmelt(phy_genus_rel) %>%
#   group_by(Genus) %>%
#   summarise(MeanAbundance = mean(Abundance)) %>%
#   arrange(desc(MeanAbundance))
# 
# # View top 10 genera
# head(genus_abundance, 10)
# # Plot taxonomic composition at Genus level
# ggplot(psmelt(phy_genus_rel), aes(x = sample_num, y = Abundance, fill = Genus)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   labs(title = "Genus-Level Taxonomic Composition",
#        x = "Sample Number",
#        y = "Relative Abundance") +
#   scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
# # Keep only top 10 most abundant genera, others become "Other"
# top_genera <- genus_abundance$Genus[1:10]
# psmelt_data <- psmelt(phy_genus_rel) %>%
#   mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other"))
# 
# # Stacked bar plot
# ggplot(psmelt_data, aes(x = sample_num, y = Abundance, fill = Genus)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   labs(title = "Genus-Level Taxonomic Composition (Top 10 Genera)",
#        x = "Sample Number",
#        y = "Relative Abundance") +
#   scale_fill_manual(values = RColorBrewer::brewer.pal(11, "Paired"))
# 

### SAVING THE PLOTS FUNCTION

### FULL 1-9

# Define color palette function for 16 categories
get_custom_palette <- function() {
  colorRampPalette(RColorBrewer::brewer.pal(18, "Paired"))(16)
}
# Loop through each region
for (region in V_regions) {
  # Agglomerate taxa at Genus level
  phy_genus <- tax_glom(phy_regions_filtered[[region]], taxrank = "Genus")
  
  # Transform to relative abundance
  phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
  
  # Summarize relative abundances
  genus_abundance <- psmelt(phy_genus_rel) %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    arrange(desc(MeanAbundance))
  
  # Get top 15 genera
  top_genera <- genus_abundance$Genus[1:15]
  
  # Prepare data for plotting
  psmelt_data <- psmelt(phy_genus_rel) %>%
    mutate(Genus = ifelse(Genus %in% top_genera & !is.na(Genus), 
                          Genus, "Other"))
  
  # Create plot
  genus_plot <- ggplot(psmelt_data, 
                       aes(x = Sample_orig, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1),
          legend.position = "bottom") +
    labs(title = paste("Top 15 Genera -", region),
         x = "Sample Origin",
         y = "Relative Abundance") +
    scale_fill_manual(values = get_custom_palette()) +
    guides(fill = guide_legend(ncol = 3))
  
  # Save plot
  ggsave(paste0("plots/", region, "_top15_genera_1-9samples.png"), 
         plot = genus_plot,
         width = 12, 
         height = 8,
         dpi = 300)
}

### SUBSET 3-8
# Define color palette function for 16 categories
get_custom_palette <- function() {
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(16)
}
# Loop through each region
for (region in V_regions) {
  # Agglomerate taxa at Genus level
  phy_genus <- tax_glom(phy_regions_filtered_subset[[region]], taxrank = "Genus")
  
  # Transform to relative abundance
  phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
  
  # Summarize relative abundances
  genus_abundance <- psmelt(phy_genus_rel) %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    arrange(desc(MeanAbundance))
  
  # Get top 15 genera
  top_genera <- genus_abundance$Genus[1:15]
  
  # Prepare data for plotting
  psmelt_data <- psmelt(phy_genus_rel) %>%
    mutate(Genus = ifelse(Genus %in% top_genera & !is.na(Genus), 
                          Genus, "Other"))
  
  # Create plot
  genus_plot <- ggplot(psmelt_data, 
                       aes(x = Sample_orig, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1),
          legend.position = "bottom") +
    labs(title = paste("Top 15 Genera -", region),
         x = "Sample Origin",
         y = "Relative Abundance") +
    scale_fill_manual(values = get_custom_palette()) +
    guides(fill = guide_legend(ncol = 3))
  
  # Save plot
  ggsave(paste0("plots/", region, "_top15_genera.png"), 
         plot = genus_plot,
         width = 12, 
         height = 8,
         dpi = 300)
}




##########################################################
###################### FINAL TAXONOMIC COMPOSITION #######

### SUBSET 3-8
# Define color palette function for 16 categories
get_custom_palette <- function() {
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(16)
}
# Loop through each region
for (region in V_regions) {
  file_name <- paste0("plots/", "top15_genera_final_", region, ".png")
  
  # Open a PNG device to save the plot
  png(file_name, width = 800, height = 600)
  
  # Agglomerate taxa at Genus level
  phy_genus <- tax_glom(phy_final[[region]], taxrank = "Genus")

  # Transform to relative abundance
  phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
  
  # Summarize relative abundances
  genus_abundance <- psmelt(phy_genus_rel) %>%
    group_by(Genus) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    arrange(desc(MeanAbundance))
  
  # Get top 15 genera
  top_genera <- genus_abundance$Genus[1:15]
  
  # Prepare data for plotting
  psmelt_data <- psmelt(phy_genus_rel) %>%
    mutate(Genus = ifelse(Genus %in% top_genera & !is.na(Genus), 
                          Genus, "Other"))
  
  # Create plot
  genus_plot <- ggplot(psmelt_data, 
                       aes(x = Sample_orig, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1),
          legend.position = "bottom") +
    labs(title = paste("Top 15 genera -", region),
         x = "Sample origin",
         y = "Relative Abundance") +
    scale_fill_manual(values = get_custom_palette()) +
    guides(fill = guide_legend(ncol = 3))
  
  print(genus_plot)

  dev.off()
}







## TOP 15 taxa

# Agglomerate taxa at the Genus level
phy_genus <- tax_glom(phy_regions_filtered_subset[["V5V7"]], taxrank = "Genus")

# Transform counts to relative abundance
phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))

# Summarize relative abundances per genus
genus_abundance <- psmelt(phy_genus_rel) %>%
  group_by(Genus) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  arrange(desc(MeanAbundance))

# View top 15 genera
head(genus_abundance, 15)

# Plot taxonomic composition at Genus level
ggplot(psmelt(phy_genus_rel), aes(x = sample_num, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Genus-Level Taxonomic Composition",
       x = "Sample",
       y = "Relative Abundance") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2"))

# Keep only top 15 most abundant genera; all others are labeled "Other"
top_genera <- genus_abundance$Genus[1:15]
psmelt_data <- psmelt(phy_genus_rel) %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other"))

# Create a custom color palette for 16 groups (15 genera + "Other")
custom_palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(16)

# Stacked bar plot for top 15 genera
ggplot(psmelt_data, aes(x = sample_num, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Genus-Level Taxonomic Composition (Top 15 Genera)",
       x = "Sample",
       y = "Relative Abundance") +
  scale_fill_manual(values = custom_palette)


### PHYLUM


## TOP 15 taxa

# Agglomerate taxa at the Genus level
phy_genus <- tax_glom(phy_regions_filtered_subset[["V3V4"]], taxrank = "Phylum")

# Transform counts to relative abundance
phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))

# Summarize relative abundances per genus
genus_abundance <- psmelt(phy_genus_rel) %>%
  group_by(Phylum) %>%
  summarise(MeanAbundance = mean(Abundance)) %>%
  arrange(desc(MeanAbundance))

# View top 15 genera
head(genus_abundance, 15)

# Plot taxonomic composition at Genus level
ggplot(psmelt(phy_genus_rel), aes(x = sample_num, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Phylum-Level Taxonomic Composition",
       x = "Sample",
       y = "Relative Abundance") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2"))

# Keep only top 15 most abundant genera; all others are labeled "Other"
top_genera <- genus_abundance$Phylum[1:15]
psmelt_data <- psmelt(phy_genus_rel) %>%
  mutate(Phylum = ifelse(Phylum %in% top_genera, Phylum, "Other"))

# Create a custom color palette for 16 groups (15 genera + "Other")
custom_palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(16)

# Stacked bar plot for top 15 genera
ggplot(psmelt_data, aes(x = sample_num, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Phylum-Level Taxonomic Composition (Top 15 Genera)",
       x = "Sample",
       y = "Relative Abundance") +
  scale_fill_manual(values = custom_palette)
