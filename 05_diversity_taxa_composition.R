library(phyloseq)
library(ggplot2)
library(vegan)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

# Configuration ---------------------------------------------------------------
PLOT_DIR <- "plots"
DEFAULT_REGION <- "V5V7"
SAMPLE_RANGE <- 3:8

# Create plot directory
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR)

# Core Analysis Functions -----------------------------------------------------
analyze_alpha_diversity <- function(phy_obj, region = NULL, measures = c("Shannon", "InvSimpson", "Chao1")) {
  # Validate input
  if (is.null(region) && is.list(phy_obj)) {
    stop("For multiple regions, use plot_all_alpha_diversity()")
  }
  
  # Get target phyloseq object
  target_phy <- if (is.list(phy_obj)) phy_obj[[region]] else phy_obj
  
  # Calculate and plot diversity
  div_plot <- plot_richness(target_phy, x = "sample_num", measures = measures) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    xlab("Sample Number") +
    scale_x_continuous(breaks = SAMPLE_RANGE, labels = paste("Sample", SAMPLE_RANGE)) +
    ggtitle(paste("Alpha Diversity -", region))
  
  return(div_plot)
}

analyze_beta_diversity <- function(phy_obj, region = NULL, color_var = "description") {
  # Validate input
  if (is.null(region) && is.list(phy_obj)) {
    stop("For multiple regions, use plot_all_beta_diversity()")
  }
  
  # Get target phyloseq object
  target_phy <- if (is.list(phy_obj)) phy_obj[[region]] else phy_obj
  
  # Calculate diversity metrics
  bray_dist <- distance(target_phy, method = "bray")
  nmds_ord <- ordinate(target_phy, method = "NMDS", distance = bray_dist)
  
  # Create ordination plot
  ord_plot <- plot_ordination(target_phy, nmds_ord, color = color_var) +
    theme_minimal() +
    ggtitle(paste("Beta Diversity (NMDS) -", region))
  
  return(ord_plot)
}

# Batch Processing Functions --------------------------------------------------
plot_all_alpha_diversity <- function(phy_list) {
  plots <- lapply(names(phy_list), function(region) {
    plot_file <- file.path(PLOT_DIR, paste0("alpha_diversity_final_", region, ".png"))
    png(plot_file, width = 800, height = 600)
    p <- analyze_alpha_diversity(phy_list, region)
    print(p)
    dev.off()
    p
  })
  names(plots) <- names(phy_list)
  invisible(plots)
}

plot_all_beta_diversity <- function(phy_list) {
  plots <- lapply(names(phy_list), function(region) {
    plot_file <- file.path(PLOT_DIR, paste0("beta_diversity_final_", region, ".png"))
    png(plot_file, width = 800, height = 600)
    p <- analyze_beta_diversity(phy_list, region)
    print(p)
    dev.off()
    p
  })
  names(plots) <- names(phy_list)
  invisible(plots)
}

# Combined Visualization ------------------------------------------------------
plot_combined_shannon <- function(phy_list) {
  # Create combined dataframe
  combined_df <- do.call(rbind, lapply(names(phy_list), function(region) {
    ps <- phy_list[[region]]
    if (nsamples(ps) == 0) return(NULL)
    
    div_df <- estimate_richness(ps, measures = "Shannon")
    div_df$Sample <- sample_names(ps)
    div_df$Region <- region
    div_df
  }))
  
  # Create plot
  p <- ggplot(combined_df, aes(x = Region, y = Shannon, color = Region)) +
    geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.3) +
    labs(title = "Shannon Diversity Across V Regions",
         x = "Hypervariable Region",
         y = "Shannon Diversity Index") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          panel.grid.major.x = element_blank())
  
  # Save and return
  ggsave(file.path(PLOT_DIR, "combined_shannon_diversity.png"), p, width = 10, height = 6)
  p
}

# Sample Exclusion Helper -----------------------------------------------------
exclude_samples <- function(phy_obj, samples_to_remove) {
  prune_samples(!(sample_names(phy_obj) %in% samples_to_remove), phy_obj)
}

# Usage Examples --------------------------------------------------------------
if (FALSE) {  # Prevent execution when sourcing
  # Single region analysis (V5V7 example)
  analyze_alpha_diversity(phy_final[[DEFAULT_REGION]], region = DEFAULT_REGION)
  analyze_beta_diversity(phy_final[[DEFAULT_REGION]], region = DEFAULT_REGION)
  
  # Batch processing
  alpha_plots <- plot_all_alpha_diversity(phy_final)
  beta_plots <- plot_all_beta_diversity(phy_final)
  plot_combined_shannon(phy_final)
  
  # Sample exclusion example
  phy_filtered <- exclude_samples(phy_final[["V3V4"]], "3_2_V3V4")
  analyze_beta_diversity(phy_filtered, region = "V3V4")
}

#-----------------------------------------------------------
#--------------------------- Taxonomic composition ---------
#-----------------------------------------------------------

# plotting function -------------------------------------
plot_taxonomic_composition <- function(phy_obj, 
                                       tax_level = "Phylum",
                                       top_n = 15,
                                       palette_name = "Set2",
                                       plot_title = NULL,
                                       output_file = NULL) {
  # Validate taxonomic rank
  if (!tax_level %in% rank_names(phy_obj)) {
    stop(paste("Taxonomic level", tax_level, "not found in phyloseq object"))
  }
  
  # Process data
  phy_glom <- tax_glom(phy_obj, taxrank = tax_level)
  phy_rel <- transform_sample_counts(phy_glom, function(x) x / sum(x))
  
  # Get abundance data
  tax_abundance <- psmelt(phy_rel) %>%
    group_by(!!sym(tax_level)) %>%
    summarise(MeanAbundance = mean(Abundance)) %>%
    arrange(desc(MeanAbundance))
  
  # Identify top taxa
  top_taxa <- tax_abundance[[tax_level]][1:top_n]
  
  # Prepare plot data
  plot_data <- psmelt(phy_rel) %>%
    mutate(!!sym(tax_level) := ifelse(!!sym(tax_level) %in% top_taxa, 
                                      !!sym(tax_level), 
                                      "Other")) %>%
    group_by(Sample, !!sym(tax_level)) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
 
  
  # Generate plot
  p <- ggplot(plot_data, aes(x = Sample, y = Abundance, 
                             fill = !!sym(tax_level))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(16)) +
    labs(title = plot_title %||% paste("Top", top_n, tax_level, "Composition"),
         x = "Sample", 
         y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "right")
  
  # Save plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 12, height = 8)
  }
  
  return(list(plot = p, abundance_data = tax_abundance))
}

# Usage example for phylum level ----------------------------------------------
# For single region (V5V7 example)
phy_v5v7 <- phy_regions_filtered_subset[["V5V7"]]

phylum_results <- plot_taxonomic_composition(
  phy_obj = phy_v5v7,
  tax_level = "Phylum",
  top_n = 15,
  palette_name = "Paired",
  plot_title = "Phylum-Level composition in V5V7 Region",
  output_file = "plots/top15_phylum_composition_v5v7.png"
)

# Usage уxample for genus Level -----------------------------------------------
genus_results <- plot_taxonomic_composition(
  phy_obj = phy_v5v7,
  tax_level = "Genus",
  output_file = "plots/top15_genus_composition_v5v7.png"
)

# Batch process all regions --------------------------------------------------
process_all_regions <- function(phy_list, tax_level = "Phylum") {
  map(names(phy_list), ~ {
    region_name <- .x
    phy_obj <- phy_list[[.x]]
    
    if (nsamples(phy_obj) > 0) {
      plot_taxonomic_composition(
        phy_obj = phy_obj,
        tax_level = tax_level,
        output_file = paste0("plots/top15_", tolower(tax_level), 
                             "_composition_", region_name, ".png")
      )
    }
  }) %>% set_names(names(phy_list))
}

#phylum level
all_region_phylum <- process_all_regions(phy_regions_filtered_subset, "Phylum")

# genus level
all_region_genus <- process_all_regions(phy_regions_filtered_subset, "Genus")

