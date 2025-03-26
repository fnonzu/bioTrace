# Load Required Packages ------------------------------------------------------
library(phyloseq)
library(vegan)
library(ggplot2)
library(purrr)

# Configuration ---------------------------------------------------------------
V_REGIONS <- names(phy_regions_filtered_subset)  # Define regions subset or all of them
OUTPUT_DIR <- "plots"
READ_THRESHOLD <- 2000

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# Helper Functions ------------------------------------------------------------
generate_rarefaction_plot <- function(ps_object, region_name, filtered = FALSE) {
  # Prepare OTU matrix
  otu_mat <- as(otu_table(ps_object), "matrix")
  if (taxa_are_rows(ps_object)) otu_mat <- t(otu_mat)
  
  # Skip if no samples
  if (nrow(otu_mat) == 0) return(NULL)
  
  # Generate plot filename
  suffix <- ifelse(filtered, "filtered", "raw")
  file_name <- file.path(OUTPUT_DIR, paste0("rarefaction_", region_name, "_", suffix, ".png"))
  
  # Create plot
  png(file_name, width = 800, height = 600)
  rarecurve(otu_mat, step = 100, sample = min(rowSums(otu_mat)),
            col = rainbow(nrow(otu_mat)), label = FALSE,
            main = paste("Rarefaction сurve", region_name, ifelse(filtered, "(Filtered)", "")),
            xlab = "Number of sequences", ylab = "Observed species")
  
  # Add legend
  legend_labels <- sample_names(ps_object)
  legend("bottomright", legend = legend_labels, col = rainbow(length(legend_labels)),
         lty = 1, bty = "n", cex = 0.8, pt.cex = 0.7)
  
  dev.off()
}

filter_low_read_samples <- function(ps_object, threshold = READ_THRESHOLD) {
  sample_read_sums <- sample_sums(ps_object)
  keep_samples <- names(which(sample_read_sums >= threshold))
  excluded <- setdiff(sample_names(ps_object), keep_samples)
  
  list(
    filtered = prune_samples(keep_samples, ps_object) %>% prune_taxa(taxa_sums(.) > 0, .),
    excluded = excluded
  )
}

# Main analysis ------------------------------------------------------
# Generate raw data rarefaction curves
walk(V_REGIONS, ~ generate_rarefaction_plot(phy_regions_filtered_subset[[.]], .))

# Process filtered data
filter_results <- map(phy_regions_filtered_subset, filter_low_read_samples)

# Generate filtered rarefaction curves
iwalk(filter_results, ~ {
  if (nsamples(.x$filtered) > 0) {
    generate_rarefaction_plot(.x$filtered, .y, filtered = TRUE)
  }
})

# Create final filtered phyloseq object list
phy_final <- map(filter_results, "filtered") %>% compact()

# reporting ---------------------------------------------------------
excluded_samples <- map(filter_results, "excluded")

cat("\n=== Sample summary ===\n")
iwalk(excluded_samples, ~ {
  if (length(.x) > 0) {
    cat(glue::glue("{.y}: Excluded {length(.x)} samples - {paste(.x, collapse = ', ')}\n"))
  } else {
    cat(glue::glue("{.y}: No samples excluded\n"))
  }
})

# Sample Count Comparison -----------------------------------------------------
sample_counts <- tibble(
  Region = V_REGIONS,
  Original = map_dbl(phy_regions_filtered_subset, nsamples),
  Filtered = map_dbl(phy_final, ~ ifelse(is.null(.), 0, nsamples(.)))
)

print("\nSample count comparison:")
print(sample_counts)