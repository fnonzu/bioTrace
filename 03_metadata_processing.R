# Load Required Packages ------------------------------------------------------
library(tidyverse)
library(rio)
library(ggfortify)
library(vegan)
library(gridExtra)

# Configuration ---------------------------------------------------------------
DATA_PATH <- "~/p952/20240918_Habkern_compiled_soil_data_v2.xlsx"
METAL_COLS <- c('208  Pb  [ He ]', '63  Cu  [ He ]', 
                '121  Sb  [ He ]', '123  Sb  [ He ]')

# Data Import and Cleaning ----------------------------------------------------
import_trace_metals <- function() {
  # Import and initial processing
  trace_metals_data <- import_list(DATA_PATH)
  raw_metals <- trace_metals_data[["Trace metals"]]
  
  # Set column names and clean data
  raw_metals %>%
    set_names(slice(., 1)) %>%
    slice(-1:-2) %>%              # Remove header rows
    slice(8:25) %>%               # Select relevant sample rows
    select(-2:-3) %>%             # Remove unnecessary columns
    rename_with(~str_remove_all(., "\\s+"), everything()) %>%  # Simplify names
    rename(
      Pb = `208Pb[He]`,
      Cu = `63Cu[He]`,
      Sb_121 = `121Sb[He]`,
      Sb_123 = `123Sb[He]`
    ) %>%
    mutate(across(c(Pb, Cu, Sb_121, Sb_123), as.numeric)) %>%
    mutate(Sb = (Sb_121 + Sb_123) / 2) %>%
    select(-Sb_121, -Sb_123)
}

# Load and process data
trace_metals <- import_trace_metals()

# Exploratory Data Analysis ---------------------------------------------------
prepare_analysis_data <- function(metal_df) {
  metal_df %>%
    mutate(SampleID = str_extract(Sample, "^[^.]+")) %>%  # Create SampleID
    mutate(across(where(is.numeric), ~ scale(.x)[,1]))    # Z-score normalization
}

# Create normalized dataset
metals_normalized <- prepare_analysis_data(trace_metals)

# Visualization Functions -----------------------------------------------------
plot_metal_distributions <- function(norm_df) {
  norm_df %>%
    pivot_longer(-c(Sample, SampleID), 
                 names_to = "metal", values_to = "value") %>%
    ggplot(aes(value)) +
    geom_histogram(bins = 20, fill = "skyblue", color = "black") +
    geom_density(aes(y = after_stat(count)), color = "red") +
    facet_wrap(~metal, scales = "free") +
    labs(title = "Distribution of Normalized Metal Concentrations",
         x = "Standardized Value", y = "Frequency") +
    theme_minimal(base_size = 12)
}

plot_metal_correlations <- function(norm_df) {
  cor_matrix <- norm_df %>%
    select(where(is.numeric)) %>%
    cor(use = "complete.obs")
  
  corrplot::corrplot(cor_matrix, method = "color", type = "upper",
                     tl.col = "black", tl.srt = 45,
                     addCoef.col = "black", number.cex = 0.7,
                     title = "Trace Metal Correlation Matrix")
}

# Principal Component Analysis ------------------------------------------------
perform_pca_analysis <- function(metal_df) {
  pca_data <- metal_df %>%
    select(where(is.numeric)) %>%
    drop_na()
  
  pca_res <- prcomp(pca_data, scale. = TRUE)
  
  list(
    variance = summary(pca_res),
    biplot = autoplot(pca_res, data = metal_df, label = TRUE, 
                      loadings = TRUE, loadings.label = TRUE,
                      main = "PCA of Trace Metal Composition")
  )
}

# Execute Analysis ------------------------------------------------------------
# Generate visualizations
dist_plot <- plot_metal_distributions(metals_normalized)
corr_plot <- plot_metal_correlations(metals_normalized)

# Perform PCA
pca_results <- perform_pca_analysis(metals_normalized)

# Display results
print(dist_plot)
print(corr_plot)
print(pca_results$biplot)



############################ OTHER METADATA ############################
## PH ##
pH <- trace_metals_data[["pH"]]
names(pH)[3] <- "pH in CaCl2"
pH[["pH in CaCl2"]] <- as.numeric(pH[["pH in CaCl2"]])

average_pH <- mean(pH[1:3, "pH in CaCl2"], na.rm = TRUE)
pH[1:3, "pH in CaCl2"] <- average_pH
pH <- pH[-c(2, 3), ]
average_pH <- mean(pH[6:8, "pH in CaCl2"], na.rm = TRUE)
pH[6:8, "pH in CaCl2"] <- average_pH
pH <- pH[-c(7, 8), ]
average_pH <- mean(pH[12:14, "pH in CaCl2"], na.rm = TRUE)
pH[12:14, "pH in CaCl2"] <- average_pH
pH <- pH[-c(13, 14), ]
rownames(pH) <- NULL

## ORGANIC MATTER 
org_mat <- trace_metals_data[["Organic Matter %"]]

## GRAIN SIZE
grain <- trace_metals_data[["grain size"]]
colnames(grain)[1] <- "Sample Name"
grain <- grain[-1, ]
average_grain <- grain %>%
  filter(grepl("- Average$", `Sample Name`))  # Keep only averages

average_grain <- average_grain %>%
  mutate(
    `clay fraction (%)` = as.numeric(as.character(`clay fraction (%)`)),
    `silt fraction (%)` = as.numeric(as.character(`silt fraction (%)`)),
    `sand fraction (%)` = as.numeric(as.character(`sand fraction (%)`))
  )
average_grain <- average_grain %>% 
  mutate(sample_group = sub("^([0-9]+).*", "\\1", `Sample Name`))
average_grain <- average_grain %>%
  group_by(sample_group) %>%
  summarise(
    clay_fraction = mean(`clay fraction (%)`, na.rm = TRUE),
    silt_fraction = mean(`silt fraction (%)`, na.rm = TRUE),
    sand_fraction = mean(`sand fraction (%)`, na.rm = TRUE)
  )
average_grain$sample_group <- as.numeric(as.character(average_grain$sample_group))
rownames(average_grain) <- average_grain$sample_group # WEIRD



# Add sample numbers from row names to each data frame
pH$Sample <- as.numeric(rownames(pH))
org_mat$Sample <- as.numeric(rownames(org_mat))

# Standardize grain size data to sum to 100% per row
average_grain[, 2:4] <- average_grain[, 2:4] / rowSums(average_grain[, 2:4]) * 100

# Prepare grain data for stacked bars
grain_long <- average_grain %>% 
  mutate(Sample = as.numeric(average_grain$sample_group)) %>% 
  pivot_longer(cols = 2:4, names_to = "Grain_Size", values_to = "Percentage")

# Open a PNG device to save the plot
png("plots/metadata_plot_3inrow.png", width = 1500, height = 700)
# Create pH plot
p_pH <- ggplot(pH, aes(x = Sample, y = pH[,3])) +
  geom_line(color = "blue", linewidth = 0.8) +
  geom_point(color = "blue", size = 2) +
  labs(title = "pH in CaCl2", x = "Sample", y = "pH value") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal()

# Create organic matter plot
p_org <- ggplot(org_mat, aes(x = Sample, y = org_mat[,3])) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_point(color = "darkgreen", size = 2) +
  labs(title = "Organic matter", x = "Sample", y = "Organic matter percentage %") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal()

# Create stacked grain size plot
p_grain <- ggplot(grain_long, aes(x = factor(Sample), y = Percentage, fill = Grain_Size)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Grain composition", x = "Sample", y = "Percentage (%)") +
  scale_x_discrete(breaks = 1:12) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine all plots in a grid
grid.arrange(p_pH, p_org, p_grain, 
             ncol = 3,
             widths = c(1, 1, 1.2))  # Adjust relative widths
dev.off()
# Add sample numbers from row names to each data frame
pH$Sample <- as.numeric(rownames(pH))
org_mat$Sample <- as.numeric(rownames(org_mat))

# Standardize grain size data to sum to 100% per row
average_grain[, 2:4] <- average_grain[, 2:4] / rowSums(average_grain[, 2:4]) * 100

# Prepare grain data for stacked bars
grain_long <- average_grain %>% 
  mutate(Sample = as.numeric(average_grain$sample_group)) %>% 
  pivot_longer(cols = 2:4, names_to = "Grain_Size", values_to = "Percentage")

# Open a PNG device to save the plot
png("plots/metadata_plot.png", width = 1500, height = 900)  # Increased height

# Create pH plot
p_pH <- ggplot(pH, aes(x = Sample, y = pH[,3])) +
  geom_line(color = "blue", linewidth = 0.8) +
  geom_point(color = "blue", size = 2) +
  labs(title = "pH in CaCl₂", x = "Sample", y = "pH value") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 5, "pt"))

# Create organic matter plot
p_org <- ggplot(org_mat, aes(x = Sample, y = org_mat[,3])) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_point(color = "darkgreen", size = 2) +
  labs(title = "Organic matter", x = "Sample", y = "Organic matter (%)") +
  scale_x_continuous(breaks = 1:12, limits = c(1, 12)) +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 5, "pt"))

# Create stacked grain size plot
p_grain <- ggplot(grain_long, aes(x = factor(Sample), y = Percentage, fill = Grain_Size)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Grain composition", x = "Sample", y = "Percentage (%)") +
  scale_x_discrete(breaks = 1:12) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.margin = margin(5, 15, 5, 15, "pt"))

# Combine plots in new layout
grid.arrange(
  arrangeGrob(p_pH, p_org, ncol = 1, heights = c(1, 1)),  # Left column
  p_grain,  # Right column
  ncol = 2, 
  widths = c(1, 1.5)
)

dev.off()
