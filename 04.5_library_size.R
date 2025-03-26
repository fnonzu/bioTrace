library(phyloseq)
library(tidyverse)
library(RColorBrewer)

# Extract and transform library size data -------------------------------------
library_sizes <- phy_regions_filtered_subset %>% 
  imap_dfr(~ {
    tibble(
      Sample = sample_names(.x),
      LibrarySize = sample_sums(.x),
      Region = .y
    )
  }) %>%
  mutate(Region = fct_inorder(Region))  # Preserve original region order

# Calculate summary statistics ------------------------------------------------
library_summary <- library_sizes %>%
  group_by(Region) %>%
  summarise(across(LibrarySize,
                   list(
                     Mean = ~ mean(.x) %>% round() %>% format(big.mark = ","),
                     Median = ~ median(.x) %>% round() %>% format(big.mark = ","),
                     Min = ~ min(.x) %>% format(big.mark = ","),
                     Max = ~ max(.x) %>% format(big.mark = ",")
                   ),
                   .names = "{.fn}"
  ))

print(library_summary)

# Create visualization --------------------------------------------------------
# Set color palette
region_palette <- brewer.pal(n = length(phy_regions_filtered_subset), name = "Set2")

library_boxplot <- library_sizes %>%
  ggplot(aes(x = Region, y = LibrarySize, fill = Region)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  scale_fill_manual(values = region_palette) +
  scale_y_continuous(
    labels = scales::comma,
    trans = "log10",
    breaks = scales::trans_breaks("log10", function(x) 10^x)
  ) +
  labs(
    title = "Sequencing counts distribution by amplification region",
    x = "Hypervariable region",
    y = "Counts number (log10 scale)",
    caption = "Points represent individual samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

# Display plot
print(library_boxplot)