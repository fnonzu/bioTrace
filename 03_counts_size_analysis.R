# Extract library sizes for each region
library_sizes <- lapply(names(phy_regions_filtered_subset), function(region) {
  data.frame(
    Sample = sample_names(phy_regions_filtered[[region]]),
    LibrarySize = sample_sums(phy_regions_filtered[[region]]),
    Region = region
  )
}) %>%
  bind_rows()

# Summary statistics
library_summary <- library_sizes %>%
  group_by(Region) %>%
  summarise(
    MeanSize = mean(LibrarySize),
    MedianSize = median(LibrarySize),
    MinSize = min(LibrarySize),
    MaxSize = max(LibrarySize)
  )
print(library_summary)

# Boxplot for library sizes by region
library_boxplot <- ggplot(library_sizes, aes(x = Region, y = LibrarySize, fill = Region)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("V Region") + ylab("Counts") +
  labs(title = "Number of Counts by V Region") +
  theme(legend.position = "none")
print(library_boxplot)


library_sample_bar <- ggplot(library_sizes, aes(x = Sample, y = LibrarySize, fill = Region)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Sample") +
  ylab("Counts") +
  labs(title = "Number of Counts by Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(library_sample_bar)
