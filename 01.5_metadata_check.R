library(rio) 
library(dplyr)
#############################
##### DATA IMPORT ###########
#############################
# reading data from all sheets 
trace_metals_data <- import_list("~/p952/20240918_Habkern_compiled_soil_data_v2.xlsx") 

trace_metals <- trace_metals_data[["Trace metals"]]

colnames(trace_metals) <- trace_metals[1, ]
trace_metals <- trace_metals[-2, ]
trace_metals <- trace_metals[8:25, ]
trace_metals <- trace_metals[ , -2:-3]

trace_metals <- trace_metals %>%
  rename(
    Pb = '208  Pb  [ He ]',
    Cu = '63  Cu  [ He ]',
    Sb_121 = '121  Sb  [ He ]',
    Sb_123 = '123  Sb  [ He ]'
  )
trace_metals <- as.data.frame(trace_metals)

#average Sb
trace_metals <- trace_metals %>%
  mutate(
    Sb_121 = as.numeric(as.character(Sb_121)),
    Sb_123 = as.numeric(as.character(Sb_123))
  ) %>%
  dplyr::mutate(Sb = (Sb_121 + Sb_123)/2) %>%
  dplyr::select(-Sb_121, -Sb_123)

trace_metals <- trace_metals %>%
  rename(
    '208  Pb  [ He ]' = Pb,
    '63  Cu  [ He ]' = Cu,
    '121-123  Sb  [ He ]' = Sb,
  )

rownames(trace_metals_metadata) <- NULL

############################################ NOT USED, maybe some parts have useful transformation
library(dplyr)      
library(tidyr)      
library(ggplot2)    
library(vegan)      
library(ggfortify)  
### Step 1: Exploratory Data Analysis (EDA) with Normalization ###

# Create a new column 'SampleID'
trace_metals <- trace_metals %>%
  mutate(SampleID = sub("\\..*", "", Sample))

# Convert metal columns to numeric (columns 2 to ncol-1, excluding the newly added SampleID)
trace_metals[, 2:(ncol(trace_metals)-1)] <- lapply(trace_metals[, 2:(ncol(trace_metals)-1)], as.numeric)

# Aggregate the data by SampleID (average the trace metal values)
agg_metals <- trace_metals %>%
  group_by(SampleID) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

# Normalize the metal concentrations (standardization: center and scale)
# We exclude the 'SampleID' column when normalizing
agg_metals_norm <- agg_metals
agg_metals_norm[,-1] <- scale(agg_metals_norm[,-1])

# Inspect the normalized aggregated data
print(head(agg_metals_norm))
summary(agg_metals_norm)

# Calculate and view the correlation matrix among the normalized metals
# (Excluding the SampleID column)
cor_matrix <- cor(agg_metals_norm[,-1], use = "pairwise.complete.obs")
print(cor_matrix)

# Reshape normalized data for visualization: histogram and boxplots of each metal
metals_long <- pivot_longer(agg_metals_norm, cols = -SampleID, names_to = "Metal", values_to = "Value")

# Histogram for each trace metal
ggplot(metals_long, aes(x = Value)) + 
  geom_histogram(bins = 20, fill = "skyblue", color = "black") + 
  facet_wrap(~Metal, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Normalized Trace Metals", 
       x = "Standardized Concentration", y = "Frequency")

# Boxplot for each trace metal
ggplot(metals_long, aes(x = Metal, y = Value)) +
  geom_boxplot(fill = "lightgreen") +
  theme_minimal() +
  labs(title = "Boxplot of Normalized Trace Metals", 
       x = "Metal", y = "Standardized Concentration")


### Step 2: Principal Component Analysis (PCA) ###

# Run PCA on the aggregated trace metal data
# Exclude the 'SampleID' column and scale the data
pca_res <- prcomp(agg_metals[,-1:-2], scale. = TRUE)

# View a summary of the PCA results (variance explained, loadings, etc.)
summary(pca_res)

# Base R scree plot to inspect variance explained by each component
screeplot(pca_res, type = "lines", main = "Scree Plot of PCA")

# Base R biplot to visualize samples and metal loadings
biplot(pca_res, scale = 0)

# Alternatively, use ggfortify for a nicer PCA plot with sample labels
autoplot(pca_res, data = agg_metals, label = TRUE, loadings = TRUE, loadings.label = TRUE, 
         main = "PCA Biplot of Trace Metals")



############## NORMILIZE TRACE METALS - USED

## transform metadata for z-scores to fit it in and also make it numeric
trace_metals_metadata_scaled <- trace_metals_metadata %>%
  mutate(across(c(Pb, Cu, Sb), as.numeric)) %>% 
  mutate(across(c(Pb, Cu, Sb), ~ as.numeric(scale(.x))))

# View the standardized data
head(trace_metals_metadata_scaled)

