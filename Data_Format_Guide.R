# ============================================================================
# DATA FORMAT EXAMPLES FOR MS/MS PLS-DA ANALYSIS
# ============================================================================

# EXAMPLE 1: Feature Matrix Format (most common from xcms, MZmine, etc.)
# -----------------------------------------------------------------------
# CSV file: feature_matrix.csv
# 
# Structure: Rows = features (m/z or metabolites), Columns = samples
#
#              ,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6,Sample7,Sample8,Sample9,Sample10
# Feature_001  ,1234.5 ,1456.2 ,1389.7 ,1298.4 ,1367.9 ,2456.1 ,2389.4 ,2512.3 ,2401.8 ,2478.6
# Feature_002  ,567.8  ,601.3  ,589.2  ,578.9  ,595.1  ,423.4  ,445.7  ,438.9  ,429.2  ,441.5
# Feature_003  ,8901.2 ,8756.4 ,8823.9 ,8912.7 ,8834.5 ,9123.8 ,9234.1 ,9178.6 ,9201.3 ,9189.7
# ...

# EXAMPLE 2: Sample Annotation Format
# ------------------------------------
# CSV file: sample_annotation.csv
#
# SampleID,Condition,Batch,InjectionOrder
# Sample1,Control,1,1
# Sample2,Control,1,2
# Sample3,Control,1,3
# Sample4,Control,1,4
# Sample5,Control,1,5
# Sample6,Treatment,1,6
# Sample7,Treatment,1,7
# Sample8,Treatment,1,8
# Sample9,Treatment,1,9
# Sample10,Treatment,1,10

# ============================================================================
# HOW TO PREPARE YOUR DATA FROM DIFFERENT SOURCES
# ============================================================================

# FROM XCMS (R package)
# ---------------------
library(xcms)

# After running xcms peak detection and alignment:
# xdata <- your xcms object after fillChromPeaks()

# Extract feature matrix
feature_matrix <- featureValues(xdata, value = "into")  # "into" = integrated intensity
# Or use "maxo" for maximum intensity

# Feature matrix is already in correct format (features × samples)
write.csv(feature_matrix, "feature_matrix.csv")

# Create annotation from xcms sample data
sample_info <- pData(xdata)
annotation <- data.frame(
  SampleID = rownames(sample_info),
  Condition = sample_info$class  # or whatever your grouping column is
)
write.csv(annotation, "sample_annotation.csv", row.names = FALSE)


# FROM MZmine
# -----------
# Export your feature table from MZmine as CSV
# Make sure to:
# 1. Include sample columns
# 2. Export with feature ID/m/z as first column
# 
# Then in R:
mzmine_data <- read.csv("mzmine_export.csv", header = TRUE)

# Extract only the peak area/height columns for your samples
# Assuming columns 10-19 are your sample intensity values
feature_matrix <- mzmine_data[, 10:19]
rownames(feature_matrix) <- mzmine_data$ID  # or use m/z values

write.csv(feature_matrix, "feature_matrix.csv")


# FROM MS-DIAL
# ------------
# Export alignment result
# Select your samples and export as text/CSV
# Similar process to MZmine

msdial_data <- read.delim("msdial_alignment.txt", header = TRUE)

# Extract intensity columns (adjust column numbers based on your export)
sample_cols <- grep("Sample.*Height", colnames(msdial_data))  # or "Area"
feature_matrix <- msdial_data[, sample_cols]
rownames(feature_matrix) <- paste0("mz_", msdial_data$Average.Mz)

write.csv(feature_matrix, "feature_matrix.csv")


# FROM METABOLANALYST FORMAT
# --------------------------
# If your data is in MetaboAnalyst format (samples as rows):
metabolanalyst_data <- read.csv("metabolanalyst_data.csv", row.names = 1)

# Transpose it
feature_matrix <- t(metabolanalyst_data)

write.csv(feature_matrix, "feature_matrix.csv")


# ============================================================================
# CREATING ANNOTATION FILE
# ============================================================================

# Method 1: Manually create annotation
# -------------------------------------
sample_names <- colnames(feature_matrix)

annotation <- data.frame(
  SampleID = sample_names,
  Condition = c(rep("Control", 5), rep("Treatment", 5)),  # Adjust to your design
  Batch = rep(1, 10),  # If you have batch information
  InjectionOrder = 1:10  # Sequence order
)

write.csv(annotation, "sample_annotation.csv", row.names = FALSE)


# Method 2: Create from sample naming pattern
# --------------------------------------------
# If your samples are named like: "Ctrl_1", "Ctrl_2", ..., "Treat_1", etc.
sample_names <- colnames(feature_matrix)

# Extract condition from sample names
conditions <- ifelse(grepl("Ctrl|Control", sample_names), "Control", "Treatment")

annotation <- data.frame(
  SampleID = sample_names,
  Condition = conditions
)

write.csv(annotation, "sample_annotation.csv", row.names = FALSE)


# ============================================================================
# QUALITY CHECK YOUR DATA BEFORE PLS-DA
# ============================================================================

# Load your prepared data
feature_matrix <- read.csv("feature_matrix.csv", row.names = 1)
annotation <- read.csv("sample_annotation.csv")

# Check dimensions
cat("Features:", nrow(feature_matrix), "\n")
cat("Samples:", ncol(feature_matrix), "\n")
cat("\nSample names match:\n")
print(all(colnames(feature_matrix) == annotation$SampleID))

# Check for missing values
missing_prop <- sum(is.na(feature_matrix)) / (nrow(feature_matrix) * ncol(feature_matrix))
cat("\nProportion of missing values:", round(missing_prop * 100, 2), "%\n")

# Check data distribution
hist(log2(as.matrix(feature_matrix) + 1), 
     main = "Distribution of Feature Intensities (log2)",
     xlab = "log2(Intensity + 1)",
     breaks = 50)

# Check for outliers with simple PCA
library(ggplot2)
pca_check <- prcomp(t(log2(feature_matrix + 1)), scale. = TRUE)
pca_df <- as.data.frame(pca_check$x[, 1:2])
pca_df$Condition <- annotation$Condition

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = rownames(pca_df))) +
  geom_point(size = 3) +
  geom_text_repel() +
  theme_bw() +
  labs(title = "Quick PCA Check")

# If everything looks good, proceed with PLS-DA analysis!

# ============================================================================
# COMPLETE WORKFLOW EXAMPLE
# ============================================================================

# 1. Prepare data (choose one method above)
feature_matrix <- read.csv("feature_matrix.csv", row.names = 1)
annotation <- read.csv("sample_annotation.csv")

# 2. Run the PLS-DA analysis
source("MSMS_PLSDA_Analysis.R")  # Use the comprehensive script
# OR
source("MSMS_PLSDA_QuickStart.R")  # Use the quick start script

# 3. Interpret results:
#    - Check the PLS-DA plot for separation
#    - Look at cross-validation accuracy
#    - Examine VIP scores for important features
#    - Investigate top discriminant features in your data
