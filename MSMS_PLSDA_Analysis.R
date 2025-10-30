# ============================================================================
# PLS-DA Analysis for MS/MS Data (Small Sample Size)
# Adapted for 10 samples in 2 conditions
# ============================================================================

# Load required libraries
# Install if needed: 
# install.packages(c("ggplot2", "ggrepel", "mixOmics", "caret"))
library(ggplot2)
library(ggrepel)
library(mixOmics)  # For PLS-DA
library(caret)     # For cross-validation

# ============================================================================
# 1. LOAD MS/MS DATA
# ============================================================================

# Option A: If you have a feature matrix from xcms or similar
# Columns = samples, Rows = features (m/z or metabolites)
msms_data <- read.csv("feature_matrix.csv", row.names = 1, header = TRUE)

# Option B: If you have a transposed matrix (samples as rows)
# msms_data <- read.csv("feature_matrix.csv", row.names = 1, header = TRUE)
# msms_data <- t(msms_data)  # Transpose if needed

# Load sample annotation
# Should have columns: SampleID, Condition
annotation <- read.csv("sample_annotation.csv", header = TRUE)

# Example of creating annotation if you don't have a file:
# annotation <- data.frame(
#   SampleID = colnames(msms_data),
#   Condition = c(rep("Control", 5), rep("Treatment", 5))
# )

cat("MS/MS data dimensions:\n")
cat("Features:", nrow(msms_data), "\n")
cat("Samples:", ncol(msms_data), "\n\n")

# ============================================================================
# 2. DATA PREPROCESSING
# ============================================================================

# Match annotation to MS/MS data
annotation <- annotation[match(colnames(msms_data), annotation$SampleID), ]

# Remove samples without annotation (if any)
if(any(is.na(annotation$SampleID))) {
  keep_samples <- !is.na(annotation$SampleID)
  msms_data <- msms_data[, keep_samples]
  annotation <- annotation[keep_samples, ]
}

# Ensure Condition is a factor
annotation$Condition <- factor(annotation$Condition)

cat("Sample distribution:\n")
print(table(annotation$Condition))
cat("\n")

# ============================================================================
# 3. DATA FILTERING AND NORMALIZATION
# ============================================================================

# Remove features with too many missing values (>50%)
missing_threshold <- 0.5
missing_prop <- apply(msms_data, 1, function(x) sum(is.na(x))/length(x))
msms_data_filtered <- msms_data[missing_prop < missing_threshold, ]

cat("Features after removing high missing values:", nrow(msms_data_filtered), "\n")

# Impute remaining missing values with half minimum
# (common approach for MS data)
for(i in 1:nrow(msms_data_filtered)) {
  if(any(is.na(msms_data_filtered[i, ]))) {
    min_val <- min(msms_data_filtered[i, ], na.rm = TRUE)
    msms_data_filtered[i, is.na(msms_data_filtered[i, ])] <- min_val / 2
  }
}

# Log2 transformation (typical for MS intensity data)
msms_log <- log2(msms_data_filtered + 1)

# Median normalization (optional but recommended)
median_values <- apply(msms_log, 2, median, na.rm = TRUE)
overall_median <- median(median_values)
msms_normalized <- sweep(msms_log, 2, median_values - overall_median, FUN = "-")

cat("Data normalized\n\n")

# ============================================================================
# 4. FEATURE SELECTION - Keep Most Variable Features
# ============================================================================

# Calculate variance for each feature
feature_variance <- apply(msms_normalized, 1, var, na.rm = TRUE)

# For small datasets, use top 100-500 features
# Adjust based on your total number of features
n_features <- min(500, nrow(msms_normalized))
top_features <- names(sort(feature_variance, decreasing = TRUE)[1:n_features])
msms_filtered <- msms_normalized[top_features, ]

cat("Using top", n_features, "most variable features\n\n")

# Transpose for analysis (samples as rows, features as columns)
msms_t <- t(msms_filtered)

# ============================================================================
# 5. PLS-DA (Partial Least Squares Discriminant Analysis)
# ============================================================================

cat("=== PERFORMING PLS-DA ===\n\n")

# For 2 conditions with 10 samples, use ncomp = 2-3
# (fewer components for small sample size)
ncomp_use <- min(3, length(unique(annotation$Condition)))

# Perform PLS-DA
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = ncomp_use)

# Extract component scores
plsda_scores <- as.data.frame(plsda_result$variates$X)
colnames(plsda_scores) <- paste0("Comp", 1:ncol(plsda_scores))
plsda_scores$Condition <- annotation$Condition
plsda_scores$SampleID <- annotation$SampleID

# Calculate variance explained
plsda_var <- plsda_result$prop_expl_var$X * 100

cat("Variance explained by PLS-DA components:\n")
print(round(plsda_var, 2))
cat("\n")

# ============================================================================
# 6. CROSS-VALIDATION (Important for small sample sizes!)
# ============================================================================

cat("=== PERFORMING CROSS-VALIDATION ===\n\n")

# Leave-One-Out Cross-Validation (LOOCV) - best for small datasets
plsda_cv <- perf(plsda_result, 
                 validation = "loo",  # Leave-one-out
                 progressBar = FALSE)

cat("Classification error rates (LOOCV):\n")
print(plsda_cv$error.rate)
cat("\n")

# Overall error rate
overall_error <- plsda_cv$error.rate$overall[, ncomp_use]
cat("Overall classification error (Component", ncomp_use, "):", 
    round(overall_error, 3), "\n")
cat("Classification accuracy:", round((1 - overall_error) * 100, 1), "%\n\n")

# ============================================================================
# 7. VISUALIZATION - PLS-DA Plot (like PCA)
# ============================================================================

# Define colors for 2 conditions
# Adjust colors as needed
group_colors <- c("#E63946", "#1D3557")  # Red and Blue
names(group_colors) <- levels(annotation$Condition)

# Main PLS-DA plot (Component 1 vs Component 2)
plsda_plot <- ggplot(plsda_scores, aes(x = Comp1, y = Comp2, color = Condition)) +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 1) +  # 95% confidence ellipses
  scale_color_manual(values = group_colors) +
  labs(
    title = "PLS-DA of MS/MS Data",
    subtitle = paste0("Classification Accuracy: ", 
                     round((1 - overall_error) * 100, 1), "%"),
    x = paste0("Component 1 (", round(plsda_var[1], 1), "%)"),
    y = paste0("Component 2 (", round(plsda_var[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(plsda_plot)
ggsave("MSMS_PLSDA_plot.png", plsda_plot, width = 10, height = 8, dpi = 300)
ggsave("MSMS_PLSDA_plot.pdf", plsda_plot, width = 10, height = 8)

# ============================================================================
# 8. PLS-DA WITH SAMPLE LABELS
# ============================================================================

plsda_plot_labeled <- ggplot(plsda_scores, 
                              aes(x = Comp1, y = Comp2, 
                                  color = Condition, label = SampleID)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(size = 3.5, max.overlaps = 20, 
                  box.padding = 0.5, 
                  point.padding = 0.3) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 1) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PLS-DA with Sample Labels",
    subtitle = "MS/MS Data Analysis",
    x = paste0("Component 1 (", round(plsda_var[1], 1), "%)"),
    y = paste0("Component 2 (", round(plsda_var[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(plsda_plot_labeled)
ggsave("MSMS_PLSDA_labeled.png", plsda_plot_labeled, width = 12, height = 9, dpi = 300)

# ============================================================================
# 9. VARIABLE IMPORTANCE IN PROJECTION (VIP) SCORES
# ============================================================================

cat("=== CALCULATING VIP SCORES ===\n\n")

# Extract VIP scores for identifying important features
vip_scores <- vip(plsda_result)

# Get top 50 most important features
n_top <- min(50, nrow(vip_scores))
top_vip_features <- vip_scores[order(vip_scores, decreasing = TRUE)[1:n_top], , drop = FALSE]

# Create dataframe with feature names and VIP scores
vip_df <- data.frame(
  Feature = rownames(top_vip_features),
  VIP_Score = top_vip_features,
  row.names = NULL
)

# Save to CSV
write.csv(vip_df, "MSMS_PLSDA_top_features.csv", row.names = FALSE)

cat("Top 10 discriminant features:\n")
print(head(vip_df, 10))
cat("\n")

# ============================================================================
# 10. VIP SCORE PLOT
# ============================================================================

# Plot top 20 VIP scores
vip_plot_data <- head(vip_df, 20)
vip_plot_data$Feature <- factor(vip_plot_data$Feature, 
                                 levels = vip_plot_data$Feature)

vip_plot <- ggplot(vip_plot_data, aes(x = VIP_Score, y = Feature)) +
  geom_bar(stat = "identity", fill = "#457B9D", alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Top 20 Discriminant Features (VIP Scores)",
    subtitle = "VIP > 1 indicates important features",
    x = "VIP Score",
    y = "Feature (m/z or metabolite)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9)
  )

print(vip_plot)
ggsave("MSMS_VIP_scores.png", vip_plot, width = 10, height = 8, dpi = 300)

# ============================================================================
# 11. LOADING PLOT (Component 1 vs Component 2)
# ============================================================================

# Extract loadings
loadings <- as.data.frame(plsda_result$loadings$X)
colnames(loadings) <- paste0("Comp", 1:ncol(loadings))
loadings$Feature <- rownames(loadings)

# Plot top features in loading space
# Identify features with high absolute loadings
loadings$distance <- sqrt(loadings$Comp1^2 + loadings$Comp2^2)
top_loadings <- loadings[order(loadings$distance, decreasing = TRUE)[1:20], ]

loading_plot <- ggplot(loadings, aes(x = Comp1, y = Comp2)) +
  geom_point(alpha = 0.3, size = 2, color = "grey60") +
  geom_point(data = top_loadings, aes(x = Comp1, y = Comp2), 
             color = "#E63946", size = 3, alpha = 0.8) +
  geom_text_repel(data = top_loadings, 
                  aes(label = Feature), 
                  size = 3, max.overlaps = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  labs(
    title = "PLS-DA Loading Plot",
    subtitle = "Top 20 discriminant features highlighted",
    x = "Loading on Component 1",
    y = "Loading on Component 2"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  )

print(loading_plot)
ggsave("MSMS_loading_plot.png", loading_plot, width = 10, height = 8, dpi = 300)

# ============================================================================
# 12. OPTIONAL: 3-COMPONENT PLOT (if you used ncomp >= 3)
# ============================================================================

if(ncomp_use >= 3) {
  # 3D plot using plotly
  library(plotly)
  
  plot_3d <- plot_ly(plsda_scores, 
                     x = ~Comp1, y = ~Comp2, z = ~Comp3,
                     color = ~Condition, 
                     colors = group_colors,
                     text = ~SampleID,
                     type = "scatter3d", 
                     mode = "markers",
                     marker = list(size = 8)) %>%
    layout(
      title = "3D PLS-DA Plot",
      scene = list(
        xaxis = list(title = paste0("Comp 1 (", round(plsda_var[1], 1), "%)")),
        yaxis = list(title = paste0("Comp 2 (", round(plsda_var[2], 1), "%)")),
        zaxis = list(title = paste0("Comp 3 (", round(plsda_var[3], 1), "%)"))
      )
    )
  
  print(plot_3d)
  
  # Save as HTML
  htmlwidgets::saveWidget(plot_3d, "MSMS_PLSDA_3D.html")
}

# ============================================================================
# 13. SAVE RESULTS SUMMARY
# ============================================================================

# Create results summary
results_summary <- data.frame(
  Metric = c("Total Features", 
             "Features Used", 
             "Number of Samples",
             "Number of Components",
             "Variance Comp1 (%)",
             "Variance Comp2 (%)",
             "Classification Accuracy (%)",
             "Overall Error Rate"),
  Value = c(nrow(msms_data),
            nrow(msms_filtered),
            nrow(msms_t),
            ncomp_use,
            round(plsda_var[1], 2),
            round(plsda_var[2], 2),
            round((1 - overall_error) * 100, 2),
            round(overall_error, 3))
)

write.csv(results_summary, "MSMS_PLSDA_summary.csv", row.names = FALSE)

# Save component scores
write.csv(plsda_scores, "MSMS_PLSDA_scores.csv", row.names = FALSE)

# Save loadings
write.csv(loadings, "MSMS_PLSDA_loadings.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output files saved:\n")
cat("  - MSMS_PLSDA_plot.png/pdf\n")
cat("  - MSMS_PLSDA_labeled.png\n")
cat("  - MSMS_VIP_scores.png\n")
cat("  - MSMS_loading_plot.png\n")
cat("  - MSMS_PLSDA_top_features.csv\n")
cat("  - MSMS_PLSDA_summary.csv\n")
cat("  - MSMS_PLSDA_scores.csv\n")
cat("  - MSMS_PLSDA_loadings.csv\n")

# ============================================================================
# NOTES FOR SMALL SAMPLE SIZES (n=10):
# ============================================================================
# 1. Use Leave-One-Out Cross-Validation (LOOCV) instead of k-fold
# 2. Use fewer components (2-3 maximum)
# 3. Feature selection is critical - use only most variable features
# 4. Consider more conservative significance thresholds
# 5. Report cross-validation accuracy along with the plot
# 6. Be cautious about overfitting - validate on independent data if possible
# ============================================================================
