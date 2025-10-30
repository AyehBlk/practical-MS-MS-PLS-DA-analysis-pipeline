# ============================================================================
# QUICK START: PLS-DA for MS/MS Data (10 samples, 2 conditions)
# ============================================================================

library(ggplot2)
library(ggrepel)
library(mixOmics)

# 1. LOAD YOUR DATA
# -----------------
# Your feature matrix: rows = features (m/z), columns = samples
msms_data <- read.csv("your_feature_matrix.csv", row.names = 1, header = TRUE)

# Your sample annotation
annotation <- data.frame(
  SampleID = colnames(msms_data),
  Condition = c(rep("Control", 5), rep("Treatment", 5))  # Adjust to your groups
)

# 2. PREPROCESSING
# ----------------
# Log transform and normalize
msms_log <- log2(msms_data + 1)
msms_t <- t(msms_log)  # Transpose: samples as rows

# 3. PLS-DA
# ---------
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = 2)

# Get scores
plsda_scores <- as.data.frame(plsda_result$variates$X)
colnames(plsda_scores) <- c("Comp1", "Comp2")
plsda_scores$Condition <- annotation$Condition
plsda_scores$SampleID <- annotation$SampleID

# Variance explained
plsda_var <- plsda_result$prop_expl_var$X * 100

# 4. CROSS-VALIDATION
# -------------------
plsda_cv <- perf(plsda_result, validation = "loo", progressBar = FALSE)
accuracy <- (1 - plsda_cv$error.rate$overall[, 2]) * 100

# 5. PLOT
# -------
group_colors <- c("Control" = "#E63946", "Treatment" = "#1D3557")

plsda_plot <- ggplot(plsda_scores, aes(x = Comp1, y = Comp2, color = Condition)) +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 1) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PLS-DA of MS/MS Data",
    subtitle = paste0("Accuracy: ", round(accuracy, 1), "%"),
    x = paste0("Component 1 (", round(plsda_var[1], 1), "%)"),
    y = paste0("Component 2 (", round(plsda_var[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right"
  )

print(plsda_plot)
ggsave("PLSDA_plot.png", plsda_plot, width = 10, height = 8, dpi = 300)

# 6. TOP DISCRIMINANT FEATURES
# -----------------------------
vip_scores <- vip(plsda_result)
top_features <- head(sort(vip_scores, decreasing = TRUE), 20)
print(top_features)
