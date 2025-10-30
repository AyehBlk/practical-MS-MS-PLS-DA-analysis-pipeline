# MS/MS Data Analysis with PLS-DA in R

## 📋 Table of Contents
- [Overview](#overview)
- [What is PLS-DA?](#what-is-pls-da)
- [Why Use PLS-DA for MS/MS Data?](#why-use-pls-da-for-msms-data)
- [Requirements](#requirements)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Quick Start Guide](#quick-start-guide)
- [Detailed Workflow](#detailed-workflow)
- [Understanding the Output](#understanding-the-output)
- [Interpreting Results](#interpreting-results)
- [Important Considerations for Small Sample Sizes](#important-considerations-for-small-sample-sizes)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [References](#references)

---

## Overview

This repository contains R scripts for performing **Partial Least Squares Discriminant Analysis (PLS-DA)** on MS/MS (tandem mass spectrometry) data. The analysis is specifically optimized for small sample sizes (e.g., 10 samples in 2 conditions) and produces publication-quality visualizations similar to PCA plots but with supervised classification.

### 📁 Files Included

1. **MSMS_PLSDA_Analysis.R** - Comprehensive analysis pipeline (recommended)
2. **MSMS_PLSDA_QuickStart.R** - Minimal quick-start script
3. **Data_Format_Guide.R** - Data formatting examples from various MS platforms
4. **README.md** - This file

---

## What is PLS-DA?

### Basic Concept

**PLS-DA (Partial Least Squares Discriminant Analysis)** is a supervised multivariate statistical method that:

1. **Finds patterns** in high-dimensional data (many features, few samples)
2. **Maximizes separation** between predefined groups (e.g., Control vs Treatment)
3. **Identifies important features** that discriminate between groups
4. **Creates low-dimensional visualizations** (like PCA, but supervised)

### PLS-DA vs PCA

| Feature | PCA | PLS-DA |
|---------|-----|--------|
| **Type** | Unsupervised | Supervised |
| **Goal** | Maximize variance | Maximize group separation |
| **Uses group info** | No | Yes |
| **Better for classification** | No | Yes |
| **Risk of overfitting** | Lower | Higher (needs validation) |

### When to Use Each

- **Use PCA** when:
  - Exploring data without prior knowledge
  - Detecting outliers
  - Checking batch effects
  - No predefined groups

- **Use PLS-DA** when:
  - You have predefined groups
  - You want to classify samples
  - You need to identify discriminant features
  - You want to maximize group separation

---

## Why Use PLS-DA for MS/MS Data?

### MS/MS Data Characteristics

Mass spectrometry data has unique properties that make PLS-DA particularly suitable:

1. **High dimensionality**: Thousands of features (m/z values, metabolites)
2. **Small sample size**: Often limited to 10-30 samples due to cost
3. **Noise and missing values**: Instrument variability and detection limits
4. **Collinear features**: Related metabolites or isotopes
5. **Need for biomarker discovery**: Identifying discriminant features

### PLS-DA Advantages for MS/MS

✅ **Handles high-dimensional data** - Works well with features >> samples  
✅ **Robust to noise** - Focuses on discriminant patterns  
✅ **Feature selection** - VIP scores identify important biomarkers  
✅ **Visual interpretation** - Easy-to-understand plots  
✅ **Classification performance** - Can predict sample class  

---

## Requirements

### Software
- **R** (version ≥ 4.0.0)
- **RStudio** (recommended but optional)

### R Packages

```r
# Core packages
install.packages(c("ggplot2", "ggrepel", "caret"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")

# Optional packages
install.packages(c("plotly", "pheatmap", "htmlwidgets"))
```

### System Requirements
- **RAM**: Minimum 4 GB (8 GB recommended)
- **Storage**: ~500 MB for packages and data
- **OS**: Windows, macOS, or Linux

---

## Installation

### Step 1: Install R and RStudio

1. Download R from [CRAN](https://cran.r-project.org/)
2. Download RStudio from [Posit](https://posit.co/download/rstudio-desktop/)

### Step 2: Install Required Packages

Open R or RStudio and run:

```r
# Install CRAN packages
install.packages(c("ggplot2", "ggrepel", "caret"))

# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install mixOmics
BiocManager::install("mixOmics")

# Verify installation
library(mixOmics)
library(ggplot2)
library(ggrepel)
```

### Step 3: Download Scripts

Download all `.R` files to a working directory on your computer.

---

## Data Preparation

### Required Data Files

You need **two files** to run the analysis:

#### 1. Feature Matrix (`feature_matrix.csv`)

**Format**: Rows = features (m/z or metabolites), Columns = samples

```
              ,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6,Sample7,Sample8,Sample9,Sample10
Feature_001   ,1234.5 ,1456.2 ,1389.7 ,1298.4 ,1367.9 ,2456.1 ,2389.4 ,2512.3 ,2401.8 ,2478.6
Feature_002   ,567.8  ,601.3  ,589.2  ,578.9  ,595.1  ,423.4  ,445.7  ,438.9  ,429.2  ,441.5
Feature_003   ,8901.2 ,8756.4 ,8823.9 ,8912.7 ,8834.5 ,9123.8 ,9234.1 ,9178.6 ,9201.3 ,9189.7
```

**Key points**:
- First column: Feature IDs (e.g., m/z values, metabolite names)
- Each subsequent column: One sample
- Values: Peak intensities/areas/heights
- Missing values: Use `NA` or leave empty

#### 2. Sample Annotation (`sample_annotation.csv`)

**Format**: Each row = one sample

```
SampleID,Condition,Batch,InjectionOrder
Sample1,Control,1,1
Sample2,Control,1,2
Sample3,Control,1,3
Sample4,Control,1,4
Sample5,Control,1,5
Sample6,Treatment,1,6
Sample7,Treatment,1,7
Sample8,Treatment,1,8
Sample9,Treatment,1,9
Sample10,Treatment,1,10
```

**Required columns**:
- `SampleID`: Must match column names in feature matrix
- `Condition`: Group labels (e.g., Control, Treatment, Diseased, Healthy)

**Optional columns**:
- `Batch`: Batch information
- `InjectionOrder`: Run sequence
- Any other metadata

### Data Sources

See `Data_Format_Guide.R` for detailed instructions on preparing data from:

- **xcms** (R package)
- **MZmine** (Java-based software)
- **MS-DIAL** (Windows software)
- **MetaboAnalyst** (web-based tool)
- **Compound Discoverer** (Thermo)
- **Progenesis QI** (Waters)
- **Custom pipelines**

---

## Quick Start Guide

### Minimal Working Example

```r
# Load libraries
library(ggplot2)
library(ggrepel)
library(mixOmics)

# 1. LOAD DATA
msms_data <- read.csv("feature_matrix.csv", row.names = 1)
annotation <- read.csv("sample_annotation.csv")

# 2. PREPROCESS
msms_log <- log2(msms_data + 1)
msms_t <- t(msms_log)  # Transpose: samples as rows

# 3. PLS-DA
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = 2)

# Get scores
plsda_scores <- as.data.frame(plsda_result$variates$X)
colnames(plsda_scores) <- c("Comp1", "Comp2")
plsda_scores$Condition <- annotation$Condition

# Variance explained
plsda_var <- plsda_result$prop_expl_var$X * 100

# 4. CROSS-VALIDATION
plsda_cv <- perf(plsda_result, validation = "loo", progressBar = FALSE)
accuracy <- (1 - plsda_cv$error.rate$overall[, 2]) * 100

# 5. PLOT
ggplot(plsda_scores, aes(x = Comp1, y = Comp2, color = Condition)) +
  geom_point(size = 5) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PLS-DA of MS/MS Data",
    subtitle = paste0("Accuracy: ", round(accuracy, 1), "%"),
    x = paste0("Component 1 (", round(plsda_var[1], 1), "%)"),
    y = paste0("Component 2 (", round(plsda_var[2], 1), "%)")
  ) +
  theme_bw()

# 6. TOP FEATURES
vip_scores <- vip(plsda_result)
top_features <- head(sort(vip_scores, decreasing = TRUE), 20)
print(top_features)
```

**Run time**: ~30 seconds for 10 samples × 1000 features

---

## Detailed Workflow

### Step-by-Step Analysis Pipeline

The complete analysis consists of 13 steps:

#### 1. Load Data
```r
msms_data <- read.csv("feature_matrix.csv", row.names = 1, header = TRUE)
annotation <- read.csv("sample_annotation.csv", header = TRUE)
```

#### 2. Data Preprocessing

**Check dimensions**:
```r
cat("Features:", nrow(msms_data), "\n")
cat("Samples:", ncol(msms_data), "\n")
```

**Match annotation to data**:
```r
annotation <- annotation[match(colnames(msms_data), annotation$SampleID), ]
```

**Ensure factor levels**:
```r
annotation$Condition <- factor(annotation$Condition)
```

#### 3. Data Filtering

**Remove features with excessive missing values**:
```r
missing_threshold <- 0.5  # Remove if >50% missing
missing_prop <- apply(msms_data, 1, function(x) sum(is.na(x))/length(x))
msms_filtered <- msms_data[missing_prop < missing_threshold, ]
```

**Impute remaining missing values**:
```r
# Method: Half-minimum (conservative approach for MS data)
for(i in 1:nrow(msms_filtered)) {
  if(any(is.na(msms_filtered[i, ]))) {
    min_val <- min(msms_filtered[i, ], na.rm = TRUE)
    msms_filtered[i, is.na(msms_filtered[i, ])] <- min_val / 2
  }
}
```

#### 4. Data Transformation

**Log2 transformation** (standard for MS intensity data):
```r
msms_log <- log2(msms_filtered + 1)
```

**Why log transform?**
- Reduces heteroscedasticity
- Makes data more normally distributed
- Reduces influence of highly abundant features

#### 5. Normalization

**Median normalization**:
```r
median_values <- apply(msms_log, 2, median, na.rm = TRUE)
overall_median <- median(median_values)
msms_normalized <- sweep(msms_log, 2, median_values - overall_median, FUN = "-")
```

**Why normalize?**
- Corrects for technical variation
- Makes samples comparable
- Reduces batch effects

#### 6. Feature Selection

**Keep most variable features**:
```r
feature_variance <- apply(msms_normalized, 1, var, na.rm = TRUE)
n_features <- min(500, nrow(msms_normalized))  # Use top 500 features
top_features <- names(sort(feature_variance, decreasing = TRUE)[1:n_features])
msms_filtered <- msms_normalized[top_features, ]
```

**Why feature selection?**
- Reduces noise
- Prevents overfitting
- Focuses on informative features
- Speeds up computation

#### 7. Transpose Data

**PLS-DA requires samples as rows**:
```r
msms_t <- t(msms_filtered)
```

#### 8. Perform PLS-DA

```r
ncomp_use <- 2  # Use 2 components for small datasets
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = ncomp_use)
```

**Component selection**:
- **2 groups**: Use 2-3 components
- **3+ groups**: Use up to min(groups, 5) components
- **Small n**: Never use more components than samples-per-group

#### 9. Extract Results

**Component scores** (for plotting):
```r
plsda_scores <- as.data.frame(plsda_result$variates$X)
colnames(plsda_scores) <- paste0("Comp", 1:ncol(plsda_scores))
plsda_scores$Condition <- annotation$Condition
```

**Variance explained**:
```r
plsda_var <- plsda_result$prop_expl_var$X * 100
```

**Loadings** (feature contributions):
```r
loadings <- as.data.frame(plsda_result$loadings$X)
```

#### 10. Cross-Validation

**Leave-One-Out Cross-Validation** (LOOCV):
```r
plsda_cv <- perf(plsda_result, 
                 validation = "loo",  # Critical for n=10!
                 progressBar = FALSE)
```

**Extract accuracy**:
```r
overall_error <- plsda_cv$error.rate$overall[, ncomp_use]
accuracy <- (1 - overall_error) * 100
```

**Why LOOCV?**
- Most reliable for small sample sizes
- Uses n-1 samples for training, 1 for testing
- Repeated n times
- Provides realistic performance estimate

#### 11. Visualization

**Main PLS-DA plot**:
```r
ggplot(plsda_scores, aes(x = Comp1, y = Comp2, color = Condition)) +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 1) +
  labs(
    title = "PLS-DA of MS/MS Data",
    subtitle = paste0("Accuracy: ", round(accuracy, 1), "%"),
    x = paste0("Component 1 (", round(plsda_var[1], 1), "%)"),
    y = paste0("Component 2 (", round(plsda_var[2], 1), "%)")
  ) +
  theme_bw()
```

**Additional visualizations**:
- PLS-DA with sample labels
- VIP score bar plot
- Loading plot
- 3D plot (if using 3 components)

#### 12. Feature Importance

**VIP (Variable Importance in Projection) scores**:
```r
vip_scores <- vip(plsda_result)
top_vip <- head(sort(vip_scores, decreasing = TRUE), 50)
```

**VIP interpretation**:
- **VIP > 1**: Important for discrimination
- **VIP > 1.5**: Very important
- **VIP < 0.8**: Less important

#### 13. Export Results

```r
# Save plots
ggsave("MSMS_PLSDA_plot.png", width = 10, height = 8, dpi = 300)

# Save scores
write.csv(plsda_scores, "MSMS_PLSDA_scores.csv", row.names = FALSE)

# Save top features
write.csv(top_vip, "MSMS_top_features.csv")
```

---

## Understanding the Output

### Generated Files

After running the analysis, you'll have:

#### 📊 Figures

1. **MSMS_PLSDA_plot.png/pdf**
   - Main PLS-DA visualization
   - Shows sample separation
   - Includes ellipses (95% confidence)
   - Displays accuracy

2. **MSMS_PLSDA_labeled.png**
   - Same as above with sample labels
   - Useful for identifying outliers

3. **MSMS_VIP_scores.png**
   - Bar plot of top 20 discriminant features
   - VIP > 1 threshold line

4. **MSMS_loading_plot.png**
   - Shows feature contributions
   - Highlights top discriminant features

5. **MSMS_PLSDA_3D.html** (if using 3 components)
   - Interactive 3D plot
   - Open in web browser

#### 📁 Data Files

1. **MSMS_PLSDA_scores.csv**
   - Component scores for each sample
   - Used for the plot

2. **MSMS_PLSDA_loadings.csv**
   - Feature loadings on each component
   - Shows how features contribute

3. **MSMS_PLSDA_top_features.csv**
   - Top 50 discriminant features
   - Ranked by VIP score

4. **MSMS_PLSDA_summary.csv**
   - Summary statistics
   - Variance explained
   - Classification accuracy

---

## Interpreting Results

### 1. PLS-DA Plot Interpretation

#### Good Separation Example
```
         Treatment
    ●●●●●●●●●●
   ●●●●●●●●●●●
  ●●●●●●●●●●●●
        |
        |_________ Control
              ●●●●●●●●
             ●●●●●●●●●
            ●●●●●●●●●
```
**Indicators of good separation**:
- ✅ Clear clusters for each condition
- ✅ Minimal overlap between ellipses
- ✅ Accuracy > 80%
- ✅ PC1 captures most variance

#### Poor Separation Example
```
     ●●●●●●●●●●●●●
    ●●●●●●●●●●●●●●●
   ●●●●●●●●●●●●●●●●
    ●●●●●●●●●●●●●●
     ●●●●●●●●●●●
    (Mixed clusters)
```
**Indicators of poor separation**:
- ❌ Overlapping clusters
- ❌ Accuracy < 60%
- ❌ No clear pattern
- ❌ Possible confounding factors

### 2. Classification Accuracy

**Interpretation guide**:

| Accuracy | Interpretation | Action |
|----------|---------------|--------|
| **>90%** | Excellent separation | Strong biomarkers likely |
| **80-90%** | Good separation | Promising results |
| **70-80%** | Moderate separation | Some discrimination |
| **60-70%** | Weak separation | Limited discrimination |
| **<60%** | Poor separation | Check data quality |

**Important notes**:
- For n=10, accuracy >80% is very good
- LOOCV accuracy is realistic (not inflated)
- Check if similar to permutation tests

### 3. Variance Explained

**Example**: Component 1 (45.2%), Component 2 (23.1%)

**What it means**:
- PC1 explains 45.2% of total variance
- PC2 explains additional 23.1%
- Together: 68.3% of variance

**Guidelines**:
- **>50%** (combined): Good representation
- **<30%** (combined): May need more components
- PC1 usually explains most variance

### 4. VIP Scores

**Variable Importance in Projection**:

```
Feature_123: VIP = 2.1  ← Very important
Feature_456: VIP = 1.3  ← Important
Feature_789: VIP = 0.7  ← Less important
```

**Interpretation**:
- **VIP > 1.5**: Strong discriminant features (high priority)
- **VIP 1.0-1.5**: Moderate importance
- **VIP < 1.0**: Low importance

**Use VIP scores to**:
1. Identify potential biomarkers
2. Select features for further validation
3. Guide targeted MS/MS experiments
4. Prioritize metabolite identification

### 5. Ellipses in PLS-DA Plot

**95% Confidence Ellipses**:
- Show expected sample distribution
- Based on multivariate normal assumption
- Larger ellipses = more within-group variability
- Smaller ellipses = more consistent group

**What if samples fall outside?**
- Possible outliers
- Check for technical issues
- Verify sample annotation
- May be biologically interesting

### 6. Sample Clustering Patterns

**Patterns to look for**:

1. **Perfect separation**: Groups completely distinct
   - Strong biological signal
   - High confidence in classification

2. **Partial overlap**: Some mixing
   - Biological heterogeneity
   - Shared features between groups

3. **Outlier samples**: Samples far from group
   - Technical errors
   - Mislabeling
   - Unique biological state

4. **Batch effects**: Clustering by batch not condition
   - Need better normalization
   - Include batch correction

---

## Important Considerations for Small Sample Sizes

### Special Guidelines for n=10

When working with small sample sizes (like 10 samples total), special care is needed:

#### 1. Cross-Validation is MANDATORY

```r
# ALWAYS use Leave-One-Out (LOO) for small n
plsda_cv <- perf(plsda_result, validation = "loo")
```

**Why?**
- Prevents overfitting
- Provides realistic accuracy
- Uses maximum training data
- Standard practice for n<20

**Never skip this step!**

#### 2. Limit Number of Components

```r
# For n=10 (5 per group):
# Use maximum 2-3 components
ncomp_use <- 2
```

**Rule of thumb**:
- n < 15: Use 2 components
- n = 15-30: Use 2-3 components
- n > 30: Use 3-5 components

**Why?**
- Prevents overfitting
- Each component uses degrees of freedom
- Need sufficient samples per component

#### 3. Feature Selection is Critical

```r
# Use only most variable features
n_features <- min(500, nrow(msms_data))
```

**Guidelines**:
- Features should be << samples
- Use 100-500 most variable features
- Remove low-variance features
- Focus on informative signals

#### 4. Report Uncertainty

Always report:
- ✅ Cross-validation accuracy
- ✅ Standard error (if available)
- ✅ Confidence intervals
- ✅ Sample size per group

#### 5. Validate Results

**Recommended validation approaches**:

1. **Permutation testing**:
```r
# Test if separation is better than random
n_perm <- 1000
perm_results <- vector("numeric", n_perm)

for(i in 1:n_perm) {
  shuffled_conditions <- sample(annotation$Condition)
  pls_perm <- plsda(msms_t, shuffled_conditions, ncomp = 2)
  perm_cv <- perf(pls_perm, validation = "loo")
  perm_results[i] <- 1 - perm_cv$error.rate$overall[, 2]
}

# Calculate p-value
real_accuracy <- (1 - overall_error)
p_value <- sum(perm_results >= real_accuracy) / n_perm
cat("Permutation p-value:", p_value, "\n")
```

2. **Independent validation set**:
- Collect new samples
- Apply trained model
- Check classification accuracy

3. **Targeted validation**:
- Validate top VIP features
- Use targeted MS/MS
- Confirm with standards

#### 6. Interpretation Caution

**Be conservative when**:
- Accuracy is 60-75% (modest)
- Large overlap in PLS-DA plot
- High within-group variability
- Limited biological replicates

**Be confident when**:
- Accuracy > 85% with LOOCV
- Clear visual separation
- VIP scores > 2 for top features
- Results consistent with biology

---

## Troubleshooting

### Common Issues and Solutions

#### Problem 1: Error - "Samples and conditions don't match"

**Error message**:
```
Error in plsda(msms_t, annotation$Condition) : 
  Length of Y must match number of rows in X
```

**Solution**:
```r
# Check dimensions
cat("Data samples:", nrow(msms_t), "\n")
cat("Annotation samples:", nrow(annotation), "\n")

# Verify names match
all(rownames(msms_t) == annotation$SampleID)

# Fix if needed
annotation <- annotation[match(rownames(msms_t), annotation$SampleID), ]
```

---

#### Problem 2: Poor Separation (Accuracy < 60%)

**Possible causes**:
1. Groups are actually similar
2. Too much noise in data
3. Batch effects
4. Need more features

**Solutions**:

**A. Check data quality**:
```r
# Plot total intensity per sample
total_intensity <- colSums(msms_data, na.rm = TRUE)
plot(total_intensity, main = "Sample Quality Check")
```

**B. Try different preprocessing**:
```r
# More stringent filtering
missing_threshold <- 0.3  # Stricter (30% instead of 50%)

# Different normalization
# Quantile normalization
library(preprocessCore)
msms_norm <- normalize.quantiles(as.matrix(msms_log))

# Pareto scaling
msms_pareto <- scale(msms_log, center = TRUE, scale = sqrt(apply(msms_log, 2, sd)))
```

**C. Remove outliers**:
```r
# Quick PCA to check
pca_check <- prcomp(t(msms_log), scale. = TRUE)
plot(pca_check$x[, 1:2])
# Identify and remove outliers, then re-run
```

**D. Check for batch effects**:
```r
# Color by batch in PCA
pca_df$Batch <- annotation$Batch
ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(Batch))) +
  geom_point(size = 3)
```

---

#### Problem 3: Too Many Missing Values

**Error**:
```
Warning: More than 70% missing values in feature matrix
```

**Solutions**:

**A. Less stringent filtering**:
```r
# Allow more missing values
missing_threshold <- 0.7  # Instead of 0.5
```

**B. Different imputation**:
```r
# K-nearest neighbors imputation
library(impute)
msms_imputed <- impute.knn(as.matrix(msms_data))$data

# Mean imputation
for(i in 1:nrow(msms_data)) {
  msms_data[i, is.na(msms_data[i,])] <- mean(msms_data[i,], na.rm = TRUE)
}
```

**C. Feature-wise imputation**:
```r
# Group-specific imputation
for(cond in unique(annotation$Condition)) {
  samples_in_group <- annotation$SampleID[annotation$Condition == cond]
  for(i in 1:nrow(msms_data)) {
    group_values <- msms_data[i, samples_in_group]
    if(any(is.na(group_values))) {
      msms_data[i, samples_in_group][is.na(group_values)] <- 
        min(group_values, na.rm = TRUE) / 2
    }
  }
}
```

---

#### Problem 4: Installation Issues

**mixOmics won't install**:

```r
# Try this sequence:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mixOmics", force = TRUE)

# If that fails, try specific version:
BiocManager::install("mixOmics", version = "3.18")

# Or from GitHub:
install.packages("devtools")
devtools::install_github("mixOmicsTeam/mixOmics")
```

**Package loading errors**:
```r
# Check package status
sessionInfo()

# Update all packages
update.packages(ask = FALSE)

# Reinstall specific package
remove.packages("mixOmics")
BiocManager::install("mixOmics")
```

---

#### Problem 5: Memory Issues

**Error**: "Cannot allocate vector of size..."

**Solutions**:

**A. Increase memory (if possible)**:
```r
# Check current memory
memory.limit()  # Windows only

# Increase memory
memory.limit(size = 16000)  # Windows: 16 GB
```

**B. Use fewer features**:
```r
# More aggressive feature selection
n_features <- 200  # Instead of 500

# Only use features with variance > threshold
var_threshold <- quantile(feature_variance, 0.75)
keep_features <- feature_variance > var_threshold
```

**C. Process in chunks**:
```r
# For very large datasets
chunk_size <- 1000
feature_chunks <- split(1:nrow(msms_data), 
                       ceiling(1:nrow(msms_data) / chunk_size))

for(chunk in feature_chunks) {
  # Process each chunk
  msms_chunk <- msms_data[chunk, ]
  # ... analysis ...
}
```

---

#### Problem 6: Ellipses Don't Appear

**Issue**: `stat_ellipse()` fails

**Solutions**:

**A. Check group sizes**:
```r
table(annotation$Condition)
# Need at least 3 samples per group for ellipses
```

**B. Alternative ellipse method**:
```r
# Use different type
stat_ellipse(level = 0.95, type = "t")  # Student's t
stat_ellipse(level = 0.95, type = "norm")  # Normal (default)
stat_ellipse(level = 0.95, type = "euclid")  # Euclidean
```

**C. Remove ellipses**:
```r
# Plot without ellipses
ggplot(plsda_scores, aes(x = Comp1, y = Comp2, color = Condition)) +
  geom_point(size = 5)
  # No stat_ellipse() call
```

---

#### Problem 7: VIP Scores Missing

**Error**: "Function 'vip' not found"

**Solution**:

**A. Load correct package**:
```r
# VIP is in mixOmics
library(mixOmics)
vip_scores <- vip(plsda_result)
```

**B. Manual VIP calculation**:
```r
# If vip() doesn't work, calculate manually
vip_manual <- function(plsda_obj) {
  W <- plsda_obj$loadings$X
  Q <- plsda_obj$loadings$Y
  p <- nrow(W)
  h <- ncol(W)
  
  vip_scores <- rep(0, p)
  for(i in 1:p) {
    w_sum <- sum(W[i, ]^2)
    vip_scores[i] <- sqrt(p * w_sum / h)
  }
  
  names(vip_scores) <- rownames(W)
  return(vip_scores)
}

vip_scores <- vip_manual(plsda_result)
```

---

#### Problem 8: Colors Don't Match Groups

**Issue**: Wrong colors assigned to groups

**Solution**:

```r
# Explicitly set color mapping
group_colors <- c("Control" = "#E63946",    # Red
                 "Treatment" = "#1D3557")   # Blue

# Check factor levels
levels(annotation$Condition)

# Ensure color names match exactly
names(group_colors) <- levels(annotation$Condition)

# Use in plot
scale_color_manual(values = group_colors)
```

---

## FAQ

### General Questions

**Q1: How many samples do I need for PLS-DA?**

**A:** Minimum recommendations:
- Absolute minimum: 6 samples (3 per group)
- Recommended minimum: 10 samples (5 per group)
- Good: 20-30 samples (10-15 per group)
- Ideal: >30 samples per group

For n=10, PLS-DA is appropriate but requires careful cross-validation.

---

**Q2: PLS-DA vs PCA - which should I use?**

**A:** Use both!

1. **Start with PCA**:
   - Check data quality
   - Identify outliers
   - Look for batch effects
   - See if natural clustering exists

2. **Then use PLS-DA**:
   - Maximize group separation
   - Identify discriminant features
   - Classify samples
   - Find biomarkers

---

**Q3: What's a good classification accuracy?**

**A:** Depends on context:

| Context | Good Accuracy |
|---------|--------------|
| Small n (<15) | >75% |
| Medium n (15-30) | >80% |
| Large n (>30) | >85% |
| Clinical biomarkers | >90% |

Always compare to:
- Random chance (50% for 2 groups)
- Permutation test results
- Independent validation

---

**Q4: How many features should I use?**

**A:** General rules:

- Features << samples (ideally features < samples/2)
- For n=10: Use 100-500 most variable features
- For n=30: Use 500-2000 features
- For n=100: Can use all features if needed

---

**Q5: Should I normalize my MS/MS data?**

**A:** Yes! Always normalize MS/MS data.

Recommended normalization methods:
1. **Median normalization** (recommended in script)
2. **Total Ion Current (TIC)** normalization
3. **Probabilistic Quotient Normalization (PQN)**
4. **Quantile normalization**

Also apply:
- Log2 transformation (reduces heteroscedasticity)
- Centering and/or scaling (for PLS-DA)

---

### Technical Questions

**Q6: What is cross-validation and why is it important?**

**A:** Cross-validation tests model performance on unseen data.

**Without cross-validation**:
- Model memorizes training data
- Accuracy is artificially high (overfitting)
- Results don't generalize

**With cross-validation** (LOOCV for n=10):
- Each sample tested independently
- Realistic accuracy estimate
- Detects overfitting

For n=10, LOOCV is essential!

---

**Q7: What are VIP scores?**

**A:** VIP (Variable Importance in Projection) scores measure feature importance.

**Interpretation**:
- VIP > 1: Important for separation
- VIP > 1.5: Very important
- VIP < 1: Less important

**Use VIP scores to**:
- Identify potential biomarkers
- Select features for validation
- Prioritize metabolite identification

---

**Q8: My accuracy is 100% - is this good?**

**A:** **Warning!** 100% accuracy with small n often indicates overfitting.

**Check**:
1. Did you use cross-validation?
2. Are you using too many components?
3. Are you using too many features?
4. Is there a data leak? (e.g., batch = condition)

**If cross-validated accuracy is 100%**:
- Groups might be truly distinct
- Run permutation test to verify
- Validate on independent data

---

**Q9: How do I handle batch effects?**

**A:** Several approaches:

1. **Include batch in model** (if balanced):
```r
# Design matrix with batch
design <- model.matrix(~ Condition + Batch, data = annotation)
```

2. **Batch correction before PLS-DA**:
```r
library(sva)
batch <- annotation$Batch
condition <- annotation$Condition
modcombat <- model.matrix(~condition)
msms_combat <- ComBat(dat = msms_log, batch = batch, mod = modcombat)
```

3. **Remove batch effect in design**:
- Balance conditions across batches
- Randomize sample order
- Include QC samples

---

**Q10: Can I use PLS-DA for >2 groups?**

**A:** Yes! PLS-DA works with multiple groups.

**Modifications needed**:
```r
# For 3 groups
annotation$Condition <- factor(annotation$Condition, 
                              levels = c("Group1", "Group2", "Group3"))

# Use more components
ncomp_use <- min(number_of_groups, 5)

# Pairwise comparisons for VIP scores
```

---

### Data-Specific Questions

**Q11: My data is from different MS platforms - can I combine them?**

**A:** Generally **not recommended** due to:
- Different mass accuracy
- Different resolution
- Platform-specific biases

**If you must combine**:
1. Normalize separately first
2. Use only common features
3. Apply batch correction
4. Validate thoroughly

Better: Analyze separately and compare results.

---

**Q12: Should I remove zeros or missing values?**

**A:** Handle carefully - zeros may be biological or technical.

**True zeros** (below detection limit):
- Impute with half-minimum
- Use in analysis

**Missing values** (not measured):
- Impute or remove feature
- Don't confuse with biological zeros

**Rule**:
- Remove features with >50% missing
- Impute remaining

---

**Q13: Can I use PLS-DA for quantitative predictions?**

**A:** PLS-DA is for classification, not quantification.

For quantitative outcomes, use:
- **PLS (not PLS-DA)**: Continuous response
- **Linear regression**: Simple relationships
- **Random Forest**: Complex relationships

---

**Q14: How do I report PLS-DA results in a paper?**

**A:** Include these elements:

1. **Methods**:
   - Data preprocessing steps
   - Number of features used
   - Number of components
   - Cross-validation method
   - Software versions

2. **Results**:
   - PLS-DA plot with ellipses
   - Variance explained per component
   - Cross-validated accuracy ± SE
   - Number of samples per group
   - Top discriminant features (VIP > 1)

3. **Statistics**:
   - Permutation test p-value
   - Confusion matrix
   - Feature validation (if available)

**Example sentence**:
"PLS-DA with leave-one-out cross-validation showed clear separation between Control and Treatment groups (accuracy: 85.0%, p<0.001 by permutation test), with the first two components explaining 68.3% of total variance."

---

**Q15: What's the difference between PLS-DA and OPLS-DA?**

**A:** 

| Feature | PLS-DA | OPLS-DA |
|---------|--------|---------|
| **Orthogonal filtering** | No | Yes |
| **Interpretation** | More complex | Easier |
| **Predictive components** | Multiple | One |
| **Orthogonal components** | Mixed | Separated |
| **Small n** | Good | Also good |

**OPLS-DA** separates:
- Predictive variation (related to Y)
- Orthogonal variation (unrelated to Y)

**Use OPLS-DA when**:
- You want clearer interpretation
- You have confounding variation
- Available in your software

Both are valid - PLS-DA is more common.

---

## References

### Key Papers

1. **PLS-DA Method**:
   - Barker, M., & Rayens, W. (2003). "Partial least squares for discrimination." *Journal of Chemometrics*, 17(3), 166-173.

2. **Cross-Validation**:
   - Westerhuis, J. A., et al. (2008). "Assessment of PLSDA cross validation." *Metabolomics*, 4(1), 81-89.

3. **VIP Scores**:
   - Chong, I. G., & Jun, C. H. (2005). "Performance of some variable selection methods when multicollinearity is present." *Chemometrics and Intelligent Laboratory Systems*, 78(1-2), 103-112.

4. **MS/MS Data Analysis**:
   - Gromski, P. S., et al. (2015). "A tutorial review: Metabolomics and partial least squares-discriminant analysis." *Analytica Chimica Acta*, 879, 10-23.

### Software Documentation

1. **mixOmics Package**:
   - Website: http://mixomics.org/
   - Vignettes: http://mixomics.org/methods/
   - GitHub: https://github.com/mixOmicsTeam/mixOmics

2. **ggplot2**:
   - Website: https://ggplot2.tidyverse.org/
   - Book: "ggplot2: Elegant Graphics for Data Analysis" by Hadley Wickham

3. **Bioconductor**:
   - Website: https://www.bioconductor.org/
   - Workflows: https://www.bioconductor.org/packages/release/BiocViews.html#___Workflow

### Online Resources

1. **MetaboAnalyst**: https://www.metaboanalyst.ca/
   - Web-based metabolomics analysis
   - Includes PLS-DA
   - Good for comparison

2. **Bioconductor Forum**: https://support.bioconductor.org/
   - Ask questions
   - Search solutions

3. **R for Data Science**: https://r4ds.had.co.nz/
   - Learn R basics
   - Data visualization

### Recommended Workflows

1. **Metabolomics Workbench**: https://www.metabolomicsworkbench.org/
   - Standard protocols
   - Example datasets

2. **XCMS Online**: https://xcmsonline.scripps.edu/
   - LC-MS/MS processing
   - Statistical analysis

---

## Citation

If you use these scripts in your research, please cite:

```
Rohart F, Gautier B, Singh A, Lê Cao KA (2017). mixOmics: An R package for 
'omics feature selection and multiple data integration. PLoS Computational 
Biology 13(11): e1005752.
```

For the complete MS/MS workflow:
```
Ayeh Bolouki (2025). PLS-DA Analysis Pipeline for MS/MS Data. 
GitHub: [https://github.com/AyehBlk/PLSDA-MSMS-Analysis/]
```

---

## Support and Contact

### Getting Help

1. **Check this README**: Most questions answered here
2. **Check troubleshooting section**: Common issues and solutions
3. **Review example code**: Complete working examples provided
4. **mixOmics documentation**: http://mixomics.org/
5. **Bioconductor support**: https://support.bioconductor.org/

### Reporting Issues

If you find bugs or have suggestions:

1. **Check if already reported**
2. **Provide reproducible example**
3. **Include**:
   - R version
   - Package versions (`sessionInfo()`)
   - Error message
   - Minimal code to reproduce

### Contributing

Contributions welcome! Please:
1. Fork repository
2. Create feature branch
3. Test thoroughly
4. Submit pull request

---

## License

This code is provided under the MIT License.

```
MIT License

Copyright (c) 2024

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## Acknowledgments

This pipeline was developed for MS/MS metabolomics and lipidomics analysis, adapted from gene expression analysis workflows. Special thanks to:

- The mixOmics team for excellent multivariate analysis tools
- The Bioconductor community for MS/MS data infrastructure
- Contributors to ggplot2 and tidyverse packages

---

## Version History

### Version 1.0.0 (2024)
- Initial release
- Basic PLS-DA workflow
- Optimized for small sample sizes (n=10)
- Comprehensive documentation
- Example data and scripts

### Planned Features
- Integration with xcms output
- Automated permutation testing
- Additional visualization options
- Batch effect correction module
- Multi-group comparison workflows

---

## Quick Reference Card

### Essential Commands

```r
# Load data
msms_data <- read.csv("feature_matrix.csv", row.names = 1)
annotation <- read.csv("sample_annotation.csv")

# Preprocess
msms_log <- log2(msms_data + 1)
msms_t <- t(msms_log)

# PLS-DA
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = 2)

# Cross-validate
plsda_cv <- perf(plsda_result, validation = "loo")

# Get accuracy
accuracy <- (1 - plsda_cv$error.rate$overall[, 2]) * 100

# Extract scores
plsda_scores <- as.data.frame(plsda_result$variates$X)

# Get VIP scores
vip_scores <- vip(plsda_result)

# Plot
ggplot(plsda_scores, aes(x = X1, y = X2, color = annotation$Condition)) +
  geom_point(size = 5) + stat_ellipse()
```

### File Structure
```
project/
├── feature_matrix.csv          # Your MS/MS data
├── sample_annotation.csv       # Sample information
├── MSMS_PLSDA_Analysis.R      # Main script
├── MSMS_PLSDA_QuickStart.R    # Quick version
└── results/                    # Output folder
    ├── MSMS_PLSDA_plot.png
    ├── MSMS_PLSDA_scores.csv
    └── MSMS_top_features.csv
```

---

**Last Updated**: October 2024  
**Author**: [Your Name]  
**Contact**: [Your Email]

---

*Happy analyzing! 🔬📊*
