# Quick Start Guide - MS/MS PLS-DA Analysis

##  5-Minute Start

### Step 1: Install Packages (2 minutes)

```r
# Install CRAN packages
install.packages(c("ggplot2", "ggrepel", "caret"))

# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")
```

### Step 2: Prepare Your Data (1 minute)

You need **two CSV files**:

**1. feature_matrix.csv** - Your MS/MS intensities
```
              ,Sample1,Sample2,Sample3,Sample4,Sample5
mz_100.05     ,1234.5 ,1456.2 ,1389.7 ,1298.4 ,1367.9
mz_200.10     ,567.8  ,601.3  ,589.2  ,578.9  ,595.1
mz_300.15     ,8901.2 ,8756.4 ,8823.9 ,8912.7 ,8834.5
```

**2. sample_annotation.csv** - Your sample groups
```
SampleID,Condition
Sample1,Control
Sample2,Control
Sample3,Treatment
Sample4,Treatment
Sample5,Treatment
```

### Step 3: Run Analysis (2 minutes)

**Option A: Quick Start Script**
```r
source("MSMS_PLSDA_QuickStart.R")
```

**Option B: Copy-Paste This:**
```r
library(ggplot2)
library(mixOmics)

# Load data
msms_data <- read.csv("feature_matrix.csv", row.names = 1)
annotation <- read.csv("sample_annotation.csv")

# Preprocess
msms_log <- log2(msms_data + 1)
msms_t <- t(msms_log)

# PLS-DA
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = 2)

# Cross-validate
plsda_cv <- perf(plsda_result, validation = "loo", progressBar = FALSE)
accuracy <- (1 - plsda_cv$error.rate$overall[, 2]) * 100

# Get scores
plsda_scores <- as.data.frame(plsda_result$variates$X)
plsda_scores$Condition <- annotation$Condition

# Plot
ggplot(plsda_scores, aes(x = X1, y = X2, color = Condition)) +
  geom_point(size = 5) +
  stat_ellipse(level = 0.95) +
  labs(title = paste0("PLS-DA (Accuracy: ", round(accuracy, 1), "%)")) +
  theme_bw()
```

**Done!** You should see a PLS-DA plot with your data.

---

## 📊 What You Get

After running the script:
-  PLS-DA scores plot (saved as PNG)
-  Cross-validated accuracy
-  Top 20 discriminant features (VIP scores)
-  Sample predictions

---

##  Troubleshooting

**Error: "cannot find function 'plsda'"**
```r
# Install mixOmics
BiocManager::install("mixOmics")
library(mixOmics)
```

**Error: "object not found"**
- Check your file names match exactly
- Make sure files are in working directory
- Use `getwd()` to check current directory

**Low accuracy (<70%)**
- Need more samples (try n≥20 total)
- Groups may be too similar biologically
- Try preprocessing data differently

**All points in one cluster**
- Groups might not be separable
- Check if you labeled samples correctly
- Try running PCA first to check for patterns

---

##  Next Steps

1. **For comprehensive analysis**: Use `MSMS_PLSDA_Analysis.R`
2. **Need to format data?**: See `Data_Format_Guide.R`
3. **Want to understand more?**: Read full [README.md](README.md)

---

##  Tips

- **Start with 2 components** (ncomp = 2) - easiest to visualize
- **Use leave-one-out CV** for small samples (n<30)
- **Check VIP scores > 1** for important features
- **Save your results** before closing R
- **Document your parameters** for reproducibility

---

##  Common Use Cases

**Comparing two conditions** (e.g., Disease vs Healthy)
```r
annotation <- data.frame(
  SampleID = colnames(msms_data),
  Condition = c(rep("Healthy", 5), rep("Disease", 5))
)
```

**Three or more groups**
```r
annotation <- data.frame(
  SampleID = colnames(msms_data),
  Condition = c(rep("Control", 4), rep("Treatment1", 3), rep("Treatment2", 3))
)
```

**Very small sample** (n=6, 3 per group)
```r
# Use leave-one-out CV (automatic for small n)
plsda_cv <- perf(plsda_result, validation = "loo")

# Consider using only 1 component
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = 1)
```

---

##  Need Help?

1. Check [Troubleshooting section](README.md#troubleshooting)
2. Read [FAQ](README.md#faq) 
3. Review [mixOmics documentation](http://mixomics.org/)
4. Open an issue on GitHub

---

**Good luck with your analysis!**
**lets make free science for everybody around the world.
