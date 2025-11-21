# MS/MS Data Analysis with PLS-DA in R

![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-active-success)
![Bioconductor](https://img.shields.io/badge/Bioconductor-mixOmics-orange)

**Practical PLS-DA analysis pipeline for MS/MS metabolomics data, optimized for small sample sizes (n=10-30).**

![Application](https://img.shields.io/badge/Application-Metabolomics-purple) ![Method](https://img.shields.io/badge/Method-PLS--DA-orange) ![Platform](https://img.shields.io/badge/Platform-MS/MS-red)

---

##  Overview

This repository provides ready-to-use R scripts for performing **Partial Least Squares Discriminant Analysis (PLS-DA)** on tandem mass spectrometry (MS/MS) data. Specifically designed for **small sample sizes** commonly encountered in metabolomics studies.

###  Key Features

-  **Optimized for small n** - Works with as few as 10 samples (5 per group)
-  **Multiple entry points** - Comprehensive, quick-start, and data formatting scripts
-  **Publication-ready plots** - Beautiful ggplot2 visualizations
-  **Cross-validation** - Leave-one-out CV for small sample sizes
-  **Biomarker discovery** - VIP scores for feature ranking
-  **Platform-agnostic** - Works with xcms, MZmine, MS-DIAL, MetaboAnalyst, and more
-  **Comprehensive documentation** - Complete README with examples and troubleshooting

---

##  What is PLS-DA?

**PLS-DA (Partial Least Squares Discriminant Analysis)** is a supervised machine learning method that:

-  Maximizes separation between predefined groups (e.g., Control vs Treatment)
-  Reduces high-dimensional data to interpretable components
-  Identifies discriminant features (biomarkers)
-  Creates visualizations similar to PCA but with better classification

### PLS-DA vs PCA

| Feature | PCA | PLS-DA |
|---------|-----|--------|
| **Type** | Unsupervised | Supervised |
| **Goal** | Maximize variance | Maximize group separation |
| **Uses group info** |  No |  Yes |
| **Better for classification** |  No |  Yes |
| **Biomarker discovery** | Limited |  Excellent |

---

##  Quick Start

### Prerequisites

```r
# Install required packages
install.packages(c("ggplot2", "ggrepel", "caret"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")
```

### Minimal Example (30 seconds)

```r
library(ggplot2)
library(mixOmics)

# 1. Load data
msms_data <- read.csv("feature_matrix.csv", row.names = 1)
annotation <- read.csv("sample_annotation.csv")

# 2. Preprocess
msms_log <- log2(msms_data + 1)
msms_t <- t(msms_log)

# 3. PLS-DA
plsda_result <- plsda(msms_t, annotation$Condition, ncomp = 2)

# 4. Cross-validate
plsda_cv <- perf(plsda_result, validation = "loo")
accuracy <- (1 - plsda_cv$error.rate$overall[, 2]) * 100

# 5. Plot
plsda_scores <- as.data.frame(plsda_result$variates$X)
plsda_scores$Condition <- annotation$Condition

ggplot(plsda_scores, aes(x = X1, y = X2, color = Condition)) +
  geom_point(size = 5) +
  stat_ellipse(level = 0.95) +
  labs(title = paste0("PLS-DA (Accuracy: ", round(accuracy, 1), "%)")) +
  theme_bw()
```

**That's it!** You now have a PLS-DA plot with cross-validated accuracy.

---

##  Repository Structure

```
MSMS-PLSDA-Analysis/
├── MSMS_PLSDA_Analysis.R       # Comprehensive analysis pipeline
├── MSMS_PLSDA_QuickStart.R     # Minimal working example
├── Data_Format_Guide.R         # Data preparation from various platforms
├── README.md                   # Complete documentation
├── LICENSE                     # MIT License
└── examples/                   # Example data (optional)
    ├── feature_matrix.csv
    └── sample_annotation.csv
```

---

##  Documentation

### Three Scripts for Different Needs

1. **MSMS_PLSDA_Analysis.R** - Full pipeline including:
   - Data loading and quality control
   - Missing value imputation
   - Normalization (median, log transformation)
   - Feature selection (variance filtering)
   - PLS-DA with cross-validation
   - Multiple visualization types
   - VIP score calculation
   - Results export

2. **MSMS_PLSDA_QuickStart.R** - Minimal script:
   - Just the essentials
   - ~70 lines of code
   - Perfect for testing or simple analyses

3. **Data_Format_Guide.R** - Data preparation:
   - Examples from xcms, MZmine, MS-DIAL
   - Format conversion helpers
   - Annotation creation

### Complete README

See [full README](README.md) for:
-  Detailed methodology explanation
-  Step-by-step workflow
-  Output interpretation guide
-  Important considerations for small sample sizes
-  Troubleshooting common issues
-  FAQ (15+ questions answered)
-  References and citations

---

##  When to Use This

Perfect for:
- **Metabolomics studies** with limited samples
- **LC-MS/MS or GC-MS/MS** data analysis
- **Biomarker discovery** in biological samples
- **Quality control** of MS methods
- **Pilot studies** before larger experiments
- **Teaching/learning** multivariate analysis

Works with data from:
- xcms (R package)
- MZmine (Java software)
- MS-DIAL (Windows software)
- Compound Discoverer (Thermo)
- Progenesis QI (Waters)
- MetaboAnalyst exports
- Custom pipelines

---

##  Example Output

The analysis generates:
- **PLS-DA scores plot** with confidence ellipses
- **Cross-validation accuracy** metrics
- **VIP scores** for biomarker ranking
- **Loadings plots** showing feature contributions
- **Sample predictions** with probabilities
- **CSV exports** for further analysis

---

## 💡 Why This Implementation?

### Designed for Real-World MS/MS Data

 **Handles small sample sizes** - Uses leave-one-out cross-validation  
 **Missing value aware** - Proper imputation strategies  
 **Normalized for MS data** - Log transformation + median normalization  
 **Feature selection** - Variance filtering to reduce noise  
 **Publication-ready** - High-quality ggplot2 visualizations  
 **Well-documented** - Every step explained  

### Built on Proven Tools

- **mixOmics** - Leading multivariate analysis package
- **ggplot2** - Professional data visualization
- **caret** - Machine learning framework
- **Bioconductor** - Established bioinformatics infrastructure

---

##  Requirements

- **R** ≥ 4.0.0
- **Packages**: ggplot2, ggrepel, mixOmics, caret
- **RAM**: 4 GB minimum (8 GB recommended)
- **Input data**: Feature matrix + sample annotation (CSV format)

---

##  Learning Resources

### For Beginners
- Start with **MSMS_PLSDA_QuickStart.R**
- Read the [What is PLS-DA?](#what-is-pls-da) section
- Check [Quick Start Guide](#quick-start)

### For Advanced Users
- Use **MSMS_PLSDA_Analysis.R** for full control
- See [Detailed Workflow](README.md#detailed-workflow)
- Customize parameters in the script

### Understanding Results
- [Interpreting Results](README.md#interpreting-results) section
- [Understanding the Output](README.md#understanding-the-output) section
- [FAQ](README.md#faq) for common questions

---

##  Troubleshooting

**Common issues:**

| Problem | Solution |
|---------|----------|
| Error: "cannot find package" | Install required packages (see Prerequisites) |
| Low accuracy (<70%) | Check sample size, increase to n≥20 if possible |
| "Singular matrix" error | Reduce number of features (variance filtering) |
| All samples in one cluster | Groups may not be separable (biological issue) |
| Missing values | Check Data_Format_Guide.R for imputation |

See [full troubleshooting guide](README.md#troubleshooting) for more.

---

##  Citation

If you use these scripts in your research, please cite:

**mixOmics package:**
```
Rohart F, Gautier B, Singh A, Lê Cao KA (2017). mixOmics: An R package for 
'omics feature selection and multiple data integration. PLoS Computational 
Biology 13(11): e1005752.
```

**This repository:**
```
Ayeh Bolouki (2025). PLS-DA Analysis Pipeline for MS/MS Data. 
GitHub: https://github.com/[AyehBlk]/MSMS-PLSDA-Analysis
```

---

##  Contributing

Contributions welcome! Areas for improvement:
- Additional preprocessing options
- More visualization types
- Integration with other MS platforms
- Extended documentation
- Example datasets

Please open an issue or submit a pull request.

---

##  License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

**You are free to:**
-  Use for academic research
-  Use for commercial projects  
-  Modify and distribute
-  Include in your own work

---

##  Support

**Getting Help:**
1. Check the [FAQ section](README.md#faq)
2. Review [Troubleshooting](README.md#troubleshooting)
3. Read [mixOmics documentation](http://mixomics.org/)
4. Open an [issue](../../issues) on GitHub

---

## 👤 Author

**Ayeh Bolouki**
- GitHub: [@AyehBlk](https://github.com/AyehBlk)
- Role: Computational Biologist / Bioinformatician

---

##  Related Projects

- [PLS-DA from Scratch](https://github.com/your-username/PLSDA-Implementation) - Complete PLS-DA algorithm implementation


---

##  Star History

If you find this useful, please consider giving it a star! It helps others discover the project.

---

<p align="center">
  Made with ❤️ - Let's make free science for everybody around the world
</p>

<p align="center">
  <sub>Practical tools for metabolomics researchers everywhere</sub>
</p>

---

**Last Updated**: October 2024 
**Status**: Active Development
