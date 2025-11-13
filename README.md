---
title: "README"
author: "Siyuan.Guo"
date: "2025-11-11"
output: html_document
editor_options: 
  markdown: 
    wrap: 72
---

## 1. Tool Overview & Purpose

`m6APrediction` is a user-friendly, efficient R package designed
specifically for predicting m6A methylation sites in RNA sequences. It
integrates multi-dimensional biological features and automatic DNA 5-mer
sequence encoding, supporting two core prediction scenarios:
**single-sample rapid analysis** and **batch-sample high-throughput
processing**.

Leveraging the `randomForest` algorithm, the package automatically
handles input feature preprocessing (e.g., factor encoding for
categorical variables) and outputs both predicted m6A methylation
probabilities (ranging from 0 to 1) and binary classification statuses
("Positive" or "Negative"). It is tailored for researchers in
epigenetics, RNA biology, and molecular genetics to streamline m6A site
screening workflows without requiring extensive computational expertise.

### Core Functions

-   `prediction_single()`: Predicts m6A status for a single sample using
    individual feature inputs, returning a named vector with prediction
    results.
-   `prediction_multiple()`: Batch-predicts m6A status for multiple
    samples using a feature data frame, returning an extended data frame
    that includes original features and prediction outputs.
-   Internal utility: `dna_encoding()` (non-exported): Automatically
    encodes DNA 5-mer sequences into factor-type features, eliminating
    the need for manual preprocessing.

## 2. Installation Guide

### Dependencies

The package requires two core dependencies, which are automatically
installed during setup (no manual intervention needed): -
`randomForest`: Core algorithm for model training and prediction. -
`stats`: Provides the `predict()` function (built-in R package; no
additional installation required).

### Install from GitHub (Recommended)

If you have uploaded the package to GitHub, users can install it
directly using `devtools` or `remotes`:

```{r}
# Install devtools if not already installed
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("[ttt66235]/m6APrediction")

# Alternative: Install via remotes
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("[ttt66235]/m6APrediction")
```

## 3. Quick Start Examples

Before running the examples, ensure the package's extdata folder
contains the following files: rf_fit.rds: Pre-trained randomForest
model. m6A_input_example.csv: Example multi-sample feature data (must
include 7 required columns).

### 3.1 Load the Package

```{r}
library(m6APrediction)
```

### 3.2 Single Sample Prediction (prediction_single())

Input individual feature values for a single sample to obtain prediction
results:

```{r}
# Load the pre-trained model (included in the package)
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

# Run single-sample prediction
single_result <- prediction_single(
  ml_fit = ml_fit,
  gc_content = 0.6,                  # Numeric: GC content (0-1 range)
  RNA_type = "mRNA",                 # Character: RNA type (options: mRNA/lincRNA/lncRNA/pseudogene)
  RNA_region = "CDS",                # Character: RNA region (options: CDS/intron/3'UTR/5'UTR)
  exon_length = 12,                  # Numeric: Exon length (positive integer)
  distance_to_junction = 50,         # Numeric: Distance to junction (non-negative value)
  evolutionary_conservation = 0.8,   # Numeric: Evolutionary conservation (0-1 range)
  DNA_5mer = "ATGAT",                # Character: 5-mer DNA sequence (A/T/C/G only; 5 characters)
  positive_threshold = 0.5           # Numeric: Classification threshold (0-1; default = 0.5)
)

# View results (predicted probability + classification status)
print(single_result)
```

``` text
predicted_m6A_prob predicted_m6A_status 
            0.723                  Positive 
```

### 3.3 Multiple Samples Prediction (prediction_multiple())

Batch-predict m6A status for multiple samples using a feature data
frame:

```{r}
# Load model and example multi-sample data
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))

# Inspect input data structure (verify 7 required columns exist)
head(feature_df[, c("gc_content", "RNA_type", "RNA_region", "DNA_5mer")])

# Run batch prediction
multiple_result <- prediction_multiple(
  ml_fit = ml_fit,
  feature_df = feature_df,
  positive_threshold = 0.5
)

# View key results (original features + prediction columns)
head(multiple_result[, c("RNA_type", "RNA_region", "predicted_m6A_prob", "predicted_m6A_status")])
```

#### Example Output：

| RNA_type | RNA_region | predicted_m6A_prob | predicted_m6A_status |
|----------|------------|--------------------|----------------------|
| mRNA     | CDS        | 0.68               | Positive             |
| lncRNA   | 3'UTR      | 0.41               | Negative             |
| lincRNA  | intron     | 0.75               | Positive             |
| mRNA     | 5'UTR      | 0.39               | Negative             |

## 4. Model Performance

The randomForest model integrated in m6APrediction was validated using
an independent test set, demonstrating strong predictive power with the
following key metrics: ROC Curve: Area Under Curve (AUC) = 0.89 PRC
Curve: Average Precision (AP) = 0.70

### Performance Visualization

Place your ROC and PRC curve images in a figures folder with the
following filenames: roc_curve.png prc_curve.png

The images will be displayed below:

![Figure 1: ROC curve of the m6A prediction model (AUC = 0.8854). A
higher AUC indicates superior ability to distinguish between positive
and negative samples.](figures/ROC_curve.png)

![Figure 2: Precision-Recall Curve (PRC) of the m6A prediction model (AP
= 0.7009). A higher AP reflects better precision for positive sample
prediction—critical for imbalanced m6A datasets.](figures/PRC_curve.png)

## 5. Important Notes

### Required Input Columns:

For prediction_multiple(), the input feature_df must contain 7 mandatory
columns: gc_content, RNA_type, RNA_region, exon_length,
distance_to_junction, evolutionary_conservation, and DNA_5mer. Missing
columns will trigger an error.

### DNA 5-mer Sequence Constraints:

Must be exactly 5 characters long, containing only A/T/C/G (no spaces,
lowercase letters, or special characters).

### Categorical Feature Levels:

RNA_type supports 4 levels ("mRNA", "lincRNA", "lncRNA", "pseudogene")
and RNA_region supports 4 levels ("CDS", "intron", "3'UTR", "5'UTR").
Other values are treated as invalid.

### Threshold Adjustment:

Modify positive_threshold based on research needs (e.g., set to 0.3 for
higher sensitivity, or 0.7 for higher specificity).

### Troubleshooting:

If installation fails, manually install dependencies first with:
install.packages(c("randomForest", "stats"))

## 6. Contact

For questions, issues, or feedback, please contact:
[Siyuan.Guo23\@student.xjtlu.edu.cn](mailto:Siyuan.Guo23@student.xjtlu.edu.cn)
