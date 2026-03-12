# Meta-analysis of SIRT6 knockout transcriptomic datasets

This directory contains a script for cross-species **meta-analysis of differential expression results** derived from SIRT6-related bulk RNA-seq datasets.

The analysis integrates differential expression results from multiple organisms to identify **genes consistently regulated by SIRT6 knockout (KO)**.

The script performs **gene-level meta-analysis across experiments and species**, estimates heterogeneity between studies, and generates several visualizations that highlight conserved transcriptional targets of SIRT6.

---

## Script overview

### `meta_analysis_script.R`

This script performs meta-analysis of differential expression results generated from the SIRT6 database.

Key features:

- Reads differential expression results aggregated across experiments
- Filters results to keep only **SIRT6 KO vs WT contrasts**
- Selects genes present in **multiple species and experiments**
- Performs **random-effects meta-analysis** for each gene using the `metafor` package
- Estimates **heterogeneity across studies** (τ² and I² statistics)
- Applies **Benjamini–Hochberg FDR correction** for multiple testing
- Identifies **significant SIRT6 transcriptional targets**
- Generates visualization plots to interpret the meta-analysis:
  - Volcano plot
  - Effect size vs heterogeneity plot
  - Forest plots for key genes
- Saves meta-analysis results and lists of significant genes

---

## Analysis logic

```text
For each gene, the script performs the following steps:

    1. Filter the dataset to retain only SIRT6 KO vs WT comparisons

    2. Identify genes present in at least:
        • 2 experiments
        • 2 species

    3. For each gene:
        • Collect log2 fold change (LFC) estimates and standard errors from all studies
        • Fit a random-effects meta-analysis model using metafor
        • Estimate the combined effect size (meta_LFC)
        • Calculate heterogeneity statistics (τ² and I²)

    4. Apply multiple testing correction (Benjamini–Hochberg FDR)

    5. Classify genes as:
        • Up-regulated
        • Down-regulated
        • Insignificant

    6. Generate visualization plots to interpret the results
```

## Requirements

### Software

- R (version ≥ 4.1 recommended)

### R packages

The following R packages are required:

- `dplyr`
- `data.table`
- `metafor`
- `ggplot2`
- `ggrepel`

Install missing packages using:

```r
install.packages(c("dplyr","data.table","metafor","ggplot2","ggrepel"))
```

---

## Input data

The script requires the following input files.

### Meta-analysis input table

`table_for_meta_analysis.csv`

This table contains differential expression results aggregated across all experiments.

Required columns:

| column | description |
|------|------|
human_gene_id | Ensembl human gene identifier |
organism | species name |
experiment_id | dataset identifier |
contrast | experimental comparison |
log2FoldChange | effect size |
lfcSE | standard error |
stratum | biological context |

---

### Ortholog mapping table

`ortholog_map_1to1.csv`

This file maps genes across species using **1:1 ortholog relationships**.

Required columns:

- `human_gene_id`
- `human_gene_symbol`

This mapping is used to annotate genes in plots.

---

### Experiment annotation table

`summary_DE_table.csv`

This table contains experiment metadata including:

- `experiment_id`
- `biological_system`

These annotations are used to describe biological context in forest plots.

---

## Usage

Run the script from the command line using **Rscript**.

### Example command

```bash
Rscript meta_analysis_script.R
```

The script reads input tables from predefined paths within the SIRT6 database.

---

## Output

The script produces several output files.

### Meta-analysis results

`meta_results_KO.csv`

This file contains gene-level meta-analysis results.

Columns include:

| column | description |
|------|------|
human_gene_id | gene identifier |
meta_LFC | combined effect size |
meta_SE | standard error |
CI_lower | lower confidence interval |
CI_upper | upper confidence interval |
pvalue | meta-analysis p-value |
FDR | adjusted p-value |
tau2 | between-study variance |
I2 | heterogeneity statistic |
k_experiments | number of experiments |
n_species | number of species |

---

### Significant genes

`meta_diffgenes_KO.csv`

Contains genes that satisfy:

- `FDR < 0.05`
- `|meta_LFC| > 0.58`

These genes represent **candidate conserved SIRT6 transcriptional targets**.

---

## Visualizations generated

The script generates three main plots.

---

### 1. Meta-analysis volcano plot

Displays gene-level meta-analysis results.

Axes:

- **x-axis:** meta log2 fold change  
- **y-axis:** −log10(FDR)

Interpretation:

- genes on the **right side** are up-regulated after SIRT6 knockout  
- genes on the **left side** are down-regulated after SIRT6 knockout  

Significant genes are highlighted and labeled.

Purpose:

Identify the strongest transcriptional responses to **SIRT6 loss**.

---

### 2. Effect size vs heterogeneity plot

Scatter plot showing:

- **x-axis:** meta_LFC  
- **y-axis:** heterogeneity (I²)

Interpretation:

| region | meaning |
|------|------|
low I² + strong effect | conserved SIRT6 targets |
high I² | species- or tissue-specific effects |

Genes with significant and conserved effects are labeled.

Purpose:

Identify **evolutionarily conserved transcriptional targets of SIRT6**.

---

### 3. Forest plots for key genes

Forest plots visualize study-level effects for selected genes:

- `H2AC7`
- `GAL`
- `ADGRE2`
- `PKP1`
- `HSPB1`
- `CCL18`

Each forest plot shows:

- effect size for each experiment
- confidence intervals
- species and biological system annotations
- combined meta-analysis estimate

Squares represent **individual study effects**.

The diamond represents the **combined meta-analysis estimate**.

Purpose:

Evaluate **consistency of gene regulation across species and experiments**.

---

## Biological interpretation

The meta-analysis identifies two main classes of SIRT6-regulated genes.

### Genes down-regulated after SIRT6 knockout

Examples:

- `H2AC7`
- `ADGRE2`
- `GAL`

These genes are associated with:

- chromatin organization
- immune signaling
- neuroendocrine regulation

---

### Genes up-regulated after SIRT6 knockout

Examples:

- `PKP1`
- `CCL18`
- `HSPB1`

These genes are involved in:

- stress response
- inflammation
- tissue remodeling

---

Overall, the results suggest that **SIRT6 promotes chromatin stability while suppressing stress and inflammatory pathways**, consistent with its known role in **genome maintenance and aging biology**.