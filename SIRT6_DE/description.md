# Differential expression analysis of SIRT6 datasets

This directory contains script for differential expression (DE) analysis
of bulk RNA-seq datasets related to SIRT6.

The analysis is designed to work with the **SIRT6 database** stored on the lab server
and performs DESeq2-based differential expression analysis in a reproducible way.

---

## Script overview.

### `DE_script_human.R`

This script performs differential expression analysis for **human (Homo sapiens)**
SIRT6 datasets.

Key features:

- Reads count matrices stored in **Parquet** format
- Uses sample metadata from the SIRT6 database
- Collapses genotypes into two groups:
  - `control` (WT)
  - `SIRT6` (KO / Het / OE / mutant)
- Stratifies experiments by **treatment** or **cell type** when applicable
- Applies conservative low-expression filtering
- Uses **DESeq2** for differential expression analysis
- Writes results in **Parquet** format for downstream meta-analysis

---

## Requirements

### Software

- R (version ≥ 4.1 recommended)

### R packages

The following R packages are required:

- `arrow`
- `dplyr`
- `tibble`
- `DESeq2`

---

## Input data

The script expects the SIRT6 database to have the following structure:

```text
SIRT6_db/
├── expression/
│   └── homo_sapiens/
│       └── <GSE_ID>.parquet
└── metadata/
    └── homo_sapiens/
        ├── samples.parquet
        ├── experiments.parquet
        └── samples_to_experiment.parquet
```

---

## Usage 

Run the script from the command line using Rscript.

### Example command

```bash
Rscript DE_script_human.R \
  --expr_path /path/to/directory/with/count/matrices/for/human \
  --meta_path /path/to/directory/with/human/metadata \
  --out_dir /path/to/output/homo_sapiens
```

---

## Output

The script creates the following directory structure in the specified output path:
```text
out_dir/
├── results/
│   └── <GSE_ID>_<stratum>_deseq2.parquet
└── logs/
    └── errors.log
```
### Output description

- results/

Contains differential expression results for each experiment and stratum.
Each file is a Parquet table with DESeq2 statistics (log2FC, p-value, adjusted p-value).
- logs/errors.log

Contains a record of experiments or strata that failed during processing,
along with error messages.



