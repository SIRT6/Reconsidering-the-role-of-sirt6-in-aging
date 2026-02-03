#!/usr/bin/env Rscript

########################################################
## Differential expression for human SIRT6 datasets ####
########################################################

library(arrow) # read/write parquet files
library(dplyr)
library(tibble)
library(DESeq2)

#########################
## Paths (arguments) ####
#########################

args <- commandArgs(trailingOnly = TRUE) # reads arguments passed to the script via a command line 

get_arg <- function(x) {
  i <- which(args == x)
  if (length(i) == 0) return(NULL)
  args[i + 1]
}

expr_path <- get_arg("--expr_path") # Count matrix
meta_path <- get_arg("--meta_path") # Meta data
out_dir <- get_arg("--out_dir") # Output dir

# Prevent from mistake
if (is.null(expr_path) || is.null(meta_path) || is.null(out_dir)) {
  stop("Usage: --expr_path --meta_path --out_dir")
}

organism <- "homo_sapiens"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), showWarnings = FALSE)

#########################
###### Load metadata ####
#########################

samples_meta <- read_parquet(
  file.path(meta_path, "samples.parquet")
) %>% column_to_rownames("sample_id")

sample2exp <- read_parquet(
  file.path(meta_path, "samples_to_experiment.parquet")
)

# Drop columns in meta consisting of NaNs only
samples_meta <- samples_meta[, colSums(is.na(samples_meta)) == 0]

#########################
###### Discover GSEs ####
#########################

# Paths to the experiment files
expr_files <- list.files(expr_path, pattern = "\\.parquet$", full.names = TRUE)
gse_ids <- sub("\\.parquet$", "", basename(expr_files))

#########################
# Filtering function ####
#########################

# Filter lowly expressed genes
# The filtering rule: gene is considered "expressed" in a group if it has at least 10 counts in enough samples of that group.
# Gene is kept if it passes these criteria in both of biological groups.
filter_counts <- function(counts, meta, min_count = 10) {
  
  keep <- rep(TRUE, nrow(counts))
  
  # Loop over biological groups 
  for (group in levels(meta$sirt6_status)) { 
    # Samples belonging to this group
    s <- rownames(meta) [meta$sirt6_status == group]
    
    # Minimum number of samples required to express the genes
    n_min <- max(2, floor(length(s) / 2)) # at least 2 or half of the samples in a group (whichever is larger)
    
    # Check expression criterion for this group
    keep_group <- rowSums(
      counts[, s, drop = FALSE] >= min_count
      ) >= n_min
    
    # Conservative rule: gene must pass in both of the groups
    keep <- keep & keep_group
  }
  
  # Return filtered count matrix
  counts[keep, , drop = FALSE]
}

#########################
########## Main loop ####
#########################

for (i in seq_along(expr_files)) { # One iteration = one GSE
  
  gse <- gse_ids[i]
  message("processing ", gse)
  
  tryCatch({
    
    # Load count matrix
    counts <- read_parquet(expr_files[i]) %>%
      column_to_rownames("gene_id")
    
    # Select only samples for this GSE
    samples <- sample2exp %>%
      filter(experiment_id == gse) %>%
      pull(sample_id)
    
    # Subsets metadata to match the count matrix
    meta <- samples_meta[samples, , drop = FALSE]
    
    ############################
    ## Derive SIRT6 status #####
    ############################
    
    meta$sirt6_status <- ifelse(meta$genotype == "WT", "control", "SIRT6") # WT -> control, anything else (KO, Het, OE) -> SIRT6
    meta$sirt6_status <- factor(meta$sirt6_status, levels = c("control", "SIRT6")) # control - reference level, SIRT6 - test group
    
    ##############################
    ## Define stratification #####
    ##############################
    
    # Split experiments by confounding factors (treatment/cell type) -> run separate DE per treatment/cell type
    
    strata <- list(all = rep(TRUE, nrow(meta))) # default: no splitting
    
    # Split the experiment by treatment 
    if ("treatment" %in% colnames(meta) && length(unique(meta$treatment)) > 1) {
      strata <- split(seq_len(nrow(meta)), meta$treatment)
    }
    
    # Split the experiment by cell type
    if ("cell_type" %in% colnames(meta) && length(unique(meta$cell_type)) > 1) {
      strata <- split(seq_len(nrow(meta)), meta$cell_type)
    }
    
    ############################
    ### Run DE per stratum #####
    ############################
    
    for (s in names(strata)) { # s is the name of stratum (DMSO, RAFi, etc.)
      
      idx <- strata[[s]] # indices of samples belonging to the stratum
      meta_s <- meta[idx, , drop = FALSE] # subset metadata to only these samples
      counts_s <- counts[, rownames(meta_s), drop = FALSE] # subset the count matrix to the same samples
      
      if (any(table(meta_s$sirt6_status) < 2)) next # if number of samples per group >=2 continue the analysis
      
      counts_f <- filter_counts(counts_s, meta_s) # filter lowly expressed genes
      if (nrow(counts_f) < 100) next # go next only if there are enough genes after filtering 
      
      # Build DESeq object
      dds <- DESeqDataSetFromMatrix(
        round(counts_f), 
        meta_s,
        design = ~ sirt6_status
      )
      
      # Run DESeq
      dds <- DESeq(dds, quiet = TRUE)
      
      # Extract results and perform DESeq2 independent filtering
      res <- results(
        dds,
        contrast = c("sirt6_status", 'SIRT6', "control") # compare SIRT6 vs control
      )
      
      # Convert DESeq2 results into a tidy data frame
      out <- as.data.frame(res) %>%
        rownames_to_column("gene_id")
      
      # Annotate results 
      out$experiment_id <- gse
      out$stratum <- s
      out$organism <- organism
      
      # Save files in a safer format (replace anything that is not letters, numbers, ., _, - with _)
      s_safe <- gsub("[^A-Za-z0-9._-]", "_", s)
      
      # Save DE results (in parquet format)
      write_parquet(
        out,
        file.path(
          out_dir,
          "results",
          paste0(gse, "_", s_safe, "_deseq2.parquet")
        )
      )
    }
    
  # Error handling (if any experiment failed)
  }, error = function(e) {
    cat(
      paste(Sys.time(), gse, e$message, "\n"),
      file = file.path(out_dir, "logs", "errors.log"),
      append = TRUE
    )
  })
}

message("Human DE analysis complete.")










