graphics.off()
rm(list = ls())

library(Seurat)
library(Matrix)
library(dplyr)
library(purrr)
library(R.utils)
library(SingleCellExperiment)
library(schard)
library(condiments)
library(SingleCellExperiment)
library(slingshot)
library(tibble)
library(dplyr)
library(scales)
library(irlba)
library(scuttle)
library(sceasy)
library(monocle3)
library(SeuratWrappers)
library(stringr)
library(batchelor)
library(scran)
library(scater)
library(tradeSeq)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(bluster)
library(biomaRt)
library(SingleR)
library(muscat)
library(tidyverse)
library(cowplot)
library(smplot2)
library(limma)
library("Libra")
library('presto')



## GSE157827 (Lau_et_al_2022) ----

sample_ids <- c("/tank/projects/public_data/Lau_et_al_2022/GSM4775561_AD1",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775562_AD2",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775563_AD4",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775564_AD5",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775565_AD6",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775566_AD8",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775567_AD9",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775568_AD10",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775569_AD13",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775570_AD19",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775571_AD20",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775572_AD21",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775573_NC3",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775574_NC7",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775575_NC11",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775576_NC12",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775577_NC14",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775578_NC15",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775579_NC16",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775580_NC17",
                "/tank/projects/public_data/Lau_et_al_2022/GSM4775581_NC18")


# Function to handle duplicate gene names
read_10x_gzipped_fixed <- function(sample_id) {
  cat("Processing:", sample_id, "\n")
  
  # Read barcodes
  barcodes <- readLines(gzfile(paste0(sample_id, "_barcodes.tsv.gz")))
  
  # Read features and handle duplicates
  features <- read.table(gzfile(paste0(sample_id, "_features.tsv.gz")), 
                         header = FALSE, stringsAsFactors = FALSE)
  
  # Check for duplicate gene names
  gene_names <- features[, 2]
  if (any(duplicated(gene_names))) {
    cat("Found", sum(duplicated(gene_names)), "duplicate gene names. Making unique...\n")
    gene_names <- make.unique(gene_names, sep = "_")
  }
  
  # Read matrix
  matrix_data <- readMM(gzfile(paste0(sample_id, "_matrix.mtx.gz")))
  
  # Set dimensions with unique gene names
  rownames(matrix_data) <- gene_names
  colnames(matrix_data) <- paste(sample_id, barcodes, sep = "_")
  
  return(matrix_data)
}

# Process each sample
seurat_objects <- list()

for (sample_id in sample_ids) {
  tryCatch({
    counts_matrix <- read_10x_gzipped_fixed(sample_id)
    
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      project = sample_id,
      min.cells = 3,
      min.features = 200
    )
    
    # Extract condition (AD or NC)
    if (grepl("_AD", sample_id)) {
      seurat_obj$condition <- "AD"
      seurat_obj$sample_number <- gsub(".*_AD(\\d+)$", "\\1", sample_id)
    } else if (grepl("_NC", sample_id)) {
      seurat_obj$condition <- "NC"
      seurat_obj$sample_number <- gsub(".*_NC(\\d+)$", "\\1", sample_id)
    }
    
    seurat_obj$sample <- basename(sample_id)  # optional: just the GSM ID
    
    seurat_objects[[sample_id]] <- seurat_obj
    cat("✓ Successfully processed", sample_id, "\n")
    
  }, error = function(e) {
    cat("✗ Error processing", sample_id, ":", e$message, "\n")
  })
}


# Merge all successful samples
if (length(seurat_objects) > 0) {
  combined_seurat <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[-1],
    project = "GSE157827_combined"
  )
  
  # Save the combined object
  saveRDS(combined_seurat, "/tank/projects/public_data/Lau_et_al_2022/GSE157827_combined_seurat.rds")
  
  cat("\n=== SUMMARY ===\n")
  cat("Total cells:", ncol(combined_seurat), "\n")
  cat("Total genes:", nrow(combined_seurat), "\n")
  cat("Samples processed:", length(seurat_objects), "/", length(sample_ids), "\n")
  cat("Output file: GSE157827_combined_seurat.rds\n")
  
} else {
  stop("No samples were processed successfully")
}


### Pseudobulk ----

sce <- readRDS('/tank/projects/public_data/Lau_et_al_2022/GSE157827_sce.rds')

assays(sce)

sce <- logNormCounts(sce)
sce <- runUMAP(sce)

unique(sce$id)

plotUMAP(sce, 
         colour_by = "celltype") +
  theme_minimal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

unique(sce$sample_number)

(sce <- prepSCE(sce, 
                kid = "celltype", # subpopulation assignments
                gid = "condition",  # group IDs (ctrl/stim)
                sid = "sample",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

levels(sce$group_id)

sce$group_id <- relevel(sce$group_id, ref = "NC")

## Pseudobulk aggregate

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))


library(limma)
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("AD-NC", levels = mm)

# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast)


# run DS analysis
#res <- pbDS(pb, method='DESeq2')

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[6])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size=5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(-2, 2) +
  coord_flip() +
  labs(title = 'OPC',
       x = "",
       y = "Log2FoldChange with 95% CI") +
  theme_minimal() +
  scale_color_manual(values = c("SIRT1" = "#E41A1C", 
                                "SIRT2" = "#377EB8",
                                "SIRT3" = "#4DAF4A", 
                                "SIRT4" = "#984EA3",
                                "SIRT5" = "#FF7F00", 
                                "SIRT6" = "coral2",
                                "SIRT7" = "#A65628")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(
      size = 15, 
      colour = ifelse(sirt$padj < 0.05, "red", "black"),
      face = ifelse(sirt$padj < 0.05, "bold", "plain")
    ),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    aspect.ratio = 1/3,
    legend.position = "none",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )




## GSE188545 (Zhang_et_al_2023) ----

sample_ids <- c("GSM5685287_AD02MTG", "GSM5685288_AD12MTG", "GSM5685289_AD30MTG",
                "GSM5685290_AD04MTG", "GSM5685291_AD16MTG", "GSM5685292_AD17MTG",
                "GSM5685293_HC14MTG", "GSM5685294_HC19MTG", "GSM5685295_HC35MTG",
                "GSM5685296_HC03MTG", "GSM5685297_HC07MTG", "GSM5685298_HC37MTG")

# Set the base directory path
base_dir <- "/tank/projects/public_data/Zhang_et_al_2023/"

# Function to handle duplicate gene names
read_10x_gzipped_fixed <- function(sample_id) {
  cat("Processing:", sample_id, "\n")
  
  # Construct full file paths
  barcodes_file <- paste0(base_dir, sample_id, "_barcodes.tsv.gz")
  features_file <- paste0(base_dir, sample_id, "_genes.tsv.gz")
  matrix_file <- paste0(base_dir, sample_id, "_matrix.mtx.gz")
  
  # Read barcodes
  barcodes <- readLines(gzfile(barcodes_file))
  
  # Read features and handle duplicates
  features <- read.table(gzfile(features_file), 
                         header = FALSE, stringsAsFactors = FALSE)
  
  # Check for duplicate gene names
  gene_names <- features[, 2]
  if (any(duplicated(gene_names))) {
    cat("Found", sum(duplicated(gene_names)), "duplicate gene names. Making unique...\n")
    gene_names <- make.unique(gene_names, sep = "_")
  }
  
  # Read matrix
  matrix_data <- readMM(gzfile(matrix_file))
  
  # Set dimensions with unique gene names
  rownames(matrix_data) <- gene_names
  colnames(matrix_data) <- paste(sample_id, barcodes, sep = "_")
  
  return(matrix_data)
}

# Process each sample
seurat_objects <- list()

for (sample_id in sample_ids) {
  tryCatch({
    counts_matrix <- read_10x_gzipped_fixed(sample_id)
    
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      project = sample_id,
      min.cells = 3,
      min.features = 200
    )
    
    # Extract condition (AD or HC) and sample details
    if (grepl("^GSM568528[7-9]|^GSM568529[0-2]", sample_id)) {  # AD samples: GSM5685287-5292
      seurat_obj$condition <- "AD"
      seurat_obj$sample_id <- gsub(".*_(AD\\d+MTG)$", "\\1", sample_id)
    } else if (grepl("^GSM568529[3-8]", sample_id)) {  # HC samples: GSM5685293-5298
      seurat_obj$condition <- "HC"
      seurat_obj$sample_id <- gsub(".*_(HC\\d+MTG)$", "\\1", sample_id)
    }
    
    seurat_obj$sample <- sample_id  # full GSM ID
    
    # Extract the numeric part for sample number
    seurat_obj$sample_number <- gsub(".*(AD|HC)(\\d+)MTG", "\\2", sample_id)
    
    seurat_objects[[sample_id]] <- seurat_obj
    cat("✓ Successfully processed", sample_id, "\n")
    
  }, error = function(e) {
    cat("✗ Error processing", sample_id, ":", e$message, "\n")
  })
}

# Merge all successful samples
if (length(seurat_objects) > 0) {
  combined_seurat <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[-1],
    project = "GSE188545_combined"
  )
  
  # Save the combined object
  saveRDS(combined_seurat, "/tank/projects/public_data/Zhang_et_al_2023/GSE188545_combined_seurat.rds")
  
  cat("\n=== SUMMARY ===\n")
  cat("Total cells:", ncol(combined_seurat), "\n")
  cat("Total genes:", nrow(combined_seurat), "\n")
  cat("Samples processed:", length(seurat_objects), "/", length(sample_ids), "\n")
  cat("AD samples:", sum(combined_seurat$condition == "AD"), "\n")
  cat("HC samples:", sum(combined_seurat$condition == "HC"), "\n")
  cat("Output file: GSE188545_combined_seurat.rds\n")
  
} else {
  stop("No samples were processed successfully")
}


sc <- readRDS("/tank/projects/public_data/Zhang_et_al_2023/GSE188545_combined_seurat.rds")

sc[["RNA"]] <- JoinLayers(sc[["RNA"]])

sceasy::convertFormat(sc, from="seurat", to="sce",
                      outFile='/tank/projects/public_data/Zhang_et_al_2023/GSE188545_sce.rds')

### Pseudobulk ----

sce <- readRDS('/tank/projects/public_data/Zhang_et_al_2023/GSE188545_sce_ann.rds')

sce@assays@data$logcounts

sce <- logNormCounts(sce)
sce <- runUMAP(sce)

unique(sce$condition)

plotUMAP(sce, 
         colour_by = "celltype") +
  theme_minimal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

unique(sce$condition)

(sce <- prepSCE(sce, 
                kid = "celltype", # subpopulation assignments
                gid = "condition",  # group IDs (ctrl/stim)
                sid = "sample_id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

levels(sce$group_id)

sce$group_id <- relevel(sce$group_id, ref = "HC")

## Pseudobulk aggregate

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))


library(limma)
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("AD-HC", levels = mm)

# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast)


# run DS analysis
#res <- pbDS(pb, method='DESeq2')

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[9])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size=5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(-2, 2) +
  coord_flip() +
  labs(title = 'OPC',
       x = "",
       y = "Log2FoldChange with 95% CI") +
  theme_minimal() +
  scale_color_manual(values = c("SIRT1" = "#E41A1C", 
                                "SIRT2" = "#377EB8",
                                "SIRT3" = "#4DAF4A", 
                                "SIRT4" = "#984EA3",
                                "SIRT5" = "#FF7F00", 
                                "SIRT6" = "coral2",
                                "SIRT7" = "#A65628")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    aspect.ratio = 1/3,
    legend.position = "none",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

## Cell type annotation ----

sce <- readRDS("/tank/projects/public_data/Zhang_et_al_2023/GSE188545_sce_ann.rds")
sce <- logNormCounts(sce) 

sce2 <- readRDS("/tank/projects/public_data/Siletti_et_al_2023/cx_mtg_sce.rds")
sce2 <- sce2[,colSums(counts(sce2)) > 0]
sce2 <- logNormCounts(sce2) 

rownames(sce)

# Connect to ENSEMBL database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get gene symbols
gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(sce2),
  mart = ensembl
)

# Create a mapping vector
symbol_mapping <- setNames(gene_symbols$hgnc_symbol, gene_symbols$ensembl_gene_id)

# Convert to uppercase and replace row names
rownames(sce2) <- toupper(symbol_mapping[rownames(sce2)])

# Handle any missing mappings (NA values)
rownames(sce2)[is.na(rownames(sce2))] <- "UNKNOWN"

pred <- SingleR(test=sce, ref=sce2, labels=sce2$cell_type)
colData(sce)$celltype <- pred$labels
table(pred$labels)

plotScoreHeatmap(pred)



sc <- readRDS("/tank/projects/public_data/Zhang_et_al_2023/GSE188545_combined_seurat.rds")

sc[["RNA"]] <- JoinLayers(sc[["RNA"]])


cat("Available layers after joining:\n")
print(SeuratObject::Layers(sc[["RNA"]]))

if (all(dim(SeuratObject::LayerData(sc[["RNA"]], layer = "data")) == 0)) {
  cat("Data layer is empty - this is normal if normalization wasn't performed\n")
}

if (all(dim(SeuratObject::LayerData(sc[["RNA"]], layer = "scale.data")) == 0)) {
  cat("Scale.data layer is empty - this is normal if scaling wasn't performed\n")
}


sce <- Seurat::as.SingleCellExperiment(sc)
sce <- logNormCounts(sce)

cat("Available metadata columns:\n")
print(colnames(colData(sce)))

cat("\nFirst few rows of metadata:\n")
print(head(colData(sce)))

cat("\nUnique sample identifiers:\n")
print(unique(sce$sample))  # or whatever your sample column is called

rownames(sce)

library(scRNAseq)

searchDatasets(
  defineTextQuery("GRCh38", field="genome") &
    (defineTextQuery("neuro%", partial=TRUE))
)[,c("name", "title", "version", "path", "genome")]

#sce2 <- fetchDataset("zhong-prefrontal-2018", "2023-12-22")

sce2 <- fetchDataset("nowakowski-cortex-2017", "2023-12-22")
unique(sce2$Laminae)

sce2 <- sce2[, colSums(counts(sce2)) > 0]

# Recompute size factors (now all will be positive)
sizeFactors(sce2) <- librarySizeFactors(sce2)

# Verify all size factors are positive
summary(sizeFactors(sce2))  # Min. should now be > 0

# Now run logNormCounts successfully
sce2 <- logNormCounts(sce2)

library(SingleR)

pred <- SingleR(test=sce, ref=sce2, labels=sce2$cell_types)
table(pred$labels)

plotScoreHeatmap(pred)



colData(sce)$celltype <- pred$labels

# UMAP with cell type colors
sce <- runUMAP(sce)

rare_labels <- c("Neurons", "Stem cells")
sce <- sce[, !pred$labels %in% rare_labels]
plotUMAP(sce, colour_by = "celltype_pred")


saveRDS(sce, '/tank/projects/public_data/Lau_et_al_2022/GSE157827_sce.rds')



sce <- readRDS('/tank/projects/public_data/Lau_et_al_2022/GSE157827_sce_ann.rds')

table(sce$celltype)

plotUMAP(sce, 
         colour_by = "celltype_pred") +
  theme_minimal() +
  labs(title = "Cell Type UMAP",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme(legend.position = "right")


# Otero-Garcia et al. 2020 ----

sce <- readRDS('/tank/projects/public_data/Otero-Garcia_et_al_2020/GSE129308_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Create a frequency table
freq_table <- as.data.frame(table(
  CellType = sce$celltype,
  Condition = sce$condition
))

# Find cell types with only 1 cell in any condition
problem_clusters <- freq_table %>%
  filter(Freq == 1) %>%
  pull(CellType) %>%
  unique()

print("Clusters with only 1 cell in any condition:")
print(problem_clusters)

if (length(problem_clusters) > 0) {
  sce_clean <- sce[, !sce$celltype %in% problem_clusters]
  
  sce_clean$celltype <- droplevels(factor(sce_clean$celltype))
  
  cat(sprintf("Removed %d problematic clusters.\n", length(problem_clusters)))
  cat("Remaining clusters:", paste(unique(sce_clean$celltype), collapse = ", "), "\n")
}

unique(sce$sample_id)
DE = run_de(sce_clean, cell_type_col = "celltype", label_col = "condition", replicate_col = "sample_id",
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)

DV = calculate_delta_variance(sce_clean, cell_type_col = "celltype", label_col = "condition", replicate_col = "sample_id")
DV.df <- as.data.frame(DV$oligodendrocyte)

de_sig <- subset(DE, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/OteroGarcia_wilcox_test_DE.csv')

DE <- read.csv('/tank/projects/ekashuk/AD/scRNAseq/OteroGarcia_wilcox_test_DE.csv', row.names = 1)
de_sig <- subset(DE, p_val_adj < 0.05)

sirt <- subset(DE, grepl('SIRT', gene))
sirt

library(tidyverse)
library(ggplot2)

# prepare
df <- sirt %>%
  mutate(
    p_adj = p_val_adj,
    negLog10Padj = -log10(p_adj + 1e-300),                   # avoid Inf
    sig = p_adj < 0.05,
    pct = pmax(AD.pct, Control.pct)                         # size = max prevalence
  ) %>%
  arrange(cell_type)

# bubble plot
ggplot(df, aes(x = cell_type, y = gene)) +
  geom_point(
    aes(
      size = pct,
      fill = avg_logFC,
      color = sig
    ),
    shape = 21,
    stroke = 0.8
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "Percent expressed") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3.5, 3.5),
    name = "avg_logFC"
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    name = "Significant"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 18, colour = "black"),
    axis.text.x  = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    #aspect.ratio   = 3/1
  ) +
  labs(
    x = "",
    y = ""
  )






sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)


p1 <- ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size = 7) +   # outlined points
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(title = "DLPFC (BA9)",
       x = NULL,
       y = "Log2 Fold Change (95% CI)") +
  scale_color_manual(values = c(
    "SIRT1" = "#FF6F91",  
    "SIRT2" = "#6FAFFF",  
    "SIRT3" = "#6FFF8F",  
    "SIRT4" = "#C66FFF",  
    "SIRT5" = "#FFB36F",  
    "SIRT6" = "#FF8F6F",  
    "SIRT7" = "#6FFFEF"   
  )) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.y  = element_text(size = 20, colour = "black"),
    axis.text.x  = element_blank(),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    #legend.position = "none",
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio   = 4/1
  )




table(sce$celltype)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$condition    
)


plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "AD")))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

p <- plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression)) +
  geom_density(aes(y = after_stat(scaled), color = CellType), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ Condition, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0,6)) +
  labs(
       color = "Cell type") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p)

## Distributions

sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$condition    
)

unique(sce$condition)

plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "AD")))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0, 2.5)) +
  scale_color_manual(values = c("Control" = "blue", "AD" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


## All genes

# Calculate mean expression across ALL genes
expr_matrix <- assay(sce, "logcounts")
mean_expr_all_genes <- colMeans(expr_matrix)

# Create metadata with mean expression
metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$condition,
  Mean_Expression = mean_expr_all_genes
)

# Filter and prepare data
plot_data_mean <- metadata %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Mean_Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "AD")))

# Filter for specific cell types
specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", 
                        "astrocyte", "central nervous system macrophage")

# Create density plot for mean expression across all genes
plot_data_mean %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Mean_Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Mean expression across all genes (logcounts)", limits = c(0, 0.3)) +
  scale_color_manual(values = c("Control" = "blue", "AD" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

sce <- runUMAP(sce)

plotUMAP(sce, 
         colour_by = "celltype") +
  theme_minimal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

unique(sce$condition)

## Average nonzero count per cell

counts_mat <- counts(sce)

# number of non-zero features per cell
nonzero_per_cell <- scater::nexprs(sce, byrow = FALSE, exprs_values = "counts")

# total counts per cell
totals_per_cell <- colSums(counts_mat)

# average non-zero count per cell = total / number of expressed features
avg_nonzero_per_cell <- totals_per_cell / nonzero_per_cell

summary(avg_nonzero_per_cell)

## Sparsity

# compute sparsity per cell
zero_rate_per_cell <- colSums(counts_mat == 0) / nrow(counts_mat)
summary(zero_rate_per_cell)



(sce <- prepSCE(sce, 
                kid = "celltype", # subpopulation assignments
                gid = "condition",  # group IDs (ctrl/stim)
                sid = "sample_id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

levels(sce$group_id)

sce$group_id <- relevel(sce$group_id, ref = "Control")

## Pseudobulk aggregate

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))


library(limma)
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))

colnames(mm) <- make.names(colnames(mm))
levels(mm) <- make.names(levels(mm))
contrast <- makeContrasts("AD-Control", levels = mm)


# run DS analysis
res <- pbDS(pb, method = 'edgeR', design = mm, contrast = contrast)


# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[2])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)


p1 <- ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size = 7) +   # outlined points
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(title = "DLPFC (BA9)",
       x = NULL,
       y = "Log2 Fold Change (95% CI)") +
  scale_color_manual(values = c(
    "SIRT1" = "#FF6F91",  
    "SIRT2" = "#6FAFFF",  
    "SIRT3" = "#6FFF8F",  
    "SIRT4" = "#C66FFF",  
    "SIRT5" = "#FFB36F",  
    "SIRT6" = "#FF8F6F",  
    "SIRT7" = "#6FFFEF"   
  )) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.y  = element_text(size = 20, colour = "black"),
    axis.text.x  = element_blank(),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    #legend.position = "none",
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio   = 4/1
  )

legend <- cowplot::get_legend(p1)


cowplot::save_plot("/tank/projects/ekashuk/AD/scRNAseq/legend.png", legend, base_height = 4, base_width = 6)


# Anderson et al. 2025 ----

metadata <- read.csv('/tank/projects/public_data/Anderson_et_al_2025/GSE214979_cell_metadata.csv')

rownames(metadata) <- metadata[[1]]  
metadata <- metadata[, -1]  

sce <- readRDS('/tank/projects/public_data/Anderson_et_al_2025/GSE214979_sce.rds')

rownames(sce@colData)
metadata_ordered <- metadata[rownames(colData(sce)), ]

colData(sce) <- DataFrame(metadata_ordered)

unique(sce$Diagnosis)


saveRDS(sce, '/tank/projects/public_data/Anderson_et_al_2025/GSE214979_sce.rds')

### Pseudobulk ----

sce <- readRDS('/tank/projects/public_data/Anderson_et_al_2025/GSE214979_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

sce <- sce[, !is.na(sce$Diagnosis)]

unique(sce$id)


DE = run_de(sce, cell_type_col = "celltype", label_col = "Diagnosis", replicate_col = "id",
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)

DE_limma = run_de(sce, cell_type_col = "celltype", label_col = "Diagnosis", replicate_col = "id",
            de_family = 'pseudobulk', de_method = 'limma', de_type = "trend", n_threads = 5)

de_sig_limma <- subset(DE_limma, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/Anderson_wilcox_test_DE.csv')

DE <- read.csv('/tank/projects/ekashuk/AD/scRNAseq/Anderson_wilcox_test_DE.csv', row.names = 1)
de_sig <- subset(DE, p_val_adj < 0.05)

sirt <- subset(DE, grepl('SIRT', gene))
sirt

library(tidyverse)
library(ggplot2)
sirt

# prepare
df <- sirt %>%
  mutate(
    p_adj = p_val_adj,
    negLog10Padj = -log10(p_adj + 1e-300),                   # avoid Inf
    sig = p_adj < 0.05,
    pct = pmax(Alzheimer.s.pct, Unaffected.pct)                         # size = max prevalence
  ) %>%
  arrange(cell_type)

# bubble plot
ggplot(df, aes(x = cell_type, y = gene)) +
  geom_point(
    aes(
      size = pct,
      fill = avg_logFC,
      color = sig
    ),
    shape = 21,
    stroke = 0.8
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "Percent expressed") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3.5, 3.5),
    name = "avg_logFC"
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    name = "Significant"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 18, colour = "black"),
    axis.text.x  = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    #aspect.ratio   = 3/1
  ) +
  labs(
    x = "",
    y = ""
  )


table(sce$celltype)

unique(sce$Diagnosis)

sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$Diagnosis    
)


plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Unaffected", "Alzheimer's")))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0, 2.5)) +
  scale_color_manual(values = c("Unaffected" = "blue", "Alzheimer's" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


## All genes

# Calculate mean expression across all genes
expr_matrix <- assay(sce, "logcounts")
mean_expr_all_genes <- colMeans(expr_matrix)

# Create metadata with mean expression
metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$Diagnosis,
  Mean_Expression = mean_expr_all_genes
)

# Filter and prepare data
plot_data_mean <- metadata %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Mean_Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Unaffected", "Alzheimer's")))

# Filter for specific cell types
specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", 
                        "astrocyte", "central nervous system macrophage")

# Create density plot for mean expression across all genes
plot_data_mean %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Mean_Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Mean expression across all genes (logcounts)", limits = c(0, 0.3)) +
  scale_color_manual(values = c("Unaffected" = "blue", "Alzheimer's" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

## Peak identification

library(mclust)

# Function to identify high/low expression cells for each cell type and condition
identify_peak_cells <- function(cell_type_data) {
  # Fit Gaussian mixture model with 2 components (assuming bimodal distribution)
  gmm <- Mclust(cell_type_data$Mean_Expression, G = 2)
  
  # Assign cells to peaks based on which Gaussian they belong to
  cell_type_data$Peak <- ifelse(gmm$classification == 1, "Low_Expression", "High_Expression")
  cell_type_data$Peak_Probability <- apply(gmm$z, 1, max)  # Probability of assignment
  
  return(cell_type_data)
}

identify_peak_cells <- function(cell_type_data, prob_threshold = 0.8) {
  # Fit Gaussian mixture model with 2 components
  gmm <- Mclust(cell_type_data$Mean_Expression, G = 2)
  
  # Assign cells to peaks
  cell_type_data$Peak <- ifelse(gmm$classification == 1, "Low_Expression", "High_Expression")
  cell_type_data$Peak_Probability <- apply(gmm$z, 1, max)
  
  # Mark cells that don't meet probability threshold as "Uncertain"
  cell_type_data$Peak_Strict <- ifelse(
    cell_type_data$Peak_Probability >= prob_threshold,
    cell_type_data$Peak,
    "Uncertain"
  )
  
  return(cell_type_data)
}

# Apply to each cell type and condition combination
peak_annotated_data <- plot_data_mean %>%
  filter(CellType %in% specific_celltypes) %>%
  group_by(CellType, Condition) %>%
  group_modify(~ identify_peak_cells(.x)) %>%
  ungroup()


# Create a new column in the metadata
peak_info <- peak_annotated_data %>%
  dplyr::(CellID, CellType, Condition, Mean_Expression, Peak, Peak_Probability)

# Add to colData of sce
colData(sce)$Expression_Group <- NA
colData(sce)$Expression_Group[match(peak_info$CellID, colnames(sce))] <- peak_info$Peak

colData(sce)$Peak_Probability <- NA
colData(sce)$Peak_Probability[match(peak_info$CellID, colnames(sce))] <- peak_info$Peak_Probability

# Create a combined grouping variable for plotting
colData(sce)$Group_Plot <- ifelse(
  !is.na(colData(sce)$Expression_Group),
  paste(colData(sce)$Diagnosis, colData(sce)$Expression_Group, sep = "_"),
  NA
)

sce <- runUMAP(sce)

colData(sce)$Condition_Expression <- ifelse(
  !is.na(colData(sce)$Expression_Group),
  paste(colData(sce)$Diagnosis, colData(sce)$Expression_Group, sep = "_"),
  NA
)

high_expr_cells <- sce[, 
                       colData(sce)$Condition_Expression %in% c("Unaffected_High_Expression", "Alzheimer's_High_Expression")
]

umap_coords <- reducedDim(sce, "UMAP")
colnames(umap_coords) <- c("UMAP1", "UMAP2")

plot_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Condition_Expression = colData(sce)$Condition_Expression
)

ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = subset(plot_df, 
                           !Condition_Expression %in% c("Unaffected_High_Expression", "Alzheimer's_High_Expression") |
                             is.na(Condition_Expression)),
             color = "grey80", alpha = 0.5, size = 0.5) +
  geom_point(data = subset(plot_df, 
                           Condition_Expression %in% c("Unaffected_High_Expression", "Alzheimer's_High_Expression")),
             aes(color = Condition_Expression), 
             alpha = 0.8, size = 1) +
  scale_color_manual(
    values = c(
      "Unaffected_High_Expression" = "blue",
      "Alzheimer's_High_Expression" = "red"
    ),
    name = ""
  ) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size = 15, colour = 'black'),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.position = "right",
    title = element_text(size = 15),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )


assays(sce)

sce@assays@data$logcounts


unique(sce$Diagnosis)

plotUMAP(sce, 
         colour_by = "Diagnosis") +
  theme_minimal() +
  labs(
       x = "UMAP 1",
       y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


## QC

sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

## Average nonzero count per cell

counts_mat <- counts(sce)

# number of non-zero features per cell
nonzero_per_cell <- scater::nexprs(sce, byrow = FALSE, exprs_values = "counts")

# total counts per cell
totals_per_cell <- colSums(counts_mat)

# average non-zero count per cell = total / number of expressed features
avg_nonzero_per_cell <- totals_per_cell / nonzero_per_cell

summary(avg_nonzero_per_cell)

## Sparsity

# compute sparsity per cell
zero_rate_per_cell <- colSums(counts_mat == 0) / nrow(counts_mat)
summary(zero_rate_per_cell)



(sce <- prepSCE(sce, 
                kid = "celltype", # subpopulation assignments
                gid = "Status",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

levels(sce$group_id)

sce$group_id <- relevel(sce$group_id, ref = "Ctrl")

## Pseudobulk aggregate

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))



# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))


library(limma)
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("AD-Ctrl", levels = mm)

# run DS analysis
res_DS <- pbDS(pb, design = mm, contrast = contrast, method='edgeR')

tbl <- res_DS$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[6])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)


ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size = 7) +   # outlined points
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(title = "DLPFC",
       x = NULL,
       y = "Log2 Fold Change (95% CI)") +
  scale_color_manual(values = c(
    "SIRT1" = "#FF6F91",  
    "SIRT2" = "#6FAFFF",  
    "SIRT3" = "#6FFF8F",  
    "SIRT4" = "#C66FFF",  
    "SIRT5" = "#FFB36F",  
    "SIRT6" = "#FF8F6F", 
    "SIRT7" = "#6FFFEF"  
  )) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.y  = element_text(size = 20, colour = "black"),
    axis.text.x  = element_blank(),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio   = 4/1
  )

# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))


# filter FDR < 5%, abs(logFC) > 0.58 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.58)
  dplyr::arrange(u, p_adj.loc)
})

# nb. of DS genes & % of total by cluster
n_de <- vapply(tbl_fil, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

# view top hits in each cluster
top <- bind_rows(lapply(tbl_fil, top_n, 10, p_adj.loc))
format(top[, -ncol(top)], digits = 2)


## Cell-level viz.: Violin plots
plotExpression(sce[, sce$cluster_id == "neuron"],
               features = c('SIRT1', 'SIRT2', 'SIRT3', 'SIRT4', 'SIRT5', 'SIRT6', 'SIRT7'),
               x = "sample_id", colour_by = "group_id", ncol = 3) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotExpression(sce[, sce$cluster_id == "neuron"],
               features = c('SIRT1', 'SIRT2', 'SIRT3', 'SIRT4', 'SIRT5', 'SIRT6', 'SIRT7'),
               x = "group_id", colour_by = "group_id", ncol = 3) +  # Changed x to "group_id"
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Calculating expression frequencies

frq <- calcExprFreqs(sce, assay = "counts", th = 0)
# one sheet per cluster
assayNames(frq)

gids <- levels(sce$group_id)
frq10 <- vapply(as.list(assays(frq)), 
                function(u) apply(u[, gids] > 0.1, 1, any), 
                logical(nrow(sce)))
t(head(frq10))

tbl_fil2 <- lapply(tbl_fil2, function(x) {
  if (is.null(x)) {
    data.frame(gene = character(), stringsAsFactors = FALSE) 
  } else {
    x
  }
})

n_de <- vapply(tbl_fil2, nrow, numeric(1))
p_de <- format(n_de / nrow(sce) * 100, digits = 3)
data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)

resDS(sce, res, bind = "row", frq = frq)

library(UpSetR)
de_gs_by_k <- map(tbl_fil, "gene")
upset(fromList(de_gs_by_k))



res <- as.data.frame(tbl[6])

head(res)

# Morabito et al 2021 ----

sce <- readRDS('/tank/projects/public_data/Morabito_et_al_2021/GSE174367_sce.rds')


metadata <- read.csv('/tank/projects/public_data/Morabito_et_al_2021/GSE174367_snRNA-seq_cell_meta.csv')

sce$name

rownames(metadata) <- metadata[[1]]  
metadata <- metadata[, -1]  

rownames(sce@colData)
metadata_ordered <- metadata[rownames(colData(sce)), ]

colData(sce) <- DataFrame(metadata_ordered)

unique(sce$Cell.Type)

saveRDS(sce, '/tank/projects/public_data/Morabito_et_al_2021/GSE174367_sce.rds')

### Pseudobulk ----
sce <- readRDS('/tank/projects/public_data/Morabito_et_al_2021/GSE174367_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

sce <- sce[, !is.na(sce$Diagnosis)]

unique(sce$SampleID)

colData(sce)$replicate <- as.factor(sce$SampleID)


DE = run_de(sce, cell_type_col = "celltype", label_col = "Diagnosis", 
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)

de_sig <- subset(DE, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/Morabito_wilcox_test_DE.csv')

DE <- read.csv('/tank/projects/ekashuk/AD/scRNAseq/Morabito_wilcox_test_DE.csv', row.names = 1)
de_sig <- subset(DE, p_val_adj < 0.05)

sirt <- subset(DE, grepl('SIRT', gene))
sirt

library(tidyverse)
library(ggplot2)
sirt

# prepare
df <- sirt %>%
  mutate(
    p_adj = p_val_adj,
    negLog10Padj = -log10(p_adj + 1e-300),                   # avoid Inf
    sig = p_adj < 0.05,
    pct = pmax(AD.pct, Control.pct)                         # size = max prevalence
  ) %>%
  arrange(cell_type)

# bubble plot
ggplot(df, aes(x = cell_type, y = gene)) +
  geom_point(
    aes(
      size = pct,
      fill = avg_logFC,
      color = sig
    ),
    shape = 21,
    stroke = 0.8
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "Percent expressed") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3.5, 3.5),
    name = "avg_logFC"
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    name = "Significant"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 18, colour = "black"),
    axis.text.x  = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    #aspect.ratio   = 3/1
  ) +
  labs(
    x = "",
    y = ""
  )



table(sce$celltype)

sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$Diagnosis    
)

unique(sce$Diagnosis)
plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "AD")))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0, 2.5)) +
  scale_color_manual(values = c("Control" = "blue", "AD" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


## All genes

# Calculate mean expression across ALL genes
expr_matrix <- assay(sce, "logcounts")
mean_expr_all_genes <- colMeans(expr_matrix)

# Create metadata with mean expression
metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$Diagnosis,
  Mean_Expression = mean_expr_all_genes
)

# Filter and prepare data
plot_data_mean <- metadata %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Mean_Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "AD")))

# Filter for specific cell types
specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", 
                        "astrocyte", "central nervous system macrophage")

# Create density plot for mean expression across all genes
plot_data_mean %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Mean_Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Mean expression across all genes (logcounts)", limits = c(0, 0.3)) +
  scale_color_manual(values = c("Control" = "blue", "AD" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


sce <- runUMAP(sce)

plotUMAP(sce, 
         colour_by = "celltype") +
  theme_minimal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

unique(sce$Age)

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

sce <- sce[, !is.na(sce$Diagnosis)]

unique(sce$Diagnosis)

## Average nonzero count per cell

counts_mat <- counts(sce)

# number of non-zero features per cell
nonzero_per_cell <- scater::nexprs(sce, byrow = FALSE, exprs_values = "counts")

# total counts per cell
totals_per_cell <- colSums(counts_mat)

# average non-zero count per cell = total / number of expressed features
avg_nonzero_per_cell <- totals_per_cell / nonzero_per_cell

summary(avg_nonzero_per_cell)

## Sparsity

# compute sparsity per cell
zero_rate_per_cell <- colSums(counts_mat == 0) / nrow(counts_mat)
summary(zero_rate_per_cell)


(sce <- prepSCE(sce, 
                kid = "celltype", # subpopulation assignments
                gid = "Diagnosis",  # group IDs (ctrl/stim)
                sid = "SampleID",   # sample IDs (ctrl/stim.1234)
                drop = F))  # drop all other colData columns

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

levels(sce$group_id)

#sce$group_id <- relevel(sce$group_id, ref = "79")

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

## Pseudobulk aggregate

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))

# run DS analysis
#res <- pbDS(pb, method='DESeq2')

sce$group_id
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("AD-Control", levels = mm)

# run DS analysis
res <- pbDS(pb, method='edgeR', design = mm, contrast = contrast)

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[5])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size = 7) +   # outlined points
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(title = "PFC",
       x = NULL,
       y = "Log2 Fold Change (95% CI)") +
  scale_color_manual(values = c(
    "SIRT1" = "#FF6F91",  
    "SIRT2" = "#6FAFFF",  
    "SIRT3" = "#6FFF8F",  
    "SIRT4" = "#C66FFF",  
    "SIRT5" = "#FFB36F",  
    "SIRT6" = "#FF8F6F",  
    "SIRT7" = "#6FFFEF"   
  )) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.y  = element_text(size = 20, colour = "black"),
    axis.text.x  = element_blank(),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio   = 4/1
  )

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)


res.df <- as.data.frame(tbl[7])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$lfcSE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$lfcSE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size=5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = 'OPC',
       x = "",
       y = "Log2FoldChange with 95% CI") +
  theme_minimal() +
  scale_color_manual(values = c("SIRT1" = "#E41A1C", 
                                "SIRT2" = "#377EB8",
                                "SIRT3" = "#4DAF4A", 
                                "SIRT4" = "#984EA3",
                                "SIRT5" = "#FF7F00", 
                                "SIRT6" = "coral2",
                                "SIRT7" = "#A65628")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    aspect.ratio = 1/3,
    legend.position = "none",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


# Gerrits et al. 2021 ----

### Pseudobulk ----
sce <- readRDS('/tank/projects/public_data/Gerrits_et_al_2021/Gerrits_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Create a frequency table
freq_table <- as.data.frame(table(
  CellType = sce$celltype,
  Condition = sce$donor_status
))

# Find cell types with only 1 cell in any condition
problem_clusters <- freq_table %>%
  filter(Freq == 1) %>%
  pull(CellType) %>%
  unique()

print("Clusters with only 1 cell in any condition:")
print(problem_clusters)

if (length(problem_clusters) > 0) {
  sce_clean <- sce[, !sce$celltype %in% problem_clusters]
  
  sce_clean$celltype <- droplevels(factor(sce_clean$celltype))
  
  cat(sprintf("Removed %d problematic clusters.\n", length(problem_clusters)))
  cat("Remaining clusters:", paste(unique(sce_clean$celltype), collapse = ", "), "\n")
}

sce_oc <- sce[, sce$subregion == 'occipital cortex']
sce_oct <- sce[, sce$subregion == 'occipitotemporal cortex']

unique(sce$sample)

DE = run_de(sce_oct, cell_type_col = "celltype", label_col = "donor_status", replicate_col = "sample",
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)

de_sig <- subset(DE, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/Gerrits_oct_wilcox_test_DE.csv')

DE <- read.csv('/tank/projects/ekashuk/AD/scRNAseq/Gerrits_oct_wilcox_test_DE.csv', row.names = 1)
de_sig <- subset(DE, p_val_adj < 0.05)

sirt <- subset(DE, grepl('SIRT', gene))
sirt

library(tidyverse)
library(ggplot2)
colnames(sirt)

# prepare
df <- sirt %>%
  mutate(
    p_adj = p_val_adj,
    negLog10Padj = -log10(p_adj + 1e-300),                   # avoid Inf
    sig = p_adj < 0.05,
    pct = pmax(Alzheimer.s.disease.pct, Healthy.pct)                         # size = max prevalence
  ) %>%
  arrange(cell_type)

# bubble plot
ggplot(df, aes(x = cell_type, y = gene)) +
  geom_point(
    aes(
      size = pct,
      fill = avg_logFC,
      color = sig
    ),
    shape = 21,
    stroke = 0.8
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "Percent expressed") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3.5, 3.5),
    name = "avg_logFC"
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    name = "Significant"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 18, colour = "black"),
    axis.text.x  = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    #aspect.ratio   = 3/1
  ) +
  labs(
    x = "",
    y = ""
  )

table(sce_oc$cell_type)

sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce_oc)]

expr_matrix <- assay(sce_oct, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

unique(sce$donor_status)

metadata <- data.frame(
  CellID = colnames(sce_oct),
  CellType = sce_oct$celltype,        
  Condition = sce_oct$donor_status   
)


plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Healthy", "Alzheimer’s disease")))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0, 2.5)) +
  scale_color_manual(values = c("Healthy" = "blue", "Alzheimer’s disease" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

sce <- runUMAP(sce)

unique(sce$celltype)

plotUMAP(sce, 
         colour_by = "cell_type") +
  theme_minimal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

## Distributions

sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce_oc)]

expr_matrix <- assay(sce_oct, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce_oc),
  CellType = sce_oc$celltype,        
  Condition = sce_oc$donor_status   
)

unique(sce_oc$donor_status)
plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("Healthy", "Alzheimer’s disease")))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression)) +
  geom_density(aes(y = after_stat(scaled), color = Gene), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0, 2.5)) +
  labs(
    color = "Cell type") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


## All genes

# Calculate mean expression across ALL genes
expr_matrix <- assay(sce_oc, "logcounts")
mean_expr_all_genes <- colMeans(expr_matrix)

# Create metadata with mean expression
metadata <- data.frame(
  CellID = colnames(sce_oc),
  CellType = sce_oc$celltype,        
  Condition = sce_oc$donor_status,
  Mean_Expression = mean_expr_all_genes
)

unique(metadata$Condition)

# Filter and prepare data
plot_data_mean <- metadata %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Mean_Expression)) %>%
  mutate(Condition = factor(Condition))

# Filter for specific cell types
specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", 
                        "astrocyte", "central nervous system macrophage")

unique(plot_data_mean$Condition)

# Create density plot for mean expression across all genes

plot_data_mean %>%
  filter(CellType %in% specific_celltypes, !is.na(Condition)) %>%
  ggplot(aes(x = Mean_Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Mean expression across all genes (logcounts)", limits = c(0, 0.3)) +
  scale_color_manual(values = c("Healthy" = "blue", "Alzheimer’s disease" = "red")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )






unique(sce$subregion)

sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

sce_oc <- sce[, sce$subregion == 'occipital cortex']
sce_oct <- sce[, sce$subregion == 'occipitotemporal cortex']

## Average nonzero count per cell

counts_mat <- counts(sce_oct)

# number of non-zero features per cell
nonzero_per_cell <- scater::nexprs(sce_oct, byrow = FALSE, exprs_values = "counts")

# total counts per cell
totals_per_cell <- colSums(counts_mat)

# average non-zero count per cell = total / number of expressed features
avg_nonzero_per_cell <- totals_per_cell / nonzero_per_cell

summary(avg_nonzero_per_cell)

## Sparsity

# compute sparsity per cell
zero_rate_per_cell <- colSums(counts_mat == 0) / nrow(counts_mat)
summary(zero_rate_per_cell)



(sce_oct <- prepSCE(sce_oct, 
                kid = "celltype", # subpopulation assignments
                gid = "donor_status",  # group IDs (ctrl/stim)
                sid = "sample_ID",   # sample IDs (ctrl/stim.1234)
                drop = F))  # drop all other colData columns

nk <- length(kids <- levels(sce_oct$cluster_id))
ns <- length(sids <- levels(sce_oct$sample_id))
names(kids) <- kids; names(sids) <- sids

levels(sce_oct$group_id)

#sce$group_id <- relevel(sce$group_id, ref = "79")

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

## Pseudobulk aggregate

pb <- aggregateData(sce_oct,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))

# run DS analysis
#res <- pbDS(pb, method='DESeq2')

library(limma)

ei <- metadata(sce_oct)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))


colnames(mm) <- make.names(colnames(mm))

colnames(mm)
contrast <- makeContrasts("Alzheimer.s.disease - Healthy", levels = mm)

# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast)


# run DS analysis
#res <- pbDS(pb, method='DESeq2')

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[8])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

unique(sce$subregion)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size = 7) +   
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(title = "Occipitotemporal cortex",
       x = NULL,
       y = "Log2 Fold Change (95% CI)") +
  scale_color_manual(values = c(
    "SIRT1" = "#FF6F91",  
    "SIRT2" = "#6FAFFF",  
    "SIRT3" = "#6FFF8F", 
    "SIRT4" = "#C66FFF",  
    "SIRT5" = "#FFB36F",  
    "SIRT6" = "#FF8F6F", 
    "SIRT7" = "#6FFFEF" 
  )) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.y  = element_text(size = 20, colour = "black"),
    axis.text.x  = element_blank(),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio   = 4/1
  )


# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)


res.df <- as.data.frame(tbl[7])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$lfcSE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$lfcSE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size=5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = 'OPC',
       x = "",
       y = "Log2FoldChange with 95% CI") +
  theme_minimal() +
  scale_color_manual(values = c("SIRT1" = "#E41A1C", 
                                "SIRT2" = "#377EB8",
                                "SIRT3" = "#4DAF4A", 
                                "SIRT4" = "#984EA3",
                                "SIRT5" = "#FF7F00", 
                                "SIRT6" = "coral2",
                                "SIRT7" = "#A65628")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    aspect.ratio = 1/3,
    legend.position = "none",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


# Jakel et al. 2019 ----

sce <- readRDS('/tank/projects/public_data/Otero-Garcia_et_al_2020/GSE129308_sce_ann.rds')

unique(sce$condition)
table(sce$celltype)

assays(sce)

sce <- runUMAP(sce)

plotUMAP(sce, 
         colour_by = "celltype") +
  theme_minimal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

unique(sce$condition)

(sce <- prepSCE(sce, 
                kid = "celltype", # subpopulation assignments
                gid = "condition",  # group IDs (ctrl/stim)
                sid = "sample_id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

levels(sce$group_id)

sce$group_id <- relevel(sce$group_id, ref = "Control")

## Pseudobulk aggregate

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))


library(limma)
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))

colnames(mm) <- make.names(colnames(mm))
levels(mm) <- make.names(levels(mm))
contrast <- makeContrasts("AD-Control", levels = mm)


# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast)


# run DS analysis
#res <- pbDS(pb, method='DESeq2')

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[2])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size=5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(-2, 2) +
  coord_flip() +
  labs(title = 'Neurons',
       x = "",
       y = "Log2FoldChange with 95% CI") +
  theme_minimal() +
  scale_color_manual(values = c("SIRT1" = "#E41A1C", 
                                "SIRT2" = "#377EB8",
                                "SIRT3" = "#4DAF4A", 
                                "SIRT4" = "#984EA3",
                                "SIRT5" = "#FF7F00", 
                                "SIRT6" = "coral2",
                                "SIRT7" = "#A65628")) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    aspect.ratio = 1/3,
    legend.position = "none",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


h5ad_file <- '/tank/projects/public_data/Sun_et_al_2022/human_brain_CV_SunN_2022_10x/processedData/Sun_2022.h5ad'
convertFormat(h5ad_file, from="anndata", to="seurat",
              outFile='/tank/projects/public_data/Sun_et_al_2022/human_brain_CV_SunN_2022_10x/processedData/Sun_2022.rds')

# Leng et al.  ----

sce <- readRDS('/tank/projects/public_data/Leng_et_al_2021/EntorhinalCortex_Leng_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Create a frequency table
freq_table <- as.data.frame(table(
  CellType = sce$celltype,
  Condition = sce$BraakStage
))

# Find cell types with only 1 cell in any condition
problem_clusters <- freq_table %>%
  dplyr::filter(Freq == 1) %>%
  pull(CellType) %>%
  unique()

print("Clusters with only 1 cell in any condition:")
print(problem_clusters)

if (length(problem_clusters) > 0) {
  sce_clean <- sce[, !sce$celltype %in% problem_clusters]
  
  sce_clean$celltype <- droplevels(factor(sce_clean$celltype))
  
  cat(sprintf("Removed %d problematic clusters.\n", length(problem_clusters)))
  cat("Remaining clusters:", paste(unique(sce_clean$celltype), collapse = ", "), "\n")
}

unique(sce$BraakStage)

sce_stage_0_3 <- sce[, sce$BraakStage %in% c(0, 2)]

sce_stage_0_3$BraakStage <- factor(sce_stage_0_3$BraakStage, levels = c(2, 0))


DE = run_de(sce_stage_0_3, cell_type_col = "celltype", label_col = "BraakStage", replicate_col = "SampleID",
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)

de_sig <- subset(DE, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/wilcox_res/Leng_ec_0_2_wilcox_test_DE.csv')

DE <- read.csv('/tank/projects/ekashuk/AD/scRNAseq/Leng_sfg_0_6_wilcox_test_DE.csv', row.names = 1)
de_sig <- subset(DE, p_val_adj < 0.05)

sirt <- subset(DE, grepl('SIRT', gene))

library(tidyverse)
library(ggplot2)
colnames(sirt)

# prepare
df <- sirt %>%
  mutate(
    p_adj = p_val_adj,
    negLog10Padj = -log10(p_adj + 1e-300),                   # avoid Inf
    sig = p_adj < 0.05,
    pct = pmax(X6.pct, X0.pct)                         # size = max prevalence
  ) %>%
  arrange(cell_type)

# bubble plot
ggplot(df, aes(x = cell_type, y = gene)) +
  geom_point(
    aes(
      size = pct,
      fill = avg_logFC,
      color = sig
    ),
    shape = 21,
    stroke = 0.8
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "Percent expressed") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-5, 5),
    name = "avg_logFC"
  ) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    name = "Significant"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 18, colour = "black"),
    axis.text.x  = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    #aspect.ratio   = 3/1
  ) +
  labs(
    x = "",
    y = ""
  )


table(sce$celltype)
table(sce$BraakStage)

sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$BraakStage    
)


plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression)) %>%
  mutate(Condition = factor(Condition, levels = c("0", "2", "6")))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

p <- plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression)) +
  geom_density(aes(y = after_stat(scaled), color = CellType), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ Condition, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0, 6)) +
  labs(
    color = "Cell type") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

## Distributions 
sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")

available_sirt_genes <- sirt_genes[sirt_genes %in% rownames(sce)]

expr_matrix <- assay(sce, "logcounts")
sirt_expr <- as.matrix(expr_matrix[available_sirt_genes, ])

metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$BraakStage    
)


plot_data <- as.data.frame(t(sirt_expr))
plot_data$CellID <- rownames(plot_data)
plot_data <- plot_data %>%
  left_join(metadata, by = "CellID") %>%
  pivot_longer(
    cols = all_of(available_sirt_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

unique(plot_data$Condition)

plot_data <- plot_data %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Expression))

specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", "astrocyte", "central nervous system macrophage")  

plot_data %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_grid(Gene ~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Expression level (logcounts)", limits = c(0, 2.5)) +
  scale_color_manual(values = c("0" = "blue", "2" = "red", "6"="purple")) +
  labs(color = "Braak Stage") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


## All genes

# Calculate mean expression across ALL genes
expr_matrix <- assay(sce, "logcounts")
mean_expr_all_genes <- colMeans(expr_matrix)

# Create metadata with mean expression
metadata <- data.frame(
  CellID = colnames(sce),
  CellType = sce$celltype,        
  Condition = sce$BraakStage,
  Mean_Expression = mean_expr_all_genes
)

# Filter and prepare data
plot_data_mean <- metadata %>%
  filter(!is.na(CellType), !is.na(Condition), !is.na(Mean_Expression)) 

# Filter for specific cell types
specific_celltypes <- c("neuron", "oligodendrocyte", "oligodendrocyte precursor cell", 
                        "astrocyte", "central nervous system macrophage")

# Create density plot for mean expression across all genes
plot_data_mean %>%
  filter(CellType %in% specific_celltypes) %>%
  ggplot(aes(x = Mean_Expression, color = Condition)) +
  geom_density(aes(y = after_stat(scaled)), 
               linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ CellType, scales = "free") +
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage of cells") +
  scale_x_continuous(name = "Mean expression across all genes (logcounts)", limits = c(0, 0.3)) +
  scale_color_manual(values = c("0" = "blue", "2" = "red", "6"="purple"))+
  labs(color = "Braak Stage") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

assays(sce)

sce <- runUMAP(sce)

plotUMAP(sce, 
         colour_by = "celltype") +
  theme_minimal() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 15, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(size=15, colour='black'),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.position = "right",
    title = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


sce <- sce[rowSums(counts(sce) > 0) > 0, ]

sce$BraakStage <- factor(sce$BraakStage, 
                       levels = c(0, 2, 6), 
                       labels = c("Stage0", "Stage2", "Stage6"))

unique(sce$BraakStage)


dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

## Average nonzero count per cell

counts_mat <- counts(sce)

# number of non-zero features per cell
nonzero_per_cell <- scater::nexprs(sce, byrow = FALSE, exprs_values = "counts")

# total counts per cell
totals_per_cell <- colSums(counts_mat)

# average non-zero count per cell = total / number of expressed features
avg_nonzero_per_cell <- totals_per_cell / nonzero_per_cell

summary(avg_nonzero_per_cell)

## Sparsity

# compute sparsity per cell
zero_rate_per_cell <- colSums(counts_mat == 0) / nrow(counts_mat)
summary(zero_rate_per_cell)


(sce <- prepSCE(sce, 
                kid = "celltype", # subpopulation assignments
                gid = "BraakStage",  # group IDs (ctrl/stim)
                sid = "SampleID",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce$cluster_id, sce$sample_id))

levels(sce$group_id)

#sce$group_id <- relevel(sce$group_id, ref = "Control")

## Pseudobulk aggregate

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

# pseudobulks for 1st subpopulation
t(head(assay(pb)))

## Pseudobulk-level MDS plot
(pb_mds <- pbMDS(pb))

library(limma)
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))

colnames(mm) <- make.names(colnames(mm))
levels(mm) <- make.names(levels(mm))
contrast <- makeContrasts("Stage6-Stage0", levels = mm)


# run DS analysis
res <- pbDS(pb, method = 'edgeR', design = mm, contrast = contrast)

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

res.df <- as.data.frame(tbl[3])
colnames(res.df) <- gsub("^.*\\.", "", colnames(res.df))

head(res.df)

sirt <- subset(res.df, grepl('SIRT', gene))

sirt$SE <- abs(sirt$logFC) / qnorm(1 - sirt$p_val/2)

# Calculate 95% confidence intervals
sirt$CI_lower <- sirt$logFC - 1.96 * sirt$SE
sirt$CI_upper <- sirt$logFC + 1.96 * sirt$SE

sirt$padj <- p.adjust(sirt$p_val)

library(ggplot2)

ggplot(sirt, aes(x = gene, y = logFC, color = gene)) +  
  geom_point(size = 7) +   # outlined points
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  scale_y_continuous(limits = c(-2, 2.8), expand = c(0, 0)) +
  labs(title = "DLPFC (BA9)",
       x = NULL,
       y = "Log2 Fold Change (95% CI)") +
  scale_color_manual(values = c(
    "SIRT1" = "#FF6F91",  
    "SIRT2" = "#6FAFFF",  
    "SIRT3" = "#6FFF8F",  
    "SIRT4" = "#C66FFF",
    "SIRT5" = "#FFB36F", 
    "SIRT6" = "#FF8F6F",  
    "SIRT7" = "#6FFFEF"   
  )) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.y  = element_text(size = 20, colour = "black"),
    axis.text.x  = element_blank(),
    axis.title.x = element_text(size = 15, face = "bold"),
    title        = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.border   = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio   = 4/1
  )

# Sun et al. ----

sce <- readRDS('/tank/projects/public_data/Sun_et_al_2022/Sun_2022_sce_ann.rds')

table(sce$donor_status)
table(sce$brain_region)
table(sce$id)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


#sce <- readRDS('/tank/projects/public_data/Anderson_et_al_2025/GSE214979_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

keep_genes <- rowSums(assay(sce, "counts") > 0) > 0
sce <- sce[keep_genes, ]

sce_ag <- sce[, sce$brain_region == 'Angular_gyrus']
sce_at <- sce[, sce$brain_region == 'Anterior_thalamus']
sce_ec <- sce[, sce$brain_region == 'Entorhinal_cortex']
sce_hip <- sce[, sce$brain_region == 'Hippocampus']
sce_mc <- sce[, sce$brain_region == 'Midtemporal_cortex']
sce_pc <- sce[, sce$brain_region == 'Prefrontal_cortex']

DE = run_de(sce, cell_type_col = "celltype", label_col = "donor_status", replicate_col = "id",
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)


write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/Anderson_wilcox_test_DE.csv')

## Lau et al. ----

sce <- readRDS('/tank/projects/public_data/Lau_et_al_2022/GSE157827_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Create a frequency table
freq_table <- as.data.frame(table(
  CellType = sce$celltype,
  Condition = sce$condition
))

# Find cell types with only 1 cell in any condition
problem_clusters <- freq_table %>%
  dplyr::filter(Freq == 1) %>%
  pull(CellType) %>%
  unique()

print("Clusters with only 1 cell in any condition:")
print(problem_clusters)

if (length(problem_clusters) > 0) {
  sce_clean <- sce[, !sce$celltype %in% problem_clusters]
  
  sce_clean$celltype <- droplevels(factor(sce_clean$celltype))
  
  cat(sprintf("Removed %d problematic clusters.\n", length(problem_clusters)))
  cat("Remaining clusters:", paste(unique(sce_clean$celltype), collapse = ", "), "\n")
}

unique(sce$sample_number)

DE = run_de(sce, cell_type_col = "celltype", label_col = "condition", replicate_col = "sample_number",
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)

de_sig <- subset(DE, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/wilcox_res/Lau_wilcox_test_DE.csv')

## Zhang et al. ----

sce <- readRDS('/tank/projects/public_data/Zhang_et_al_2023/GSE188545_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Create a frequency table
freq_table <- as.data.frame(table(
  CellType = sce$celltype,
  Condition = sce$condition
))

# Find cell types with only 1 cell in any condition
problem_clusters <- freq_table %>%
  dplyr::filter(Freq == 1) %>%
  pull(CellType) %>%
  unique()

print("Clusters with only 1 cell in any condition:")
print(problem_clusters)

if (length(problem_clusters) > 0) {
  sce_clean <- sce[, !sce$celltype %in% problem_clusters]
  
  sce_clean$celltype <- droplevels(factor(sce_clean$celltype))
  
  cat(sprintf("Removed %d problematic clusters.\n", length(problem_clusters)))
  cat("Remaining clusters:", paste(unique(sce_clean$celltype), collapse = ", "), "\n")
}

unique(sce$sample_id)

DE = run_de(sce, cell_type_col = "celltype", label_col = "condition", replicate_col = "sample_id",
            de_family = 'singlecell', de_method = 'wilcox', n_threads = 5)

de_sig <- subset(DE, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/wilcox_res/Zhang_wilcox_test_DE.csv')

## Grubman et al. ----

sce <- readRDS('/tank/projects/public_data/Grubman_et_al_2019/Grubman_2019_sce_ann.rds')

## QC
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

min_cells <- 10
keep_types <- names(which(
  apply(table(sce$celltype, sce$donor_status), 1, min) >= min_cells
))

sce_filtered <- sce[, sce$celltype %in% keep_types]

DE <- run_de(sce_filtered, cell_type_col = "celltype", label_col = "donor_status",
             replicate_col = "sample_ID", de_family = 'singlecell',
             de_method = 'wilcox', n_threads = 5)

de_sig <- subset(DE, p_val_adj < 0.05)

write.csv(DE, '/tank/projects/ekashuk/AD/scRNAseq/wilcox_res/Grubman_wilcox_test_DE.csv')

sce$sample_ID

## ---- Load and prepare data ----

base_path <- "/tank/projects/ekashuk/AD/scRNAseq/wilcox_res"

files <- list.files(base_path, pattern = "*.csv", full.names = TRUE)

all_de <- map_dfr(files, function(f) {
  read_csv(f) %>%
    mutate(experiment = tools::file_path_sans_ext(basename(f)))
})

sirt_data <- all_de %>%
  dplyr::filter(gene == "SIRT6") %>%
  mutate(
    logFDR = -log10(p_val_adj),
    experiment = recode(experiment,
                        "Anderson_wilcox_test_DE"    = "Anderson: Dorsolateral prefrontal cortex",
                        "Morabito_wilcox_test_DE"    = "Morabito: Prefrontal cortex",
                        "Leng_sfg_wilcox_test_DE"    = "Leng: Superior frontal gyrus (Braak 0 vs 2)",
                        "Leng_sfg_0_6_wilcox_test_DE"= "Leng: Superior frontal gyrus (Braak 0 vs 6)",
                        "OteroGarcia_wilcox_test_DE" = "Otero-Garcia: Prefrontal cortex (BA9)",
                        "Gerrits_oc_wilcox_test_DE"  = "Gerrits: Occipital cortex",
                        "Gerrits_oct_wilcox_test_DE" = "Gerrits: Occipitotemporal cortex",
                        "Lau_wilcox_test_DE"         = "Lau: Prefrontal cortex",
                        "Leng_ec_0_2_wilcox_test_DE" = "Leng: Entorhinal cortex (Braak 0 vs 2)",
                        "Leng_ec_0_6_wilcox_test_DE" = "Leng: Entorhinal cortex (Braak 0 vs 6)",
                        "Zhang_wilcox_test_DE"       = "Zhang: Middle temporal gyrus",
                        "Grubman_wilcox_test_DE"     = "Grubman: Entorhinal cortex"
    ),
    experiment = factor(experiment),
    cell_type  = factor(cell_type)
  )

sirt_sig <- sirt_data %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  mutate(cell_type = recode(cell_type,
                            "central nervous system macrophage" = "macrophage",
                            "oligodendrocyte precursor cell"    = "OPC"
  )) %>%
  subset(!cell_type %in% c("pericyte", "fibroblast"))

cell_order <- c("astrocyte", "macrophage", "neuron", "oligodendrocyte", "OPC")

## ---- Panel geometry shared between both plots ----

n_cats <- c(astrocyte = 4, macrophage = 3, neuron = 4, oligodendrocyte = 4, OPC = 3)
unit_w <- 0.75
panel_widths <- unit(n_cats[cell_order] * unit_w, "cm")

fill_limits <- c(-0.81, 0.81)
fill_scale <- scale_fill_gradient2(
  low = "blue", mid = "white", high = "red",
  midpoint = 0, limits = fill_limits,
  name = expression(paste("average log"[2], "FC"))
)

## ---- Top plot: main cross-sectional data by region code ----

main_data <- sirt_sig %>%
  filter(!str_starts(as.character(experiment), "Leng")) %>%
  mutate(
    region_group = case_when(
      str_detect(experiment, "Anderson|Otero|Lau|Morabito") ~ "PFC",
      str_detect(experiment, "Gerrits")                      ~ "OC",
      str_detect(experiment, "Grubman")                      ~ "EC",
      str_detect(experiment, "Zhang")                        ~ "MTG",
      TRUE ~ "OTHER"
    ),
    code = case_when(
      str_detect(experiment, "Anderson") ~ "PFC1",
      str_detect(experiment, "Otero")    ~ "PFC2",
      str_detect(experiment, "Lau")      ~ "PFC3",
      str_detect(experiment, "Morabito") ~ "PFC4",
      str_detect(experiment, "Occipital cortex") & str_detect(experiment, "Gerrits") ~ "OC1",
      str_detect(experiment, "Occipitotemporal") ~ "OC2",
      str_detect(experiment, "Grubman")  ~ "EC",
      str_detect(experiment, "Zhang")    ~ "MTG",
      TRUE ~ "OTHER"
    ),
    code = factor(code, levels = c("PFC1","PFC2","PFC3","PFC4","OC1","OC2","EC","MTG")),
    cell_type = factor(cell_type, levels = cell_order)
  )

p_top <- ggplot(main_data,
                aes(x = code, y = gene, size = logFDR, fill = avg_logFC)) +
  geom_point(shape = 21) +
  fill_scale +
  scale_size_area(max_size = 8, breaks = c(10, 30, 50),
                  name = expression(paste(-log[10], "(FDR)"))) +
  facet_grid(. ~ cell_type, scales = "free_x", space = "free_x") +
  force_panelsizes(rows = unit(1, "cm"), cols = panel_widths) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(x = "", y = "")

## ---- Leng meta-analysis by Braak stage ----

leng_data_full <- sirt_data %>%
  filter(str_starts(as.character(experiment), "Leng")) %>%
  filter(!cell_type %in% c("pericyte", "fibroblast")) %>%
  filter(!is.na(p_val_adj), !is.na(avg_logFC)) %>%
  mutate(
    region = case_when(
      str_detect(experiment, "Superior frontal gyrus") ~ "SFG",
      str_detect(experiment, "Entorhinal cortex")       ~ "EC",
      TRUE ~ NA_character_
    ),
    braak_stage = case_when(
      str_detect(experiment, "0 vs 2") ~ "2 vs 0",
      str_detect(experiment, "0 vs 6") ~ "6 vs 0",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("2 vs 0", "6 vs 0")),
    cell_type_short = recode(as.character(cell_type),
                             "central nervous system macrophage" = "macrophage",
                             "oligodendrocyte precursor cell"    = "OPC",
                             .default = as.character(cell_type)
    )
  )

leng_meta_input <- leng_data_full %>%
  mutate(
    p_safe = pmin(pmax(p_val_adj, 1e-300), 1 - 1e-10),
    z      = qnorm(p_safe / 2, lower.tail = FALSE) * sign(avg_logFC),
    weight = -log10(p_safe)
  )

leng_meta <- leng_meta_input %>%
  group_by(cell_type, braak_stage) %>%
  filter(n() >= 2) %>%
  summarise(
    pooled_logFC = weighted.mean(avg_logFC, w = weight),
    logFC_min    = min(avg_logFC),
    logFC_max    = max(avg_logFC),
    z_combined   = sum(weight * z) / sqrt(sum(weight^2)),   # метод Стоуффера
    pooled_p     = 2 * pnorm(abs(z_combined), lower.tail = FALSE),
    k            = n(),
    .groups = "drop"
  ) %>%
  mutate(
    logFDR = -log10(pooled_p),
    sig_flag = pooled_p < 0.05,
    cell_type_short = recode(as.character(cell_type),
                             "central nervous system macrophage" = "macrophage",
                             "oligodendrocyte precursor cell"    = "OPC",
                             .default = as.character(cell_type)
    )
  )

low_power_types <- c("endothelial cell", "leukocyte")

leng_meta_main <- leng_meta %>%
  filter(!cell_type %in% low_power_types) %>%
  mutate(cell_type_short = factor(cell_type_short, levels = cell_order))

leng_data_main <- leng_data_full %>%
  filter(!cell_type %in% low_power_types) %>%
  mutate(cell_type_short = factor(cell_type_short, levels = cell_order))

leng_meta_supp <- leng_meta %>% filter(cell_type %in% low_power_types)
leng_data_supp <- leng_data_full %>% filter(cell_type %in% low_power_types)

region_shapes <- c(SFG = 24, EC = 25)

p_bottom <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(data = leng_data_main,
             aes(x = braak_stage, y = avg_logFC, shape = region),
             fill = "grey75", color = "grey50",
             size = 2.2, alpha = 0.65, stroke = 0.4,
             position = position_jitter(width = 0.05, height = 0)) +
  geom_line(data = leng_meta_main,
            aes(x = braak_stage, y = pooled_logFC, group = cell_type),
            color = "grey40", linewidth = 0.4) +
  geom_point(data = leng_meta_main,
             aes(x = braak_stage, y = pooled_logFC,
                 size = -log10(pooled_p), fill = pooled_logFC),
             shape = 21, stroke = 0.5) +
  geom_text(data = filter(leng_meta_main, sig_flag),
            aes(x = braak_stage, y = pooled_logFC, label = "*"),
            vjust = -1.6, size = 6, fontface = "bold") +
  scale_shape_manual(values = region_shapes, name = "Region") +
  fill_scale +
  scale_size_area(max_size = 8, breaks = c(0.5, 1.5, 3),
                  name = expression(paste(-log[10], "(pooled p, Stouffer)"))) +
  scale_y_continuous(limits = c(-1, 0.5), breaks = seq(-1, 0.5, 0.5),
                     minor_breaks = seq(-1, 0.5, 0.25)) +
  facet_wrap(~ cell_type_short, nrow = 1) +
  force_panelsizes(rows = unit(2.5, "cm"), cols = panel_widths) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(x = "Braak stage", y = expression(paste("weighted mean log"[2], "FC")),
  )

p_bottom


## ---- Assemble final figure ----

library(gtable)
library(grid)

extract_legend <- function(p) {
  gt <- ggplotGrob(p)
  legend_gt <- gtable_filter(gt, "guide-box", trim = TRUE)
  legend_gt
}

fdr_legend  <- extract_legend(p_top)
fill_legend <- extract_legend(p_bottom)

class(fdr_legend) 

grid.newpage(); grid.draw(fdr_legend)
grid.newpage(); grid.draw(fill_legend)


library(cowplot)

p_top_clean    <- p_top    + theme(legend.position = "none")
p_bottom_clean <- p_bottom + theme(legend.position = "none")

plots_column <- plot_grid(
  p_top_clean, p_bottom_clean,
  ncol = 1, rel_heights = c(1, 2.5),
  align = "v", axis = "lr"
)

legend_column <- plot_grid(
  fdr_legend, fill_legend,
  ncol = 1, rel_heights = c(1, 2.2)
)

final_figure <- plot_grid(
  plots_column, legend_column,
  ncol = 2, rel_widths = c(10, 2)
)

final_figure
ggsave("sirt6_figure_f.png", final_figure, width = 14, height = 9, dpi = 300, bg = "white")

real_legend <- fdr_legend$grobs[[1]]

class(real_legend)   

grid.newpage()
grid.draw(real_legend)

png("AD_legend_top.png", width = 3, height = 4, units = "in", res = 1200, bg = "transparent")
grid.newpage()
grid.draw(real_legend)
dev.off()

## ---- Leng
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)

## ---- 1. Базовый датасет с весами ----

sirt_all <- sirt_data %>%
  filter(!cell_type %in% c("pericyte", "fibroblast")) %>%
  mutate(
    cell_type = recode(as.character(cell_type),
                       "central nervous system macrophage" = "macrophage",
                       "oligodendrocyte precursor cell"    = "OPC",
                       .default = as.character(cell_type)
    ),
    region = case_when(
      str_detect(experiment, "Anderson|Otero")            ~ "PFC/BA9",
      str_detect(experiment, "Lau|Morabito")               ~ "PFC",
      str_detect(experiment, "Occipitotemporal")            ~ "OCT",
      str_detect(experiment, "Occipital cortex")            ~ "OC",
      str_detect(experiment, "Grubman")                     ~ "EC",
      str_detect(experiment, "Zhang")                        ~ "MTG",
      str_detect(experiment, "Superior frontal gyrus")       ~ "SFG",
      str_detect(experiment, "Leng.*Entorhinal")              ~ "EC",
      TRUE ~ NA_character_
    ),
    sig_flag = p_val_adj < 0.05
  ) %>%
  filter(!is.na(region))

sirt_weighted_input <- sirt_all %>%
  mutate(
    p_safe = pmin(pmax(p_val_adj, 1e-300), 1 - 1e-10),
    z      = qnorm(p_safe / 2, lower.tail = FALSE) * sign(avg_logFC),
    weight = -log10(p_safe)
  )

## across cell types

sirt_region_input <- sirt_all %>%
  filter(cell_type %in% c("astrocyte", "macrophage", "neuron", "oligodendrocyte", "OPC")) %>%
  filter(!str_detect(experiment, "Leng: Entorhinal cortex \\(Braak 0 vs 2\\)")) %>%
  mutate(
    p_safe = pmin(pmax(p_val_adj, 1e-300), 1 - 1e-10),
    z      = qnorm(p_safe / 2, lower.tail = FALSE) * sign(avg_logFC),
    weight = -log10(p_safe)
  )

sirt_region_input$region <- replace(sirt_region_input$region, sirt_region_input$region == "PFC/BA9", "PFC")

pooled_by_region_v2 <- sirt_region_input %>%
  group_by(region) %>%
  summarise(
    pooled_logFC = weighted.mean(avg_logFC, w = weight),
    z_combined   = sum(weight * z) / sqrt(sum(weight^2)),
    pooled_p     = 2 * pnorm(abs(z_combined), lower.tail = FALSE),
    k            = n(),
    .groups = "drop"
  ) %>%
  mutate(sig_flag = pooled_p < 0.05)

val_range <- range(pooled_by_region_v2$pooled_logFC)

ggplot(pooled_by_region_v2,
                              aes(x = fct_reorder(region, pooled_logFC), y = 1, fill = pooled_logFC)) +
  geom_tile(color = "black", linewidth = 0.3, width = 1, height = 1) +
  geom_text(aes(label = ifelse(sig_flag,
                               paste0(round(pooled_logFC, 2), "*"),
                               round(pooled_logFC, 2))),
            size = 5, color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-0.5, 0.5),
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
    name = expression(paste("weighted \n mean log"[2],"FC"))
  ) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    aspect.ratio=1/6
  ) +
  labs(x = "", y = "")


ggplot(pooled_by_region_v2,
                                   aes(x = 1, y = fct_reorder(region, pooled_logFC), fill = pooled_logFC)) +
  geom_tile(color = "black", linewidth = 0.3, width = 1, height = 1) +
  geom_text(aes(label = ifelse(sig_flag,
                               paste0(round(pooled_logFC, 2), "*"),
                               round(pooled_logFC, 2))),
            size = 5.5, color = "black", angle = 90) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black", size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "")


library(ggplot2)
library(forcats)

ggplot(pooled_by_region_v2,
                                 aes(x = 1, y = fct_reorder(region, pooled_logFC), fill = pooled_logFC)) +
  geom_tile(color = "black", linewidth = 0.3, width = 1, height = 1) +
  geom_text(aes(label = ifelse(sig_flag,
                               paste0(round(pooled_logFC, 2), "*"),
                               round(pooled_logFC, 2))),
            size = 5.5, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_y_discrete(position = "right") +   # <- переносит подписи оси Y на правую сторону
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "")



