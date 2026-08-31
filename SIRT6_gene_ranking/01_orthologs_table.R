library(data.table)
library(arrow)
library(gprofiler2)

# 1. Define directories


root_dir <- "/tank/projects/public_data/Sirt6_datasets/Expression/DE_results"
out_dir <- "D:/skoltech/SIRT6_aging/data/gene_lists"


# 2. Find and read all DE tables


de_files <- list.files(root_dir, pattern = "deseq2\\.parquet$", recursive = TRUE, full.names = TRUE)
cat("Number of DE files found:", length(de_files), "\n")

all_de_list <- lapply(de_files, read_parquet)
all_de <- rbindlist(all_de_list, fill = TRUE)

cat("Total rows in merged DE table:", nrow(all_de), "\n")


# 3. Define organisms and g:Profiler codes


organisms <- c("drosophila_melanogaster", "homo_sapiens", "macaca_fascicularis", "mus_musculus",
               "rattus_norvegicus", "sus_scrofa")

gprofiler_codes <- c(drosophila_melanogaster = "dmelanogaster",
                     homo_sapiens = "hsapiens",
                     macaca_fascicularis = "mfascicularis",
                     mus_musculus = "mmusculus",
                     rattus_norvegicus = "rnorvegicus",
                     sus_scrofa = "sscrofa")

gene_id_col <- "gene_id"


# 4. Build gene lists per organism


gene_lists <- lapply(organisms, function(org) {
  ids <- all_de[organism == org, get(gene_id_col)]
  unique(ids)
})

names(gene_lists) <- organisms

gene_counts <- data.table(organism = names(gene_lists),
                          n_unique_genes = vapply(gene_lists, length, integer(1)))

print(gene_counts)


# 5. Function: one source organism -> one target organism


get_orthologs_pair_gprofiler <- function(gene_ids,
                                         source_organism,
                                         target_organism,
                                         chunk_size = 1000) {
  
  cat("\n=============================\n")
  cat("Processing:", source_organism, "->", target_organism, "\n")
  cat("=============================\n")
  
  gene_ids <- unique(gene_ids)
  
  cat("Total unique input genes:", length(gene_ids), "\n")
  
  source_code <- gprofiler_codes[[source_organism]]
  target_code <- gprofiler_codes[[target_organism]]
  
  chunks <- split(gene_ids, ceiling(seq_along(gene_ids) / chunk_size))
  res_list <- vector("list", length(chunks))
  
  for (i in seq_along(chunks)) {
    cat("Chunk", i, "of", length(chunks), "- genes:", length(chunks[[i]]), "\n")
    
    res <- gorth(query = chunks[[i]],
                 source_organism = source_code,
                 target_organism = target_code,
                 mthreshold = Inf,
                 filter_na = FALSE)
    
    res_list[[i]] <- as.data.table(res)
  }
  
  orth <- rbindlist(res_list, fill = TRUE)
  
  orth$source_organism <- source_organism
  orth$target_organism <- target_organism
  
  setcolorder(orth, c("source_organism", "target_organism",
                      setdiff(names(orth), c("source_organism", "target_organism"))))
  
  cat("Rows returned:", nrow(orth), "\n")
  cat("Input genes returned:", uniqueN(orth$input), "\n")
  
  return(orth)
}


# 6. Build all pairwise source -> target combinations


pair_grid <- CJ(source_organism = organisms,
                target_organism = organisms)

pair_grid <- pair_grid[source_organism != target_organism]

cat("Number of pairwise organism comparisons:", nrow(pair_grid), "\n")


# 7. Run g:Profiler orthology mapping for all pairs


all_pairwise_list <- vector("list", nrow(pair_grid))

for (i in seq_len(nrow(pair_grid))) {
  
  source_org <- pair_grid$source_organism[i]
  target_org <- pair_grid$target_organism[i]
  
  orth_pair <- get_orthologs_pair_gprofiler(gene_ids = gene_lists[[source_org]],
                                            source_organism = source_org,
                                            target_organism = target_org,
                                            chunk_size = 1000)
  
  all_pairwise_list[[i]] <- orth_pair
  
  pair_file <- file.path(out_dir, paste0("orthologs_", source_org, "_to_", target_org, ".tsv"))
}

all_pairwise_orthologs <- rbindlist(all_pairwise_list, fill = TRUE)

combined_file <- file.path(out_dir, "all_pairwise_orthologs_gprofiler2.tsv")
fwrite(all_pairwise_orthologs, file = combined_file, sep = "\t")
