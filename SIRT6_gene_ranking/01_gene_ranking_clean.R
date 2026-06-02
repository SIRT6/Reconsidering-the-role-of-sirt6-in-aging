library(arrow)
library(tidyverse)
library(data.table)

# Ortologs -------------------------------------------------------

ortologs <- read_tsv("output_data/ortologs_tables/all_pairwise_orthologs_gprofiler2.tsv") #721835

root_dir <- "/tank/projects/public_data/Sirt6_datasets/Expression/DE_results"
parquet_files <- list.files(path = root_dir, pattern = "KO_vs_WT_deseq2\\.parquet$", recursive = TRUE, full.names = TRUE)
cat("Number of DE files found:", length(parquet_files), "\n")

read_de_table <- function(file_path) {
  read_parquet(file_path)
}

all_de_list <- lapply(parquet_files, read_de_table)
all_de <- rbindlist(all_de_list, fill = TRUE)

length(unique(all_de$gene_id)) # 114957

# writeLines(unique(all_de[all_de$organism == 'drosophila_melanogaster', gene_id]), './output_data/ortologs_tables/gene_lists/drosophila_melanogaster.txt')
# writeLines(unique(all_de[all_de$organism == 'homo_sapiens', gene_id]), './output_data/ortologs_tables/gene_lists/homo_sapiens.txt')
# writeLines(unique(all_de[all_de$organism == 'macaca_fascicularis', gene_id]), './output_data/ortologs_tables/gene_lists/macaca_fascicularis.txt')
# writeLines(unique(all_de[all_de$organism == 'mus_musculus', gene_id]), './output_data/ortologs_tables/gene_lists/mus_musculus.txt')
# writeLines(unique(all_de[all_de$organism == 'rattus_norvegicus', gene_id]), './output_data/ortologs_tables/gene_lists/rattus_norvegicus.txt')
# writeLines(unique(all_de[all_de$organism == 'sus_scrofa', gene_id]), './output_data/ortologs_tables/gene_lists/sus_scrofa.txt')

clean_ortologs <- function(ortologs_table) {
  ortologs_table <- ortologs_table[ortologs_table$ortholog_ensg != 'N/A',]
  
  ortologs_sm <- ortologs_table[(ortologs_table$ortholog_ensg %in% unique(all_de$gene_id)) &
                            (ortologs_table$input_ensg %in% unique(all_de$gene_id)) ,]
  
  orth_pairs <- ortologs_sm %>%
    group_by(source_organism, target_organism, input_ensg) %>%
    mutate(n_targets = n_distinct(ortholog_ensg)) %>%
    ungroup()
  
  orth_pairs <- orth_pairs %>%
    group_by(source_organism, target_organism, ortholog_ensg) %>%
    mutate(n_sources = n_distinct(input_ensg)) %>%
    ungroup()
  
  orth_unique <- orth_pairs %>%
    filter(n_targets == 1, n_sources == 1)
  
  return(orth_unique)
}

ortologs_clean <- clean_ortologs(ortologs) # 232200

priority_order <- c("homo_sapiens", "mus_musculus", "macaca_fascicularis", "rattus_norvegicus",
                    "sus_scrofa", "drosophila_melanogaster")

ortolog_tables <- list()

for (target_org in priority_order) {
  
  allowed_sources <- priority_order[seq_along(priority_order) < match(target_org, priority_order)]
  
  dt <- ortologs_clean[ortologs_clean$target_organism == target_org &
                         ortologs_clean$source_organism %in% allowed_sources,]
  
  dt <- dt %>%
    pivot_wider(id_cols = c(ortholog_ensg, ortholog_name),
                names_from = source_organism,
                values_from = input_ensg)
  
  ortolog_tables[[target_org]] <- dt
}

# All organisms advanced ranking with many ortologs -----------------------

thresholds <- seq(1000, 8000, by = 1000)
weights <- 10^(-seq_along(thresholds) + 1)

organisms <- c('drosophila_melanogaster', 'homo_sapiens', 'macaca_fascicularis', 'mus_musculus', 'rattus_norvegicus', 'sus_scrofa')

# filepath = parquet_files[[38]]
# organism = organisms[[5]]
extract_scores <- function(filepath, organism) {
  
  df <- read_parquet(filepath)
  p.adj.min <- min(df$padj[df$padj != 0], na.rm = TRUE)
  
  df <- df %>% 
    mutate(padj = ifelse(padj == 0, p.adj.min, padj),
           sign_log_padj = -log10(padj) * sign(log2FoldChange)) %>%
    select(gene_id, sign_log_padj)
  
  if (organism != 'homo_sapiens') {
    
    mapping_table <- ortolog_tables[[organism]]
    
    df <- df %>%
      left_join(mapping_table, by = c("gene_id" = "ortholog_ensg"))
    
    source_cols <- intersect(priority_order, colnames(mapping_table))
    df$gene_id <- do.call(coalesce, c(df[source_cols], list(df$gene_id)))
    # возможно имеет смысл дропнуть гены, которые не имеют mapping с другими организмами (но с текущей версией таблички ортологов это сложно)
    df <- df %>%
      select(gene_id, sign_log_padj)
  }
  return(df)
}

portrait <- list()
# i = 2
for (i in seq_along(organisms)) {
  root_dir <- paste0("/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/", organisms[i])
  parquet_files <- list.files(path = root_dir, pattern = "KO_vs_WT_deseq2\\.parquet$", recursive = TRUE, full.names = TRUE)
  
  cat("Number of files for", organisms[i], ':', length(parquet_files), '\n')
  score_list <- map(parquet_files, extract_scores, organism = organisms[i])
  
  combined <- reduce(score_list, full_join, by = "gene_id")
  colnames(combined)[-1] <- paste0("D_", seq_along(score_list))
  
  print(sapply(combined, function(x) length(unique(x))))
  # удаляем странные датасеты
  if (organisms[i] == 'mus_musculus') {
    combined <- combined %>% select(-c('D_4','D_9','D_13','D_18'))
  }
  cat("Number of files for", organisms[i], ':', length(combined[, -1]), '\n')
  
  # rank сортирует по возрастанию, поэтому если нам нужно от наиболее положительному, к наименее, то нужно * на (-1)
  rank_up <- apply(combined[, -1], 2, function(x) rank(-x, na.last = "keep", ties.method = "min"))
  rank_down <- apply(combined[, -1], 2, function(x) rank(x, na.last = "keep", ties.method = "min"))
  
  rank_up <- as.data.frame(rank_up)
  rank_down <- as.data.frame(rank_down)
  
  up_counts <- sapply(thresholds, function(k) rowSums(rank_up <= k, na.rm = T))
  down_counts <- sapply(thresholds, function(k) rowSums(rank_down <= k, na.rm = T))
  
  theory_max <- sum(length(combined[-c(1)]) * weights)
  cat("Theory max score for", organisms[i], ':', theory_max, '\n\n')
  
  scores <- as.vector((up_counts - down_counts) %*% weights)
  final_scores <- (scores * sqrt(length(combined[, -1]))) / theory_max
  portrait[[organisms[i]]] <- data.frame(gene_id = combined$gene_id,
                                         final_score = final_scores) %>% 
    arrange(desc(abs(final_score))) %>%
    rename(!!organisms[i] := final_score)
}

portrait_tmp <- reduce(portrait, full_join, by = "gene_id")

portrait_tmp$final_score <- rowSums(portrait_tmp[,organisms], na.rm=T)

portrait_all <- portrait_tmp %>% 
  select(c('gene_id', 'final_score')) %>%
  mutate(direction = case_when(final_score > 0 ~ "up", final_score < 0 ~ "down", TRUE ~ "neutral")) %>%
  arrange(desc(abs(final_score))) %>% 
  left_join(unique(ortologs_clean[ortologs_clean$target_organism == 'homo_sapiens', c('ortholog_ensg', 'ortholog_name')]), by = c("gene_id" = "ortholog_ensg"))


# Include phylogeny -------------------------------------------------------

# phylogeny 1. Получается ерунда, так как если только у одного вида есть ген, то тогда итоговый скор равен скору этого организма

library(ape)
library(picante)

tree <- read.tree(text = "((((homo_sapiens:28.8,macaca_fascicularis:28.8):58.2,(mus_musculus:13.1,rattus_norvegicus:13.1):73.9):7,sus_scrofa:94):592,drosophila_melanogaster:686);")

# evol.distinct считает, какую долю общей эволюционной истории вносит каждый вид
pw <- evol.distinct(tree, type = "fair.proportion")

# слишком сильный перекос без корня будет
phylo_weights <- sqrt(pw$w)
phylo_weights <- setNames(phylo_weights / sum(phylo_weights), pw$Species)

score_matrix <- as.matrix(portrait_tmp[, organisms])

portrait_tmp$final_score_p1 <- apply(score_matrix, 1, function(s) {
  present <- !is.na(s)
  if (!any(present)) return(NA)
  w <- phylo_weights[organisms][present]
  w <- w / sum(w)          # renormalize to sum=1 for available species only
  sum(s[present] * w)
})

# phylogeny 2 (не перенормализовываем веса, получается, что гены, которые есть у всех организмов и скоры в одном направлении - будут автоматически значимее)

portrait_tmp$final_score_p2 <- apply(score_matrix, 1, function(s) {
  present <- !is.na(s)
  if (!any(present)) return(NA)
  w <- phylo_weights[organisms][present]
  sum(s[present] * w)
})

# phylogeny 3

n_datasets <- c(homo_sapiens = 7,
                macaca_fascicularis = 7,
                mus_musculus = 17,
                rattus_norvegicus = 1,
                sus_scrofa = 1,
                drosophila_melanogaster = 2)

reliability_w <- sqrt(n_datasets)
reliability_w <- reliability_w / sum(reliability_w)

combined_w <- phylo_weights[names(reliability_w)] * reliability_w
combined_w <- combined_w / sum(combined_w)

data.frame(organism = names(combined_w),
           n_datasets = n_datasets[names(combined_w)],
           phylo = round(phylo_weights[names(combined_w)], 3),
           reliability = round(reliability_w, 3),
           combined = round(combined_w, 3)) %>% 
  arrange(desc(combined))

portrait_tmp$final_score_p3 <- apply(score_matrix, 1, function(s) {
  present <- !is.na(s)
  if (!any(present)) return(NA)
  w <- combined_w[organisms][present]
  w <- w / sum(w)
  sum(s[present] * w)
})

# phylogenia 4 (не перенормализовываем веса, получается, что гены, которые есть у всех организмов и скоры в одном направлении - будут автоматически значимее)

portrait_tmp$final_score_p4 <- apply(score_matrix, 1, function(s) {
  present <- !is.na(s)
  if (!any(present)) return(NA)
  w <- combined_w[organisms][present]
  sum(s[present] * w)
})

# phylogeny 5

# Variance-covariance matrix: diagonal = distance from root (var),
# off-diagonal = shared branch length (cov, i.e. how related species are)
C <- vcv(tree)
C_inv <- solve(C)

# reorder matrix to match your organisms vector
org_order <- c("homo_sapiens", "macaca_fascicularis", "mus_musculus",
               "rattus_norvegicus", "sus_scrofa", "drosophila_melanogaster")
C_inv <- C_inv[org_order, org_order]

score_matrix_tmp <- as.matrix(portrait_tmp[, organisms])
score_matrix_tmp[is.na(score_matrix_tmp)] <- 0  # NAs → 0

ones <- matrix(1, nrow = length(organisms))
normalizer <- as.numeric(t(ones) %*% C_inv %*% ones)

portrait_tmp$final_score_pgls <- apply(score_matrix_tmp, 1, function(s) {
  s <- matrix(s, ncol = 1)
  as.numeric(t(ones) %*% C_inv %*% s) / normalizer
})

# write_tsv(portrait_all, 'output_data/signatures/portrait_all.tsv')

# Enrichment --------------------------------------------------------------

aging <- read_tsv("input_data/gene-aging-mechanisms.tsv", skip = 1, col_names = F)
colnames(aging) <- c('symbol', 'hallmark')

aging_long <- aging %>%
  filter(hallmark != "") %>%
  mutate(hallmark = gsub("'", "", hallmark)) %>%
  separate_rows(hallmark, sep = ",") %>%
  mutate(hallmark = str_trim(hallmark))

aging_pathways <- split(aging_long$symbol, aging_long$hallmark)

scores <- portrait_all %>%
  filter(grepl('ENSG', gene_id), !is.na(ortholog_name), ortholog_name != 'N/A') %>% 
  select(ortholog_name, final_score_p4)

geneList <- scores$final_score
names(geneList) <- scores$ortholog_name

geneList <- sort(geneList, decreasing = TRUE)

library(fgsea)

fgseaRes <- fgsea(pathways = aging_pathways,
                  stats = geneList,
                  minSize = 5,
                  maxSize = 500,
                  nperm = 10000)

plot_df <- fgseaRes %>%
  arrange(NES)

ggplot(plot_df, aes(x = NES, y = reorder(pathway, NES), size = size, color = padj)) +
  labs(y = 'Aging hallmarks') +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 14))

# Organism specific analysis + Anya comparison ----------------------------------------------

mus_siggenes <- read_csv('/tank/projects/public_data/Sirt6_datasets/Expression/meta_results/meta_siggenes_KO_mouse.csv')
human_siggenes <- read_csv('/tank/projects/public_data/Sirt6_datasets/Expression/meta_results/meta_siggenes_KO_human.csv')
macaca_siggenes <- read_csv('/tank/projects/public_data/Sirt6_datasets/Expression/meta_results/meta_siggenes_KO_macaca.csv')

cat('Number of ortologs mapped to human out of 100:', 
    sum(grepl('ENSG', portrait[['mus_musculus']]$gene_id[1:100])))
cat('Number of ortologs mapped to human out of 1000:', 
    sum(grepl('ENSG', portrait[['mus_musculus']]$gene_id[1:1000])))
cat('Number of significant intersected genes out of 100:',
    length(intersect(mus_siggenes$human_gene_id, portrait[['mus_musculus']]$gene_id[1:100])))
cat('Number of significant intersected genes out of 1000:',
    length(intersect(mus_siggenes$human_gene_id, portrait[['mus_musculus']]$gene_id[1:1000])))

# writeLines(portrait[['mus_musculus']]$gene_id[1:100], 'output_data/signatures/mus_musculus_top100.txt')
# writeLines(portrait[['mus_musculus']]$gene_id[1:1000], 'output_data/signatures/mus_musculus_top1000.txt')
# writeLines((portrait[['homo_sapiens']] %>% filter(homo_sapiens > 0))$gene_id[1:100], 'output_data/signatures/homo_sapiens_top100_up.txt')
# writeLines((portrait[['homo_sapiens']] %>% filter(homo_sapiens > 0))$gene_id[1:1000], 'output_data/signatures/homo_sapiens_top1000_up.txt')
# writeLines((portrait[['homo_sapiens']] %>% filter(homo_sapiens < 0))$gene_id[1:100], 'output_data/signatures/homo_sapiens_top100_down.txt')
# writeLines((portrait[['homo_sapiens']] %>% filter(homo_sapiens < 0))$gene_id[1:1000], 'output_data/signatures/homo_sapiens_top1000_down.txt')

cat('Number of significant intersected genes out of 100:',
    length(intersect(human_siggenes$human_gene_id, portrait[['homo_sapiens']]$gene_id[1:100])))
cat('Number of significant intersected genes out of 1000:',
    length(intersect(human_siggenes$human_gene_id, portrait[['homo_sapiens']]$gene_id[1:1000])))

cat('Number of ortologs mapped to human out of 100:', 
    sum(grepl('ENSG', portrait[['macaca_fascicularis']]$gene_id[1:100])))
cat('Number of ortologs mapped to human out of 1000:', 
    sum(grepl('ENSG', portrait[['macaca_fascicularis']]$gene_id[1:1000])))
cat('Number of significant intersected genes out of 100:',
    length(intersect(macaca_siggenes$human_gene_id, portrait[['macaca_fascicularis']]$gene_id[1:100])))
cat('Number of significant intersected genes out of 1000:',
    length(intersect(macaca_siggenes$human_gene_id, portrait[['macaca_fascicularis']]$gene_id[1:1000])))

# Compare top genes with Anya ---------------------------------------------

portrait_all <- portrait_tmp %>% 
  select(c('gene_id', 'final_score', 'final_score_p1' , 'final_score_p2',
           'final_score_p3' , 'final_score_p4', 'final_score_pgls')) %>%
  mutate(direction = case_when(final_score > 0 ~ "up", final_score < 0 ~ "down", TRUE ~ "neutral")) %>%
  arrange(desc(abs(final_score_p4))) %>% 
  left_join(unique(ortologs_clean[ortologs_clean$target_organism == 'homo_sapiens', c('ortholog_ensg', 'ortholog_name')]), by = c("gene_id" = "ortholog_ensg"))

anya <- readLines('output_data/signatures_meta_KO_Anna.txt')

top30_or <- portrait_all[1:30, 'ortholog_name']
top100_or <- portrait_all[1:100, 'ortholog_name']
top1000_or <- portrait_all[1:1000, 'ortholog_name']

intersect(anya, top30_or) # 2
intersect(anya, top100_or) # 3
intersect(anya, top1000_or) # 10

# Meta-analysis plot ------------------------------------------------------

portrait_all$rank <- rank(-abs(portrait_all$final_score_p4), ties.method = "min")

meta_analysis <- read_csv('input_data/meta_KO_all_genes.csv')
meta_analysis$signLFC_FDR <- sign(meta_analysis$meta_LFC)*-log10(meta_analysis$FDR)
meta_analysis$signFDR <- sign(meta_analysis$meta_LFC)*(meta_analysis$FDR)

combined_meta <- meta_analysis %>% 
  filter(FDR < 0.05) %>% 
  select(human_gene_symbol, human_gene_id, signFDR, signLFC_FDR) %>% 
  left_join(portrait_all[c('final_score_p4', 'gene_id', 'rank')], by = c('human_gene_id' = 'gene_id'))

combined_meta$final_score_p4[combined_meta$human_gene_id == 'ENSG00000165985'] <- portrait_all$final_score_p4[portrait_all$gene_id == 'ENSMUSG00000049630']
combined_meta$rank[combined_meta$human_gene_id == 'ENSG00000165985'] <- portrait_all$rank[portrait_all$gene_id == 'ENSMUSG00000049630']

combined_meta$logrank <- log10(combined_meta$rank)

ggplot(combined_meta, aes(x = final_score_p4, y = signFDR, color = logrank)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "darkred", high = "grey70", name = "log10(rank)") +
  xlim(min(portrait_all$final_score_p4), max(portrait_all$final_score_p4)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  theme_bw() +
  labs(x = 'Ranking score', y = 'Meta-analysis signed FDR') +
  geom_text(data = filter(combined_meta, abs(final_score_p4) > 0.5),
            aes(label = human_gene_symbol),
            nudge_y = 0.005, nudge_x = -0.02)

library(ggrepel)
ggplot(combined_meta, aes(x = final_score_p4, y = signLFC_FDR, color = logrank)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "darkred", high = "grey70", name = "log10(rank)") +
  xlim(min(portrait_all$final_score_p4), max(portrait_all$final_score_p4)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  theme_bw() +
  labs(x = 'Ranking score', y = 'Meta-analysis -log10(FDR)*sign(metaLFC)') +
  geom_text_repel(data = filter(combined_meta, abs(final_score_p4) > 0.5),
                  aes(label = human_gene_symbol))
