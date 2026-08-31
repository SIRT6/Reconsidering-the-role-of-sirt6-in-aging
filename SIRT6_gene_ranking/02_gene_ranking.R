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

# ortologs_table <- ortologs
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
extract_scores <- function(filepath, organism, mapping = TRUE) {
  
  df <- read_parquet(filepath)
  p.adj.min <- min(df$padj[df$padj != 0], na.rm = TRUE)
  
  df <- df %>% 
    mutate(padj = ifelse(padj == 0, p.adj.min, padj),
           sign_log_padj = -log10(padj) * sign(log2FoldChange)) %>%
    select(gene_id, sign_log_padj)
  
  if (organism != 'homo_sapiens' && mapping == TRUE) {
    
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



library(ape)
library(picante)

tree <- read.tree(text = "((((homo_sapiens:28.8,macaca_fascicularis:28.8):58.2,(mus_musculus:13.1,rattus_norvegicus:13.1):73.9):7,sus_scrofa:94):592,drosophila_melanogaster:686);")

# evol.distinct считает, какую долю общей эволюционной истории вносит каждый вид
pw <- evol.distinct(tree, type = "fair.proportion")

# слишком сильный перекос без корня будет
phylo_weights <- sqrt(pw$w)
phylo_weights <- setNames(phylo_weights / sum(phylo_weights), pw$Species)

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

n_boot <- 100
boot_list <- vector("list", n_boot)

for (b in seq_len(n_boot)) {
  portrait <- list()
  # i = 4
  for (i in seq_along(organisms)) {
    root_dir <- paste0("/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/", organisms[i])
    parquet_files <- list.files(path = root_dir, pattern = "KO_vs_WT_deseq2\\.parquet$", recursive = TRUE, full.names = TRUE)
    
    # cat("Number of files for", organisms[i], ':', length(parquet_files), '\n')
    score_list <- map(parquet_files, extract_scores, organism = organisms[i])
    
    combined <- reduce(score_list, full_join, by = "gene_id")
    colnames(combined)[-1] <- paste0("D_", seq_along(score_list))
    
    # print(sapply(combined, function(x) length(unique(x))))
    # remove datasets which are not suiutable for analysis
    if (organisms[i] == 'mus_musculus') {
      combined <- combined %>% select(-c('D_4','D_9','D_13','D_18'))
    }
    # cat("Number of files for", organisms[i], ':', length(combined[, -1]), '\n')
    
    idx <- sample(seq_len(ncol(combined) - 1), replace = TRUE)
    
    boot_combined <- combined[, c(1, idx + 1)]
    
    # rank сортирует по возрастанию, поэтому если нам нужно от наиболее положительному, к наименее, то нужно * на (-1)
    rank_up <- apply(boot_combined[, -1], 2, function(x) rank(-x, na.last = "keep", ties.method = "min"))
    rank_down <- apply(boot_combined[, -1], 2, function(x) rank(x, na.last = "keep", ties.method = "min"))
    
    rank_up <- as.data.frame(rank_up)
    rank_down <- as.data.frame(rank_down)
    
    up_counts <- sapply(thresholds, function(k) rowSums(rank_up <= k, na.rm = T))
    down_counts <- sapply(thresholds, function(k) rowSums(rank_down <= k, na.rm = T))
    
    theory_max <- sum(length(boot_combined[-c(1)]) * weights)
    # cat("Theory max score for", organisms[i], ':', theory_max, '\n\n')
    
    scores <- as.vector((up_counts - down_counts) %*% weights)
    final_scores <- (scores * sqrt(length(boot_combined[, -1]))) / theory_max
    
    portrait[[organisms[i]]] <- data.frame(gene_id = boot_combined$gene_id,
                                           final_score = final_scores) %>%
      rename(!!organisms[i] := final_score)
  }
  
  portrait_tmp <- reduce(portrait, full_join, by = "gene_id")
  
  score_matrix <- as.matrix(portrait_tmp[, organisms])
  portrait_tmp$final_score_p4 <- apply(score_matrix, 1, function(s) {
    present <- !is.na(s)
    if (!any(present)) return(NA)
    w <- combined_w[organisms][present]
    sum(s[present] * w)
  })
  
  boot_list[[b]] <- portrait_tmp %>%
    select(gene_id, final_score_p4)
}

boot_scores_df <- bind_rows(boot_list, .id = "rep")

ci_low  <- apply(boot_scores, 1, quantile, 0.025, na.rm = TRUE)
ci_high <- apply(boot_scores, 1, quantile, 0.975, na.rm = TRUE)

ci_df <- boot_scores_df %>%
  group_by(gene_id) %>%
  summarise(ci_low  = quantile(final_score_p4, 0.025, na.rm = TRUE),
            ci_high = quantile(final_score_p4, 0.975, na.rm = TRUE),
            .groups = "drop")

portrait <- list()
# i = 4
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
    rename(!!organisms[i] := final_score)
}

portrait_tmp <- reduce(portrait, full_join, by = "gene_id")

portrait_tmp$final_score <- rowSums(portrait_tmp[,organisms], na.rm=T)

portrait_all <- portrait_tmp %>% 
  select(c('gene_id', 'final_score_p4')) %>%
  left_join(ci_df, by = "gene_id") %>%
  mutate(direction = case_when(final_score_p4 > 0 ~ "up", final_score_p4 < 0 ~ "down", TRUE ~ "neutral")) %>%
  arrange(desc(abs(final_score_p4))) %>% 
  left_join(unique(ortologs_clean[ortologs_clean$target_organism == 'homo_sapiens', c('ortholog_ensg', 'ortholog_name')]), by = c("gene_id" = "ortholog_ensg"))

# write_tsv(portrait_all, 'output_data/signatures/portrait_all.boot.tsv')
portrait_all <- read_tsv('output_data/signatures/portrait_all.boot.tsv')


# Single organisms - advanced ranking with many ortologs  -------------------

organisms <- c('homo_sapiens', 'macaca_fascicularis', 'mus_musculus')

n_boot <- 100
portrait <- list()
# i = 2
set.seed(123)
for (i in seq_along(organisms)) {
  root_dir <- paste0("/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/", organisms[i])
  parquet_files <- list.files(path = root_dir, pattern = "KO_vs_WT_deseq2\\.parquet$", recursive = TRUE, full.names = TRUE)
  
  cat("Number of files for", organisms[i], ':', length(parquet_files), '\n')
  score_list <- map(parquet_files, extract_scores, organism = organisms[i], mapping = FALSE)
  
  combined <- reduce(score_list, full_join, by = "gene_id")
  colnames(combined)[-1] <- paste0("D_", seq_along(score_list))
  
  print(sapply(combined, function(x) length(unique(x))))
  # удаляем странные датасеты
  if (organisms[i] == 'mus_musculus') {
    combined <- combined %>% select(-c('D_4','D_9','D_13','D_18'))
  }
  cat("Number of files for", organisms[i], ':', length(combined[, -1]), '\n')
  
  boot_scores <- replicate(n_boot, {

    ## sample datasets (columns) with replacement
    idx <- sample(seq_len(ncol(combined) - 1), replace = TRUE)

    boot_combined <- combined[, c(1, idx + 1)]

    rank_up <- apply(boot_combined[, -1], 2,
                     function(x) rank(-x, na.last = "keep", ties.method = "min"))
    rank_down <- apply(boot_combined[, -1], 2,
                       function(x) rank(x, na.last = "keep", ties.method = "min"))

    up_counts <- sapply(thresholds, function(k) rowSums(rank_up <= k, na.rm = TRUE))
    down_counts <- sapply(thresholds, function(k) rowSums(rank_down <= k, na.rm = TRUE))

    theory_max <- length(idx) * sum(weights)

    scores <- as.vector((up_counts - down_counts) %*% weights)
    (scores * sqrt(length(idx))) / theory_max
  })
  
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
                                         final_score = final_scores,
                                         ci_low = apply(boot_scores, 1, quantile, 0.025, na.rm = TRUE),
                                         ci_high = apply(boot_scores, 1, quantile, 0.975, na.rm = TRUE))
}

# write_tsv(portrait[["homo_sapiens"]], 'output_data/signatures/portrait_homo_sapience.boot.tsv')
# write_tsv(portrait[["macaca_fascicularis"]], 'output_data/signatures/portrait_macaca_fascicularis.boot.tsv')
# write_tsv(portrait[["mus_musculus"]], 'output_data/signatures/portrait_mus_musculus.boot.tsv')

portrait_homo <- portrait[["homo_sapiens"]]
portrait_mus <- portrait[["mus_musculus"]]
portrait_macaca <- portrait[["macaca_fascicularis"]]

portrait_homo <- read_tsv('output_data/signatures/portrait_homo_sapience.boot.tsv')
portrait_macaca <- read_tsv('output_data/signatures/portrait_macaca_fascicularis.boot.tsv')
portrait_mus <- read_tsv('output_data/signatures/portrait_mus_musculus.boot.tsv')

create_combined_meta <- function(portrait_org, meta_file) {
  
  portrait_org$rank <- rank(-abs(portrait_org$final_score), ties.method = "min")
  
  meta_path <- paste0('/tank/projects/public_data/Sirt6_datasets/Expression/meta_results/', meta_file)
  meta_analysis <- read_csv(meta_path)
  
  combined_meta <- meta_analysis %>% 
    filter(FDR < 0.05) %>% 
    select(gene_symbol, gene_id, meta_LFC, CI_lower, CI_upper) %>% 
    left_join(portrait_org[c('final_score', 'gene_id', 'rank', 'ci_low', 'ci_high')], by = 'gene_id')
  
  combined_meta$logrank <- log10(combined_meta$rank)
  return(combined_meta)
}

combined_homo <- create_combined_meta(portrait_homo, 'meta_siggenes_KO_human.csv')
combined_mus <- create_combined_meta(portrait_mus, 'meta_siggenes_KO_mouse.csv')
combined_macaca <- create_combined_meta(portrait_macaca, 'meta_siggenes_KO_macaca.csv')

# Combined meta plot ------------------------------------------------------

library(ggplot2)
library(DOSE)
library(ggrepel)

# result <- coalesce(unlist(list1), unlist(list2))
combined_macaca$gene_symbol <- coalesce(combined_macaca$gene_symbol, combined_macaca$gene_id)
xr <- range(c(combined_macaca$ci_low, combined_macaca$ci_high), na.rm = TRUE)

scatter_plot <- ggplot(combined_macaca, aes(x = final_score, y = meta_LFC, color = pmax(logrank, 1), size = 1 / pmax(logrank, 1))) +
  geom_errorbar(data = subset(combined_macaca, rank <= 300),
                aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", linewidth = 1, alpha = 0.3) +
  geom_errorbar(data = subset(combined_macaca, rank <= 300),
                aes(ymin = CI_lower, ymax = CI_upper),
                orientation = "x", linewidth = 1, alpha = 0.3) +
  scale_x_continuous(limits = c(xr[1] - 0.15, xr[2] + 0.2)) +
  # scale_x_continuous(limits = range(c(combined_macaca$ci_low, combined_macaca$ci_high), na.rm = TRUE)) +
  scale_y_continuous(limits = range(c(combined_macaca$CI_lower, combined_macaca$CI_upper), na.rm = TRUE)) +
  geom_point() +
  geom_text(data = subset(combined_macaca, ((rank <= 300) & (final_score > 0))),
            aes(label = gene_symbol),
            nudge_y = 0.25, nudge_x = -0.25, size = 4,
            check_overlap = TRUE, color = "black") +
  geom_text(data = subset(combined_macaca, ((rank <= 300) & (final_score < 0))),
            aes(label = gene_symbol),
            nudge_y = 0.25, size = 4,
            check_overlap = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_dose(20) +
  scale_color_gradient(low = "#e49235", high = "grey70", name = "Logrank\nscore", trans = "log10") +
  scale_size_continuous(range = c(1, 5), name = "Inverse\nlogrank score") +
  labs(x = "Rank-based meta-analysis scores",
       y = "Meta-regression LFC")

# ggsave("output_plot/combined_meta/macaca_siggenes.png", scatter_plot, width = 8.5, height = 6, dpi = 300)

combined_homo$gene_symbol <- coalesce(combined_homo$gene_symbol, combined_homo$gene_id)
xr <- range(c(combined_homo$ci_low, combined_homo$ci_high), na.rm = TRUE)

scatter_plot <- ggplot(combined_homo, aes(x = final_score, y = meta_LFC, color = pmax(logrank, 1), size = 1 / pmax(logrank, 1))) +
  geom_errorbar(data = subset(combined_homo, rank <= 300),
                aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", linewidth = 1, alpha = 0.3) +
  geom_errorbar(data = subset(combined_homo, rank <= 300),
                aes(ymin = CI_lower, ymax = CI_upper),
                orientation = "x", linewidth = 1, alpha = 0.3) +
  scale_x_continuous(limits = c(xr[1] - 0.1, xr[2] + 0.05)) +
  # scale_x_continuous(limits = range(c(combined_homo$ci_low, combined_homo$ci_high), na.rm = TRUE)) +
  scale_y_continuous(limits = range(c(combined_homo$CI_lower, combined_homo$CI_upper), na.rm = TRUE)) +
  geom_point() +
  geom_text(data = subset(combined_homo, rank <= 300),
            aes(label = gene_symbol),
            nudge_y = 0.15, nudge_x = 0, size = 4,
            check_overlap = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_dose(20) +
  scale_color_gradient(low = "#e49235", high = "grey70", name = "Logrank\nscore", trans = "log10") +
  scale_size_continuous(range = c(1, 5), name = "Inverse\nlogrank score") +
  labs(x = "Rank-based meta-analysis scores",
       y = "Meta-regression LFC")

# ggsave("output_plot/combined_meta/homo_siggenes.png", scatter_plot, width = 8.5, height = 6, dpi = 300)

combined_mus$gene_symbol <- coalesce(combined_mus$gene_symbol, combined_mus$gene_id)
xr <- range(c(combined_mus$ci_low, combined_mus$ci_high), na.rm = TRUE)

scatter_plot <- ggplot(combined_mus, aes(x = final_score, y = meta_LFC, color = pmax(logrank, 1), size = 1 / pmax(logrank, 1))) +
  geom_errorbar(data = subset(combined_mus, rank <= 500),
                aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", linewidth = 1, alpha = 0.3) +
  geom_errorbar(data = subset(combined_mus, rank <= 300),
                aes(ymin = CI_lower, ymax = CI_upper),
                orientation = "x", linewidth = 1, alpha = 0.3) +
  scale_x_continuous(limits = c(xr[1], xr[2])) +
  # scale_x_continuous(limits = range(c(combined_mus$ci_low, combined_mus$ci_high), na.rm = TRUE)) +
  scale_y_continuous(limits = range(c(combined_mus$CI_lower, combined_mus$CI_upper), na.rm = TRUE)) +
  geom_point() +
  geom_text(data = subset(combined_mus, rank <= 500),
            aes(label = gene_symbol),
            nudge_y = 0.05, nudge_x = 0, size = 4,
            check_overlap = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_dose(20) +
  scale_color_gradient(low = "#e49235", high = "grey70", name = "Logrank\nscore", trans = "log10") +
  scale_size_continuous(range = c(1, 5), name = "Inverse\nlogrank score") +
  labs(x = "Rank-based meta-analysis scores",
       y = "Meta-regression LFC")

ggsave("output_plot/combined_meta/mus_siggenes.png", scatter_plot, width = 8.5, height = 6, dpi = 300)

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

# Meta-analysis plot ------------------------------------------------------

portrait_all <- read_tsv('output_data/signatures/portrait_all.boot.tsv')
portrait_all$rank <- rank(-abs(portrait_all$final_score_p4), ties.method = "min")

meta_analysis <- read_csv('input_data/meta_KO_all_genes.csv')
meta_analysis$signLFC_FDR <- sign(meta_analysis$meta_LFC)*-log10(meta_analysis$FDR)
meta_analysis$signFDR <- sign(meta_analysis$meta_LFC)*(meta_analysis$FDR)

combined_meta <- meta_analysis %>% 
  filter(FDR < 0.05) %>% 
  select(human_gene_symbol, human_gene_id, signFDR, signLFC_FDR, meta_LFC, CI_lower, CI_upper) %>% 
  left_join(portrait_all[c('final_score_p4', 'gene_id', 'rank', 'ci_low', 'ci_high')], by = c('human_gene_id' = 'gene_id'))

combined_meta$final_score_p4[combined_meta$human_gene_id == 'ENSG00000165985'] <- portrait_all$final_score_p4[portrait_all$gene_id == 'ENSMUSG00000049630']
combined_meta$rank[combined_meta$human_gene_id == 'ENSG00000165985'] <- portrait_all$rank[portrait_all$gene_id == 'ENSMUSG00000049630']

combined_meta$logrank <- log10(combined_meta$rank)

# combined_meta <- read_csv("output_data/combined_meta.csv")

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
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, alpha = 0.3) +
  geom_point(size = 2) +
  scale_color_gradient(low = "darkred", high = "grey70", name = "log10(rank)") +
  # xlim(min(portrait_all$final_score_p4), max(portrait_all$final_score_p4)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  theme_bw() +
  labs(x = 'Ranking score', y = 'Meta-analysis -log10(FDR)*sign(metaLFC)') +
  geom_text_repel(data = filter(combined_meta, abs(final_score_p4) > 0.5),
                  aes(label = human_gene_symbol))

# write_csv(combined_meta, "output_data/combined_meta.boot.csv")

