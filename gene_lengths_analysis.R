graphics.off()
rm(list = ls())

suppressPackageStartupMessages({
  library(nanoparquet); library(fs); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(ggplot2); library(metafor); library(forcats)
  library(rtracklayer); library(GenomicRanges)
})

# ── 0. Config ────────────────────────────────────────────────────────────────
BASE <- "/tank/projects/public_data/Sirt6_datasets/Expression/DE_results"
ANN  <- "/tank/projects/ekashuk/annotation"
OUT  <- "/tank/projects/ekashuk/sirt6_biomarker/results/gltd"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MIN_BASEMEAN <- 20      # length ~ expression; low counts give artefactual slope
MIN_LEN      <- 1000    # bp
MIN_GENES    <- 500     # minimum genes after filters, else skip dataset
N_BINS       <- 10
EXCLUDE      <- c("GSE206513", "GSE221077")   # technical artefact + dataset without effects

# ── 1. Contrast catalog from summary ─────────────────────────────────────────
summ <- read.csv(file.path(BASE, "summary_DE_table.csv"))
if ("X" %in% names(summ)) summ$X <- NULL
summ <- as_tibble(summ) %>%
  mutate(is_OE = str_detect(contrast, "OE"),
         is_KO = str_detect(contrast, "KO"))

message("== Contrast catalog ==")
message(sprintf("Organisms: %d | KO: %d | OE: %d",
                n_distinct(summ$organism), sum(summ$is_KO), sum(summ$is_OE)))
print(count(summ, organism, is_KO, is_OE))

# ── 2. Gene length annotations per organism ──────────────────────────────────
GTF <- list(
  "Homo sapiens"            = file.path(ANN, "Homo_sapiens.GRCh38.110.gtf.gz"),
  "Mus musculus"            = file.path(ANN, "Mus_musculus.GRCm39.110.gtf.gz"),
  "Rattus norvegicus"       = file.path(ANN, "Rattus_norvegicus.mRatBN7.2.110.gtf.gz"),
  "Sus scrofa"              = file.path(ANN, "Sus_scrofa.Sscrofa11.1.110.gtf.gz"),
  "Drosophila melanogaster" = file.path(ANN, "Drosophila_melanogaster.BDGP6.46.110.gtf.gz"),
  "Macaca fascicularis"     = file.path(ANN, "Macaca_fascicularis.Macaca_fascicularis_6.0.110.gtf.gz")
)
AUTOSOMES <- list(
  "Homo sapiens"            = as.character(1:22),
  "Mus musculus"            = as.character(1:19),
  "Rattus norvegicus"       = as.character(1:20),
  "Sus scrofa"              = as.character(1:18),
  "Drosophila melanogaster" = c("2L", "2R", "3L", "3R", "4"),
  "Macaca fascicularis"     = as.character(1:20)
)

load_lengths <- function(org) {
  gtf_path <- GTF[[org]]
  if (!file.exists(gtf_path)) {
    warning(sprintf("GTF not found for %s: %s -- skipping", org, gtf_path))
    return(NULL)
  }
  g <- import(gtf_path); g <- g[g$type == "gene"]
  gm <- tibble(gene_id  = g$gene_id,
               biotype  = g$gene_biotype,
               chr      = as.character(seqnames(g)),
               span_bp  = width(g)) %>%
    filter(biotype == "protein_coding",
           chr %in% AUTOSOMES[[org]],
           span_bp >= MIN_LEN) %>%
    mutate(log_len = log10(span_bp))
  # chromosome-style diagnostic: if empty, AUTOSOMES almost certainly wrong
  if (nrow(gm) == 0)
    warning(sprintf("%s: 0 genes after chromosome filter -- check seqnames in GTF", org))
  gm
}

message("\n== Loading gene lengths ==")
gene_meta <- map(set_names(names(GTF)), load_lengths)
gene_meta <- compact(gene_meta)   # drop organisms without GTF
orgs_ready <- names(gene_meta)
walk(orgs_ready, ~ message(sprintf("  %-24s %6d protein-coding genes",
                                   .x, nrow(gene_meta[[.x]]))))

# ── 3. Slope for a single parquet ────────────────────────────────────────────
slope_one <- function(f, org) {
  d <- tryCatch(read_parquet(f), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  need <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE")
  if (!all(need %in% names(d))) return(NULL)
  
  d <- d %>%
    mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
    inner_join(gene_meta[[org]], by = "gene_id") %>%
    filter(baseMean >= MIN_BASEMEAN, is.finite(log2FoldChange), lfcSE > 0)
  if (nrow(d) < MIN_GENES) return(NULL)
  
  w    <- 1 / d$lfcSE^2
  fit  <- lm(log2FoldChange ~ log_len, d, weights = w)
  fit2 <- lm(log2FoldChange ~ log_len + log10(baseMean), d, weights = w)  # expression control
  s    <- summary(fit)$coefficients["log_len", ]
  s2   <- summary(fit2)$coefficients["log_len", ]
  
  # fraction significantly-down among long vs short genes (top/bottom decile)
  q_hi <- quantile(d$log_len, 0.9); q_lo <- quantile(d$log_len, 0.1)
  frac <- function(mask) mean(d$padj < 0.05 & d$log2FoldChange < 0 & mask, na.rm = TRUE) /
    mean(mask)
  
  tibble(
    file        = basename(f),
    organism    = org,
    n_genes     = nrow(d),
    slope       = s[[1]],  slope_se = s[[2]],  p_slope = s[[4]],
    slope_adj   = s2[[1]], slope_adj_se = s2[[2]],
    rho         = cor(d$log_len, d$log2FoldChange, method = "spearman"),
    frac_dn_long  = frac(d$log_len >= q_hi),
    frac_dn_short = frac(d$log_len <= q_lo)
  )
}

# ── 4. Find all KO parquet files for ready organisms ─────────────────────────
find_ko <- function(org) {
  d <- file.path(BASE, tolower(gsub(" ", "_", org)), "results")
  if (!dir_exists(d)) return(character(0))
  dir_ls(d, recurse = TRUE, glob = "*deseq2.parquet") %>%
    keep(~ str_detect(.x, "SIRT6-KO"))
}

job <- tibble(organism = orgs_ready) %>%
  mutate(files = map(organism, find_ko)) %>%
  unnest(files) %>%
  filter(!str_detect(files, str_c(EXCLUDE, collapse = "|")))

message(sprintf("\n== Processing %d KO files ==", nrow(job)))

res <- map2_dfr(job$files, job$organism, slope_one) %>%
  mutate(gse = str_extract(file, "GSE[0-9]+|HRA[0-9]+"),
         organism = factor(organism, levels = names(GTF))) %>%   # or ORG_ORDER, if set
  arrange(organism, slope)

message("Slopes computed: ", nrow(res))
print(res %>% select(organism, file, n_genes, slope, slope_se, p_slope, rho))

write.csv(res, file.path(OUT, "gltd_slopes_all_organisms.csv"), row.names = FALSE)

res <- read.csv("/tank/projects/ekashuk/sirt6_biomarker/results/gltd/gltd_slopes_all_organisms.csv")
print(res %>% select(organism, file, n_genes, slope, slope_se, p_slope, rho))

library(dplyr); library(stringr); library(tibble)

res <- res %>%
  mutate(
    experiment = str_extract(file, "GSE[0-9]+(,GSE[0-9]+)?|HRA[0-9]+"),
    stratum    = file %>%
      str_remove(".*?(GSE[0-9]+(,GSE[0-9]+)?|HRA[0-9]+)_") %>%
      str_remove("_SIRT6-KO.*")
  )

DROP_EXPERIMENTS <- c("GSE206513", "GSE221077")

DROP_ARMS <- tribble(
  ~experiment,   ~drop_pattern,          ~keep_comment,
  # --- Mus musculus ---
  "GSE290902",   "Ang.?II",              "keep Saline",
  "GSE129370",   "ethanol",              "keep pair-fed",
  "GSE221092",   "^\\s*(20|60)\\s*$",    "keep t0",
  # --- Homo sapiens ---
  "GSE102813",   "RAFi",                 "keep DMSO (RAFi catches both drug arms)",
  # --- Drosophila melanogaster ---
  "GSE309387",   "100.?uM",              "keep DMSO (do NOT match on TDO2 -- would drop both)"
)

SHOULD_BE_ABSENT <- c(
  "GSE212057",  # human SIRT6-OE (A549/H1299)
  "GSE213425",  # human HUVEC, genotypes only WT/NA -- no KO arm
  "GSE191320",  # drosophila: only WT/NA, age contrast (OE in summary)
  "GSE157838",  # mouse OE
  "GSE216185",  # mouse OE
  "GSE287696"   # mouse OE
)

# ── 4. Apply filters ──────────────────────────────────────────────────────────
drop_arm_rows <- res %>%
  inner_join(DROP_ARMS, by = "experiment") %>%
  filter(str_detect(stratum, drop_pattern)) %>%
  select(file, experiment, stratum, drop_pattern)

res_keep <- res %>%
  filter(!experiment %in% DROP_EXPERIMENTS) %>%
  anti_join(drop_arm_rows %>% select(file), by = "file")


res %>% filter(experiment %in% DROP_EXPERIMENTS) %>%
  select(organism, experiment, stratum, slope) %>% print()

res %>% semi_join(drop_arm_rows, by = "file") %>%
  select(organism, experiment, stratum, slope) %>% print()

res_keep %>%
  filter(experiment %in% DROP_ARMS$experiment) %>%
  select(organism, experiment, stratum, slope) %>% print()

leaked <- res %>% filter(experiment %in% SHOULD_BE_ABSENT)

res_keep <- res_keep %>%
  mutate(flag = case_when(
    experiment == "GSE102830"                ~ "neonatal constitutive KO; tissues, RIN not checked",
    experiment == "GSE64642"                 ~ "n=2/group; Early/Late = passages, nested within dataset",
    experiment == "GSE166840"                ~ "n=2/group, low power",
    experiment == "GSE236460"                ~ "only ~10 DEGs, slope unreliable",
    str_detect(experiment, "GSE130690")      ~ "mESC -- pluripotent, unusual chromatin",
    TRUE ~ NA_character_))

res_keep %>% filter(!is.na(flag)) %>%
  distinct(organism, experiment, flag) %>% print()

cat(sprintf("\nContrasts before: %d | after: %d | removed: %d\n",
            nrow(res), nrow(res_keep), nrow(res) - nrow(res_keep)))
cat("Expected removal: 2 experiments (GSE206513, GSE221077, if present)\n",
    "  + 4 mouse arms (Ang II, ethanol, t20, t60) + 2 human RAFi + 1 dros. 100uM\n")

library(dplyr); library(ggplot2); library(forcats); library(stringr)
library(tidytext); library(metafor)

safe_rma <- function(d) {
  if (nrow(d) < 2)
    return(tibble(slope = d$slope[1], ci.lb = NA_real_, ci.ub = NA_real_,
                  pval = NA_real_, I2 = NA_real_, k = nrow(d)))
  m <- rma(yi = slope, sei = slope_se, data = d, method = "REML")
  tibble(slope = as.numeric(m$b), ci.lb = m$ci.lb, ci.ub = m$ci.ub,
         pval = m$pval, I2 = m$I2, k = m$k)
}

res_keep <- res_keep %>%
  select(-any_of("organism")) %>%
  left_join(res %>% distinct(file, organism), by = "file") %>%
  mutate(organism = as.character(organism))


if (!"experiment" %in% names(res_keep))
  res_keep <- res_keep %>%
  mutate(experiment = str_extract(file, "GSE[0-9]+(,GSE[0-9]+)?|HRA[0-9]+"))


ORG_ORDER <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
               "Sus scrofa", "Drosophila melanogaster", "Macaca fascicularis")
ORG_ORDER <- intersect(ORG_ORDER, unique(res_keep$organism))

res_keep <- res_keep %>% mutate(organism = factor(organism, levels = ORG_ORDER))


pool_keep <- res_keep %>%
  group_by(organism) %>% group_modify(~ safe_rma(.x)) %>% ungroup() %>%
  mutate(organism = as.character(organism))
overall_keep <- safe_rma(res_keep) %>% mutate(organism = "All species")
pool_keep <- bind_rows(pool_keep, overall_keep)

n_exp <- res_keep %>% group_by(organism) %>%
  summarise(n_exp = n_distinct(experiment), .groups = "drop") %>%
  mutate(organism = as.character(organism))

clean_lab <- function(x) x %>%
  str_remove("_SIRT6-KO_vs_WT_deseq2\\.parquet") %>%
  str_replace_all("_", " ") %>% str_squish()

pts <- res_keep %>%
  mutate(ci.lb = slope - 1.96*slope_se,
         ci.ub = slope + 1.96*slope_se,
         label = clean_lab(file),
         organism = as.character(organism),
         type = "dataset")

dia <- pool_keep %>%
  inner_join(n_exp, by = "organism") %>%
  filter(n_exp >= 2, !is.na(ci.lb)) %>%
  transmute(organism, slope, ci.lb, ci.ub, label = "Pooled (RE)", type = "pool")

plot_df <- bind_rows(pts, dia) %>%
  mutate(organism = factor(organism, levels = ORG_ORDER)) %>%
  group_by(organism) %>%
  mutate(ord = ifelse(type == "pool", min(slope, na.rm = TRUE) - 1, slope),
         y   = reorder_within(label, ord, organism)) %>%
  ungroup()

xr <- range(c(pts$ci.lb, pts$ci.ub), na.rm = TRUE)


plot_forest_gltd <- function(es, class_order, title = NULL,
                             x_lab = "Slope: log2FC per decade of gene length",
                             xlim = NULL, alim = NULL, text_cex = 1,
                             min_k_for_pool = 2, overall_pool = TRUE,
                             open_device = TRUE, device_width = 11) {
  
  stopifnot(all(c("yi", "vi", "family", "class") %in% names(es)))
  
  class_colors <- c(
    "Homo sapiens"             = "pink2",
    "Mus musculus"             = "skyblue2",
    "Rattus norvegicus"        = "lightgreen",
    "Sus scrofa"               = "mediumpurple1",
    "Macaca fascicularis"      = "coral2",
    "Drosophila melanogaster"  = "#A65628"
  )
  
  es <- es %>%
    mutate(class = factor(class, levels = class_order)) %>%
    filter(!is.na(class)) %>%
    arrange(class, yi)
  
  family_width <- max(nchar(es$family), na.rm = TRUE)
  es$label <- sprintf(paste0("%-", family_width, "s"), es$family)
  header   <- sprintf(paste0("%-", family_width, "s"), "Dataset")
  
  present_classes <- intersect(class_order, unique(as.character(es$class)))
  
  block_idx <- split(seq_len(nrow(es)), es$class)[present_classes]
  cursor <- 0
  block_rows <- list()
  for (cl in present_classes) {
    k <- length(block_idx[[cl]])
    block_rows[[cl]] <- list(idx = block_idx[[cl]], k = k, pool = k >= min_k_for_pool)
    cursor <- cursor + k + block_rows[[cl]]$pool + 1
  }
  total_rows <- cursor - 1
  
  cursor <- total_rows
  rows_study <- setNames(rep(NA_real_, nrow(es)), NULL)
  pool_fits <- list(); pool_rows <- c(); block_ymin <- c(); block_ymax <- c()
  
  for (cl in present_classes) {
    b <- block_rows[[cl]]
    ymax_cl <- cursor + 0.5
    for (j in seq_len(b$k)) { rows_study[b$idx[j]] <- cursor; cursor <- cursor - 1 }
    if (b$pool) {
      fit <- rma(yi = yi, vi = vi, data = es[b$idx, ], method = "REML")
      pool_fits[[cl]] <- fit; pool_rows[cl] <- cursor; cursor <- cursor - 1
    }
    block_ymin[cl] <- cursor + 0.5
    block_ymax[cl] <- ymax_cl
    cursor <- cursor - 1
  }
  es$row <- rows_study
  
  overall_fit <- if (overall_pool) rma(yi = yi, vi = vi, data = es, method = "REML") else NULL
  overall_row <- if (overall_pool) min(c(unlist(pool_rows), rows_study)) - 3 else NA_real_
  stats_row1  <- overall_row - 1.3
  stats_row2  <- overall_row - 2.2
  
  if (is.null(alim)) {
    max_abs <- max(abs(c(es$yi - 1.96*sqrt(es$vi), es$yi + 1.96*sqrt(es$vi))), na.rm = TRUE)
    alim <- c(-max_abs*1.1, max_abs*1.1)
  }
  if (is.null(xlim)) xlim <- c(alim[1]*3.2, alim[2]*1.6)
  x_breaks <- pretty(alim)
  x_breaks <- x_breaks[abs(x_breaks) <= max(abs(alim))*1.001]
  
  ylim_top    <- total_rows + 3.5
  ylim_bottom <- if (overall_pool) stats_row2 - 1.2 else min(pool_rows, rows_study, na.rm = TRUE) - 1.5
  ylim <- c(ylim_bottom, ylim_top)
  
  if (open_device) {
    height_in <- max(6, total_rows * 0.28 + 3)
    if (Sys.getenv("RSTUDIO") == "1") {
      dev.new(width = device_width, height = height_in)
    } else if (capabilities("X11") || capabilities("aqua")) {
      dev.new(width = device_width, height = height_in)
    } else {
      message("No interactive graphics device -- saving gltd_forest.pdf")
      pdf("gltd_forest.pdf", width = device_width, height = height_in)
      on.exit(dev.off(), add = TRUE)
    }
  }
  
  # more margin on both sides: left for "All species (RE)"/stats,
  # right for "Drosophila melanogaster"
  par(mar = c(4, 2, 4, 17), family = "mono", font.main = 1)
  
  forest(x = es$yi, vi = es$vi, slab = es$label,
         rows = es$row, ylim = ylim,
         xlab = x_lab, main = title,
         cex.main = 1.4*text_cex, refline = 0,
         cex = 0.85*text_cex, cex.lab = 1*text_cex, cex.axis = 1*text_cex,
         header = header, alim = alim, xlim = xlim, at = x_breaks,
         annotate = FALSE)
  
  usr <- par("usr")
  
  class_bg <- setNames(adjustcolor(unname(class_colors[present_classes]), alpha.f = 0.25),
                       present_classes)
  for (cl in present_classes) {
    rect(xleft = usr[1], ybottom = block_ymin[cl], xright = usr[2], ytop = block_ymax[cl],
         col = class_bg[cl], border = NA)
  }
  with(es, points(yi, row, pch = 15, cex = 0.7*text_cex))
  with(es, segments(yi - 1.96*sqrt(vi), row, yi + 1.96*sqrt(vi), row))
  
  half_height <- 0.38
  for (cl in present_classes) {
    if (!block_rows[[cl]]$pool) next
    fit <- pool_fits[[cl]]; row <- pool_rows[cl]
    polygon(x = c(fit$ci.lb, fit$b[1], fit$ci.ub, fit$b[1]),
            y = row + c(0, half_height, 0, -half_height),
            col = unname(class_colors[cl]), border = "black", lwd = 1.1)
    text(x = usr[1], y = row, "Pooled (RE)", pos = 4, xpd = NA,
         cex = 0.75*text_cex, font = 3)
  }
  
  for (cl in present_classes) {
    mid_y <- (block_ymin[cl] + block_ymax[cl]) / 2
    text(x = usr[2] + 0.06*diff(usr[1:2]), y = mid_y, labels = cl,
         xpd = NA, adj = 0, font = 4, cex = 0.95*text_cex,
         col = adjustcolor(unname(class_colors[cl]), red.f = .6, green.f = .6, blue.f = .6))
  }
  
  if (overall_pool) {
    hh <- 0.5
    polygon(x = c(overall_fit$ci.lb, overall_fit$b[1], overall_fit$ci.ub, overall_fit$b[1]),
            y = overall_row + c(0, hh, 0, -hh),
            col = "firebrick", border = "black", lwd = 1.3)
    text(x = usr[1], y = overall_row, "All species (RE)", pos = 4, xpd = NA,
         cex = 0.9*text_cex, font = 2)
    
    text(x = usr[1], y = stats_row1, pos = 4, xpd = NA, cex = 0.8*text_cex,
         "Overall random-effects model:")
    text(x = usr[1], y = stats_row2, pos = 4, xpd = NA, cex = 0.8*text_cex,
         bquote(I^2 ~ "=" ~ .(formatC(overall_fit$I2, digits = 1, format = "f")) * "%," ~
                  tau^2 ~ "=" ~ .(formatC(overall_fit$tau2, digits = 4, format = "f")) * "," ~
                  "k =" ~ .(overall_fit$k)))
  }
  
  invisible(list(overall = overall_fit, by_class = pool_fits))
}

clean_lab <- function(x) x %>%
  str_remove("_SIRT6-KO_vs_WT_deseq2\\.parquet") %>%
  str_replace_all("_", " ") %>% str_squish()

es <- res_keep %>%
  transmute(
    family = clean_lab(file),
    class  = as.character(organism),
    yi     = slope,
    vi     = slope_se^2
  )

ORG_ORDER <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
               "Sus scrofa", "Drosophila melanogaster", "Macaca fascicularis")
ORG_ORDER <- intersect(ORG_ORDER, unique(es$class))

pdf(file = "/tank/projects/ekashuk/sirt6_biomarker/results/gltd/forest_gltd_all_species.pdf",
    width = 9, height = 12)

plot_forest_gltd(es, ORG_ORDER,
                 open_device = FALSE)

dev.off()