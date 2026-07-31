graphics.off()
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(scales)
  library(decoupleR)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(nanoparquet)
  library(fs)
  library(ggrepel)
  library(patchwork)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ── Configuration ────────────────────────────────────────────────────────────

MOUSE_DIR    <- "/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/mus_musculus/results/"
HUMAN_DIR    <- "/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/homo_sapiens/results/"
NETWORK_FILE <- "/tank/projects/ekashuk/sirt6_biomarker/interactions.csv"
ORTHO_FILE   <- "/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/ortholog_map_1to1.csv"

MINSIZE   <- 5
PADJ_THR  <- 0.05
COL_UP    <- "#d6604d"
COL_DN    <- "#4393c3"
COL_AGE   <- "#c0392b"     # for aging-related TFs

# Datasets to drop entirely (mouse KO technical artefacts)
DATASET_EXCLUDE <- c("GSE206513", "GSE221077")

# GSE213425 was catalogued KO but is a SIRT6 overexpression (HUVEC) experiment.
ARM_OVERRIDE <- c("GSE213425" = "OE")

SPECIES <- list(
  mouse = list(dir = MOUSE_DIR, ortho_org = "Mus musculus"),
  human = list(dir = HUMAN_DIR, ortho_org = NULL)
)

aging_related_tfs <- c(
  # nutrient sensing / FOXO / autophagy
  "FOXO1","FOXO3","FOXO4","TFEB",
  # cellular senescence / SASP
  "TP53","RB1","CEBPB","GATA4",
  # inflammaging: NF-kB / AP-1 / JAK-STAT
  "NFKB1","NFKB2","RELA","REL","RELB",
  "JUN","JUNB","JUND","FOS","FOSB","ATF3","STAT1","STAT3",
  # proteostasis: HSR / UPR / ISR
  "HSF1","XBP1","ATF4","ATF6","DDIT3","NFE2L2",
  # mitochondrial biogenesis
  "PPARGC1A","NRF1","TFAM","ESRRA",
  # epigenetic / growth
  "MYC","E2F1","EZH2",
  # circadian
  "ARNTL","CLOCK","NR1D1"
)


# Network (human CollecTRI) 

net_human <- read_tsv(NETWORK_FILE, show_col_types = FALSE) %>%
  transmute(
    source = source_genesymbol,
    target = target_genesymbol,
    mor = case_when(
      is_stimulation & !is_inhibition ~  1L,
      !is_stimulation &  is_inhibition ~ -1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(mor)) %>%
  distinct(source, target, .keep_all = TRUE)


# Gene-id -> human-symbol maps 

ortho <- read.csv(ORTHO_FILE, stringsAsFactors = FALSE)

list_de_files <- function(dir) as.character(dir_ls(dir, recurse = TRUE, glob = "*deseq2.parquet"))
strip_ver     <- function(x) sub("\\..*$", "", x)

build_ortho_map <- function(org) {
  ortho %>%
    filter(organism == org, orthology_type == "ortholog_one2one",
           !is.na(human_gene_symbol), human_gene_symbol != "") %>%
    distinct(organism_gene_id, .keep_all = TRUE) %>%
    { setNames(.$human_gene_symbol, .$organism_gene_id) }
}

build_human_map <- function(files) {
  ens <- unique(unlist(lapply(files, function(f)
    strip_ver(as.data.frame(read_parquet(f))$gene_id))))
  tab <- AnnotationDbi::select(org.Hs.eg.db, keys = ens,
                               columns = c("ENSEMBL", "SYMBOL"),
                               keytype = "ENSEMBL") %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  setNames(tab$SYMBOL, tab$ENSEMBL)
}

symbol_maps <- imap(SPECIES, function(cfg, sp) {
  if (is.null(cfg$ortho_org)) build_human_map(list_de_files(cfg$dir))
  else                        build_ortho_map(cfg$ortho_org)
})
iwalk(symbol_maps, ~ message(sprintf("  %s map: %d ids", .y, length(.x))))


# Filename metadata 

parse_meta <- function(path) {
  fn      <- basename(path)
  dataset <- basename(dirname(path))
  arm     <- if (grepl("SIRT6-OE", fn, ignore.case = TRUE)) "OE" else "KO"
  if (dataset %in% names(ARM_OVERRIDE)) arm <- ARM_OVERRIDE[[dataset]]
  stratum <- fn %>%
    str_remove(paste0("^", dataset, "_")) %>%
    str_remove("_SIRT6-.*$")
  tibble(dataset = dataset, stratum = stratum, arm = arm,
         label = if (stratum %in% c("", "all")) dataset
         else paste0(dataset, " (", stratum, ")"))
}


# Per-file ULM on `stat` 

process_file <- function(path, species) {
  meta <- parse_meta(path)
  if (meta$dataset %in% DATASET_EXCLUDE) return(NULL)
  
  de <- as.data.frame(read_parquet(path))
  if (!all(c("gene_id", "stat") %in% names(de))) {
    message("  SKIP (no stat col): ", basename(path)); return(NULL)
  }
  de$gene_id <- strip_ver(de$gene_id)
  de$hsym    <- symbol_maps[[species]][de$gene_id]
  
  de <- de %>% filter(!is.na(hsym), !is.na(stat), is.finite(stat))
  if (nrow(de) < 100) return(NULL)
  
  # collapse duplicate human symbols by the most extreme |stat|
  de <- de %>%
    group_by(hsym) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  mat <- matrix(de$stat, ncol = 1, dimnames = list(de$hsym, "stat"))
  
  acts <- tryCatch(
    run_ulm(mat, net_human,
            .source = "source", .target = "target", .mor = "mor",
            minsize = MINSIZE),
    error = function(e) { message("  ULM error (", meta$dataset, "): ",
                                  e$message); NULL }
  )
  if (is.null(acts) || nrow(acts) == 0) return(NULL)
  
  acts %>%
    transmute(source, score, p_value) %>%
    mutate(padj = p.adjust(p_value, "BH"),            # within-file, for stars
           species = species,
           dataset = meta$dataset, stratum = meta$stratum,
           arm = meta$arm, label = meta$label)
}

acts_all <- imap_dfr(SPECIES, function(cfg, sp) {
  files <- list_de_files(cfg$dir)
  message(sprintf("  %s: %d files", sp, length(files)))
  map_dfr(files, process_file, species = sp)
})

#Stouffer meta per (species, arm, TF)

stouffer_meta <- function(df) {
  df %>%
    mutate(zmag = qnorm(pmin(1 - p_value / 2, 1 - 1e-16)),
           z    = sign(score) * zmag) %>%
    group_by(species, arm, source) %>%
    summarise(
      k               = n(),
      score_mean      = mean(score),
      score_sd        = sd(score),
      se_score        = score_sd / sqrt(k),
      n_up            = sum(score > 0),
      n_dn            = sum(score < 0),
      frac_consistent = max(n_up, n_dn) / k,
      direction       = if_else(score_mean > 0, "up", "down"),
      z_meta          = sum(z) / sqrt(k),
      p_meta          = 2 * pnorm(-abs(z_meta)),
      .groups         = "drop"
    ) %>%
    group_by(species, arm) %>%
    mutate(padj_meta = p.adjust(p_meta, "BH")) %>%
    ungroup()
}

meta_tf <- stouffer_meta(acts_all)

top_tfs <- function(sp, ar, n = 20, min_k = 2) {
  meta_tf %>%
    filter(species == sp, arm == ar, k >= min_k) %>%
    arrange(padj_meta, desc(abs(z_meta))) %>%
    slice_head(n = n)
}

meta_tf %>% count(species, arm) %>% as.data.frame() %>% print()


# Concordance 

concordance <- function(df, xcol, ycol, expect = c("mirror", "same")) {
  expect <- match.arg(expect)
  x <- df[[xcol]]; y <- df[[ycol]]
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]; n <- length(x)
  if (n < 3) return(tibble(n = n, hits = NA, frac = NA,
                           binom_p = NA, rho = NA, spearman_p = NA))
  hit <- if (expect == "mirror") sign(x) != sign(y) else sign(x) == sign(y)
  bt  <- binom.test(sum(hit), n, 0.5)
  sp  <- suppressWarnings(cor.test(x, y, method = "spearman"))
  tibble(n = n, hits = sum(hit), frac = mean(hit),
         binom_p = bt$p.value, rho = unname(sp$estimate),
         spearman_p = sp$p.value)
}

# wide-per-arm (for KO<->OE mirror within a species)
mirror_wide <- function(sp) {
  meta_tf %>%
    filter(species == sp) %>%
    dplyr::select(source, arm, score_mean, se_score, frac_consistent, z_meta, padj_meta) %>%
    pivot_wider(names_from = arm,
                values_from = c(score_mean, se_score, frac_consistent,
                                z_meta, padj_meta)) %>%
    filter(!is.na(score_mean_KO), !is.na(score_mean_OE)) %>%
    mutate(
      mirror         = sign(score_mean_KO) != sign(score_mean_OE),
      combined_score = (abs(score_mean_KO) * frac_consistent_KO +
                          abs(score_mean_OE) * frac_consistent_OE) / 2
    )
}

# wide-per-species (for mouse<->human within an arm)
cross_wide <- function(ar) {
  meta_tf %>%
    filter(arm == ar) %>%
    dplyr::select(source, species, score_mean, se_score, frac_consistent, z_meta, padj_meta) %>%
    pivot_wider(names_from = species,
                values_from = c(score_mean, se_score, frac_consistent,
                                z_meta, padj_meta)) %>%
    filter(!is.na(score_mean_mouse), !is.na(score_mean_human))
}

mouse_mirror <- mirror_wide("mouse")
human_mirror <- mirror_wide("human")
cross_ko     <- cross_wide("KO")
cross_oe     <- cross_wide("OE")

theme_custom <- function(font.size = 12) {
  theme_bw(base_size = font.size) +
    theme(panel.grid.major = element_line(colour = "grey90"),
          panel.grid.minor = element_blank(),
          legend.key       = element_blank(),
          strip.background = element_rect(fill = "grey85", colour = "grey50"),
          plot.title       = element_text(hjust = 0.5))
}

star <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                              p < 0.05 ~ "*", TRUE ~ "")




# ── Scatter: KO vs OE (mirror), aging TFs highlighted red ────────────────────

make_mirror_scatter <- function(mir, sp, lim = NULL, size_limits = NULL) {
  if (nrow(mir) == 0) return(NULL)
  if (is.null(lim)) lim <- max(abs(c(mir$score_mean_KO, mir$score_mean_OE)),
                               na.rm = TRUE) * 1.1
  
  mir <- mir %>%
    mutate(
      is_aging = source %in% aging_related_tfs,
      is_rest  = source == "REST",
      color_group = case_when(
        is_aging ~ "Aging-related",
        mirror   ~ "Mirror",
        TRUE     ~ "Other"),
      lab      = if_else(is_aging | is_rest, source, ""),
      lab_col  = case_when(is_aging ~ "aging", is_rest ~ "rest", TRUE ~ NA_character_)
    )
  
  ggplot(mir, aes(score_mean_KO, score_mean_OE)) +
    xlim(-lim, lim) + ylim(-lim, lim) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = COL_UP, alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = COL_UP, alpha = 0.04) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = -1, intercept = 0, linetype = "dotted", color = "grey40") +
    geom_errorbar(aes(ymin = score_mean_OE - se_score_OE,
                      ymax = score_mean_OE + se_score_OE),
                  width = 0, color = "grey70", linewidth = 0.3, alpha = 0.6) +
    geom_errorbarh(aes(xmin = score_mean_KO - se_score_KO,
                       xmax = score_mean_KO + se_score_KO),
                   height = 0, color = "grey70", linewidth = 0.3, alpha = 0.6) +
    geom_point(aes(size = combined_score, fill = color_group),
               shape = 21, color = "white", stroke = 0.5, alpha = 0.9) +
    geom_text_repel(aes(label = lab, color = lab_col,
                        fontface = if_else(is_aging, "bold", "plain")),
                    size = 4, box.padding = 1, point.padding = 0.4,
                    force = 4, max.overlaps = 150,
                    min.segment.length = 0, segment.size = 0.3, seed = 1,
                    na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = c("Aging-related" = COL_AGE,
                                 "Mirror"        = "#4d4d4d",
                                 "Other"         = "grey75"),
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = c(aging = COL_AGE, rest = "black"),
                       na.translate = FALSE, guide = "none") +
    scale_size_continuous(range = c(2, 9), name = "Combined\nscore",
                          limits = size_limits,
                          guide = guide_legend(
                            override.aes = list(shape = 21, fill = "grey40",
                                                color = "black", stroke = 0.6))) +
    theme_custom() +
    theme(axis.text    = element_text(size = 15, colour = "black"),
          axis.title   = element_text(size = 15),
          legend.text  = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    labs(x = "Mean ULM score (SIRT6-KO vs WT)",
         y = "Mean ULM score (SIRT6-OE vs WT)")
}

mirror_size_limits <- range(c(mouse_mirror$combined_score, human_mirror$combined_score),
                            na.rm = TRUE)

p_sc_mouse <- make_mirror_scatter(mouse_mirror, "Mouse", size_limits = mirror_size_limits)
p_sc_human <- make_mirror_scatter(human_mirror, "Human", size_limits = mirror_size_limits)


p_sc_mouse
p_sc_human


# ── Scatter: cross-species (mouse vs human) per arm ──────────────────────────

make_cross_scatter <- function(cw, ar, lim = NULL, size_limits = NULL) {
  if (nrow(cw) == 0) return(NULL)
  if (is.null(lim)) lim <- max(abs(c(cw$score_mean_mouse, cw$score_mean_human)),
                               na.rm = TRUE) * 1.1
  cw <- cw %>%
    mutate(
      is_aging = source %in% aging_related_tfs,
      is_rest  = source == "REST",
      agree    = sign(score_mean_mouse) == sign(score_mean_human),
      mag      = abs(score_mean_mouse) + abs(score_mean_human),
      lab      = if_else(is_aging | is_rest, source, ""),
      lab_col  = case_when(is_aging ~ "aging", is_rest ~ "rest", TRUE ~ NA_character_)
    )
  
  ggplot(cw, aes(score_mean_mouse, score_mean_human)) +
    xlim(-lim, lim) + ylim(-lim, lim) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = COL_UP, alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = COL_DN, alpha = 0.04) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
    geom_point(aes(size = mag, fill = is_aging), shape = 21,
               color = "white", stroke = 0.5, alpha = 0.9) +
    geom_text_repel(aes(label = lab, color = lab_col,
                        fontface = if_else(is_aging, "bold", "plain")),
                    size = 5, box.padding = 0.6, max.overlaps = 100,
                    min.segment.length = 0, segment.size = 0.3, seed = 1,
                    na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = COL_AGE, `FALSE` = "grey70"),
                      labels = c(`TRUE` = "Aging-related", `FALSE` = "Other"),
                      name = NULL) +
    scale_color_manual(values = c(aging = COL_AGE, rest = "black"),
                       na.translate = FALSE, guide = "none") +
    scale_size_continuous(range = c(2, 9), name = "|mouse|+|human|",
                          limits = size_limits) +
    theme_custom() +
    theme(axis.text  = element_text(size = 16, colour = "black"),
          axis.title = element_text(size = 18)) +
    labs(x = "Mouse ULM score", y = "Human ULM score")
}

mag_ko <- abs(cross_ko$score_mean_mouse) + abs(cross_ko$score_mean_human)
mag_oe <- abs(cross_oe$score_mean_mouse) + abs(cross_oe$score_mean_human)
cross_size_limits <- range(c(mag_ko, mag_oe), na.rm = TRUE)

p_cross_ko <- make_cross_scatter(cross_ko, "KO", size_limits = cross_size_limits)
p_cross_oe <- make_cross_scatter(cross_oe, "OE", size_limits = cross_size_limits)

if (!is.null(p_cross_ko) && !is.null(p_cross_oe)) {
  p_cross_combo <- (p_cross_ko + ggtitle("KO arm") | p_cross_oe + ggtitle("OE arm")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right", plot.title = element_text(size = 18, hjust = 0.5))
  print(p_cross_combo)
} else {
  if (!is.null(p_cross_ko)) print(p_cross_ko)
  if (!is.null(p_cross_oe)) print(p_cross_oe)
}

