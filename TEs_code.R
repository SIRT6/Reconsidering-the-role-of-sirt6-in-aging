graphics.off()
rm(list = ls())

suppressPackageStartupMessages({
  library(nanoparquet); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(forcats); library(readr)
  library(DESeq2); library(metafor)
  library(ggplot2); library(ggrepel); library(pheatmap)
})


theme_custom <- function(font.size = 12) {
  ggplot2::theme_bw(base_size = font.size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey85", colour = "grey50"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

# Metadata
tele_dir <- "/tank/projects/ekashuk/TE/telescope_output_mouse/"
samples <- read_parquet("/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/sus_scrofa/samples.parquet")
s2e     <- read_parquet("/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/sus_scrofa/samples_to_experiment.parquet")

meta <- samples %>%
  left_join(s2e, by = "sample_id") %>%
  mutate(
    genotype  = factor(genotype,  levels = c("WT", "SIRT6-KO", "SIRT6-OE/K3Q")),
    cell_type = factor(cell_type),
    experiment_id = factor(experiment_id)
  )

print(meta, n=72)


# Build count matrix from telescope_report.tsv
read_telescope <- function(sample_id, dir) {
  f <- file.path(dir, sample_id, paste0(sample_id, "-telescope_report.tsv"))
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  df <- read.table(f, header = TRUE, sep = "\t", comment.char = "#",
                   stringsAsFactors = FALSE)
  df <- df[, c("transcript", "final_count")]
  colnames(df)[2] <- sample_id 
  df
}

count_list <- lapply(meta$sample_id, read_telescope, dir = tele_dir)
names(count_list) <- meta$sample_id
count_list <- Filter(Negate(is.null), count_list)

count_matrix <- Reduce(function(a, b) full_join(a, b, by = "transcript"), count_list) %>%
  column_to_rownames("transcript") %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  as.matrix()
count_matrix <- round(count_matrix)

cat("Matrix:", nrow(count_matrix), "TEs x", ncol(count_matrix), "samples\n")


run_deseq2_by_exp <- function(exp_id) {
  m <- meta %>%
    filter(experiment_id == exp_id,
           !is.na(genotype),
           sample_id %in% colnames(count_matrix))
  
  genotypes <- unique(as.character(m$genotype))
  if (!"WT" %in% genotypes)        { message("Skip ", exp_id, ": no WT");        return(NULL) }
  if (n_distinct(m$genotype) < 2)  { message("Skip ", exp_id, ": one genotype"); return(NULL) }
  if (nrow(m) < 4)                 { message("Skip ", exp_id, ": n=", nrow(m));  return(NULL) }
  
  m$genotype_raw  <- as.character(m$genotype)
  m$genotype      <- factor(make.names(m$genotype_raw))
  geno_map        <- setNames(m$genotype_raw, as.character(m$genotype))  
  m$genotype      <- relevel(droplevels(m$genotype), ref = "WT")
  
  if (any(table(m$genotype) < 2)) { message("Skip ", exp_id, ": <2 reps/group"); return(NULL) }
  
  mat <- count_matrix[, m$sample_id, drop = FALSE]
  cd  <- m %>% column_to_rownames("sample_id")
  
  has_ct <- n_distinct(droplevels(factor(cd$cell_type))) > 1 &&
    all(rowSums(table(cd$cell_type, cd$genotype) > 0) > 1)
  design <- if (has_ct) ~ cell_type + genotype else ~ genotype
  
  message("\n▶ ", exp_id, " | n=", nrow(m),
          " | design: ", paste(deparse(design), collapse = ""),
          " | groups: ", paste(geno_map[levels(m$genotype)], collapse = ", "))
  
  dds <- DESeqDataSetFromMatrix(countData = mat, colData = cd, design = design)
  
  min_reps <- min(table(m$genotype))
  dds <- dds[rowSums(counts(dds) >= 5) >= min_reps, ]
  dds <- DESeq(dds, quiet = TRUE)
  
  non_wt <- setdiff(levels(m$genotype), "WT")
  
  res_list <- lapply(non_wt, function(g) {
    coef_name <- paste0("genotype_", g, "_vs_WT")
    res <- if (coef_name %in% resultsNames(dds)) {
      tryCatch(lfcShrink(dds, coef = coef_name, type = "normal"),
               error = function(e) results(dds, name = coef_name))
    } else {
      results(dds, contrast = c("genotype", g, "WT"),
              alpha = 0.05, pAdjustMethod = "BH")
    }
    res %>%
      as.data.frame() %>%
      rownames_to_column("transcript") %>%
      mutate(experiment_id = exp_id,
             comparison    = paste0(geno_map[[g]], "_vs_WT"))
  })
  
  list(res = bind_rows(res_list), dds = dds)
}


# Run
experiments <- unique(as.character(meta$experiment_id))
all_results <- lapply(experiments, run_deseq2_by_exp)
names(all_results) <- experiments

# Significant summary
sig_all <- bind_rows(lapply(all_results, function(x) {
  if (is.null(x)) return(NULL)
  x$res %>% filter(!is.na(padj), padj < 0.05)
})) %>% arrange(padj)

as.data.frame(sig_all %>% dplyr::count(experiment_id, comparison, name = "sig_TEs"))

head(sig_all$transcript, 20)

all <- bind_rows(lapply(all_results, function(x) {
  if (is.null(x)) return(NULL)
  x$res
}))

write.csv(all, "/tank/projects/ekashuk/TE/deseq2_by_experiment_homo_sapiens.csv", row.names = FALSE)
write.csv(sig_all, "/tank/projects/ekashuk/TE/sig_deseq2_by_experiment_mus_musculus.csv", row.names = FALSE)


## Homo sapiens ----
library(metafor); library(dplyr); library(tidyr); library(tibble)
library(stringr); library(forcats); library(ggplot2); library(pheatmap)

select <- dplyr::select; filter <- dplyr::filter
mutate <- dplyr::mutate; rename <- dplyr::rename; count <- dplyr::count

theme_custom <- function(base = 12) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))
}

COLORS <- c(ERV = "yellowgreen", LINE = "salmon", SINE = "mediumorchid2",
            DNA = "orange", SVA = "steelblue", Other = "grey60")

# Load data
all      <- read.csv('/tank/projects/ekashuk/TE/deseq2_by_experiment_homo_sapiens.csv')
sig_all  <- read.csv('/tank/projects/ekashuk/TE/sig_deseq2_by_experiment_homo_sapiens.csv')

samples <- as.data.frame(nanoparquet::read_parquet(
  "/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/homo_sapiens/samples.parquet"))
s2e     <- as.data.frame(nanoparquet::read_parquet(
  "/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/homo_sapiens/samples_to_experiment.parquet"))

meta_human <- samples %>%
  left_join(s2e, by = "sample_id") %>%
  filter(organism == "Homo sapiens") %>%
  mutate(genotype      = factor(genotype, levels = c("WT","SIRT6-KO","SIRT6-OE/K3Q")),
         cell_type     = factor(cell_type),
         experiment_id = factor(experiment_id))

# Classification
classify_te <- function(df) {
  df %>%
    mutate(
      family = str_extract(transcript, "^[^_]+"),
      class  = case_when(
        str_detect(family, "^HERV|^HML|^HUER|^ERV|^MER|^MLT|^THE|^MST|^PRIMA|^HARLEQUIN|^LTR") ~ "ERV",
        str_detect(family, "^L1|^LINE|^L2")                                                     ~ "LINE",
        str_detect(family, "^Alu|^SINE|^MIR|^B2")                                               ~ "SINE",
        str_detect(family, "^SVA")                                                              ~ "SVA",
        str_detect(family, "^DNA|^hAT|^Charlie|^Mariner|^Tc|^piggyBac|^Harbinger|^TcMar")       ~ "DNA",
        TRUE ~ "Other"
      ),
      direction = ifelse(log2FoldChange > 0, "up", "down"),
      contrast  = ifelse(str_detect(comparison, "KO"), "SIRT6-KO", "SIRT6-OE/K3Q")
    )
}

all     <- classify_te(all)
sig_all <- classify_te(sig_all)

# Descriptive part
family_updown <- sig_all %>%
  filter(contrast == "SIRT6-KO",
         !is.na(log2FoldChange), !is.na(padj)) %>%
  group_by(family, class) %>%
  summarise(
    n_up        = sum(log2FoldChange > 0, na.rm = TRUE),
    n_down      = sum(log2FoldChange < 0, na.rm = TRUE),
    n_exp_up    = n_distinct(experiment_id[log2FoldChange > 0]),
    n_exp_down  = n_distinct(experiment_id[log2FoldChange < 0]),
    mean_lfc_up = mean(log2FoldChange[log2FoldChange > 0], na.rm = TRUE),
    mean_lfc_dn = mean(log2FoldChange[log2FoldChange < 0], na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(
    n_total  = n_up + n_down,
    pct_up   = n_up / n_total,
    net_loci = n_up - n_down
  ) %>%
  filter(n_total >= 5) %>%
  slice_max(n_total, n = 30)

COLORS <- c(
  ERV   = "#4C9F70",
  LINE  = "#E07A5F",
  SINE  = "#9B6FB0",
  DNA   = "#E8A857",
  SVA   = "#3D6B99",
  Other = "#B5B5B5"
)

class_order <- c("ERV", "LINE", "SINE", "SVA", "DNA", "Other")

family_updown %>%
  mutate(class = factor(class, levels = class_order)) %>%
  ggplot(aes(fct_reorder(family, pct_up), pct_up, fill = class)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.3) +
  geom_text(
    aes(
      label = scales::percent(pct_up, accuracy = 1),
      hjust = ifelse(pct_up >= 0.5, -0.15, 1.15)
    ),
    size = 3.8, color = "black", fontface = "bold"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 0.4, alpha = 0.7) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.1), expand = c(0, 0),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_fill_manual(values = COLORS, breaks = class_order, name = NULL) +
  labs(x = NULL, y = "% upregulated over all DE TEs") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_line(color = "grey90", linewidth = 0.3),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position       = "top",
    legend.text          = element_text(size = 15, color = "black"),
    legend.key.size       = unit(0.6, "cm"),
    axis.text.y          = element_text(size = 14, color = "black"),
    axis.text.x          = element_text(size = 14, color = "black"),
    axis.title.x          = element_text(size = 14, color = "black", margin = margin(t = 8)),
    plot.margin          = margin(12, 40, 10, 10)
  )


## Meta for proportions

es_human  <- escalc(measure = "PLO", xi = n_up, ni = n_total, data = family_updown)
res_class <- rma(yi, vi, mods = ~ class, data = es_human)
res_class
predict(res_class, transf = transf.ilogit)



# Human

make_forest_plot <- function(es, res, class_colors, class_order,
                             title = NULL,
                             x_lab = "Log-odds (up vs down), SIRT6-KO") {
  
  z <- qnorm(0.975)
  
  df <- es %>%
    mutate(
      ci.lb  = yi - z * sqrt(vi),
      ci.ub  = yi + z * sqrt(vi),
      class  = factor(class, levels = class_order),
      weight = 1 / vi,
      sig    = ci.lb > 0 | ci.ub < 0
    ) %>%
    arrange(class, yi) %>%
    mutate(row_id = row_number())
  
  n <- nrow(df)
  
  stripes <- df %>%
    group_by(class) %>%
    summarise(ymin = min(row_id) - 0.5, ymax = max(row_id) + 0.5, .groups = "drop")
  
  
  line_y      <- 0.45
  zone_bottom <- -0.9
  pooled_row  <- (line_y + zone_bottom) / 2
  half_height <- 0.4
  
  pooled <- data.frame(
    x = c(res$ci.lb, res$b[1], res$ci.ub, res$b[1]),
    y = pooled_row + c(0, half_height, 0, -half_height)
  )
  
  
  label_df <- df %>%
    mutate(label = sprintf("%6.2f [%6.2f, %6.2f]", yi, ci.lb, ci.ub))
  
  sig_df <- df %>% filter(sig)
  
  max_abs  <- max(abs(c(df$ci.lb, df$ci.ub, res$ci.lb, res$ci.ub)))
  x_breaks <- pretty(c(-max_abs, max_abs))
  x_lim_lo <- min(x_breaks) * 1.05
  
  family_x <- x_lim_lo - max_abs * 0.55
  label_x  <- max_abs * 1.25
  star_x   <- label_x + max_abs * 1.6
  xlim_hi  <- star_x + max_abs * 0.4
  
  model_text <- sprintf(
    "Random-effects model (REML): k = %d   I\u00b2 = %.1f%%   \u03c4\u00b2 = %.4f\nQ(%d) = %.1f, p %s   |   Pooled: %.2f [%.2f, %.2f]",
    res$k, res$I2, res$tau2, res$k - 1, res$QE,
    ifelse(res$QEp < 0.0001, "< .0001", paste0("= ", round(res$QEp, 4))),
    res$b[1], res$ci.lb, res$ci.ub
  )
  
  ggplot() +
    geom_rect(data = stripes,
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = class),
              alpha = 0.18) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey30", linewidth = 0.5) +
    geom_segment(data = df, aes(y = row_id, yend = row_id, x = ci.lb, xend = ci.ub, color = class),
                 linewidth = 0.4) +
    geom_point(data = df, aes(y = row_id, x = yi, size = weight, color = class), shape = 15) +
    geom_polygon(data = pooled, aes(x = x, y = y), fill = "firebrick", color = "black", linewidth = 0.3) +
    geom_text(data = df, aes(y = row_id, x = family_x, label = family),
              hjust = 0, size = 3.8, fontface = "bold", color = "black") +
    geom_text(data = label_df, aes(y = row_id, x = label_x, label = label),
              hjust = 0, size = 3.3, fontface = "bold", color = "black", family = "mono") +
    geom_text(data = sig_df, aes(y = row_id, x = star_x), label = "*",
              size = 7, fontface = "bold", color = "black") +
    annotate("text", x = label_x, y = n + 1.2,
             label = "Estimate [95% CI]", hjust = 0, fontface = "bold", size = 3.6, color = "black") +
    annotate("text", x = family_x, y = n + 1.2,
             label = "Family", hjust = 0, fontface = "bold", size = 3.6, color = "black") +
    annotate("text", x = star_x, y = n + 1.2,
             label = "Sig.", hjust = 0, fontface = "bold", size = 3.6, color = "black") +
    annotate("segment", x = family_x, xend = xlim_hi, y = line_y, yend = line_y,
             linewidth = 0.6, color = "black") +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values = class_colors) +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    coord_cartesian(
      xlim = c(family_x, xlim_hi),
      ylim = c(zone_bottom, n + 1.5),
      clip = "off"
    ) +
    labs(title = title, x = x_lab, y = NULL, caption = model_text) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid            = element_blank(),
      panel.border          = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.title             = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
      axis.text.y            = element_blank(),
      axis.ticks.y           = element_blank(),
      axis.text.x            = element_text(color = "black", size = 14, face = "bold"),
      axis.title.x           = element_text(color = "black", size = 13, face = "bold", margin = margin(t = 8)),
      plot.caption            = element_text(hjust = 0, face = "bold", size = 9.5,
                                             color = "black", margin = margin(t = 10),
                                             lineheight = 1.2),
      plot.caption.position   = "plot",
      legend.position         = "none",
      plot.margin = margin(5, 160, 5, 20)
    )
}
COLORS_HUMAN      <- c(ERV = "#4C9F70", LINE = "#E07A5F")
class_order_human <- c("ERV", "LINE")

res_human_plain <- rma(yi, vi, data = es_human)

p_forest_human <- make_forest_plot(
  es_human, res_human_plain, COLORS_HUMAN, class_order_human,
)
p_forest_human

row_height_in <- 0.22
extra_in      <- 1.8
n_rows_human  <- nrow(es_human)
fig_height_h  <- n_rows_human * row_height_in + extra_in

ggsave("forest_human.pdf", p_forest_human,
       width = 11, height = fig_height_h, device = cairo_pdf)

## Mus musculus ----
library(metafor); library(dplyr); library(tidyr); library(tibble)
library(stringr); library(forcats); library(ggplot2)
select <- dplyr::select; filter <- dplyr::filter
mutate <- dplyr::mutate; rename <- dplyr::rename; count <- dplyr::count

OUTLIERS <- c("GSE166840", "GSE168983")

# Data
all     <- read.csv('/tank/projects/ekashuk/TE/deseq2_by_experiment_mus_musculus.csv')
sig_all <- read.csv('/tank/projects/ekashuk/TE/sig_deseq2_by_experiment_mus_musculus.csv')

samples <- as.data.frame(nanoparquet::read_parquet(
  "/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/mus_musculus/samples.parquet"))
s2e     <- as.data.frame(nanoparquet::read_parquet(
  "/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/mus_musculus/samples_to_experiment.parquet"))

meta_mouse <- samples %>%
  left_join(s2e, by = "sample_id") %>%
  mutate(genotype      = factor(genotype, levels = c("WT", "SIRT6-KO", "SIRT6-OE/K3Q")),
         cell_type     = factor(cell_type),
         experiment_id = factor(experiment_id))

# Classification
classify_te <- function(df) {
  df %>%
    mutate(
      family = str_extract(transcript, "^[^_]+"),
      class  = case_when(
        str_detect(family, "^IAP|^RLTR|^RMER|^MMERGLN|^MERVL|^MERV|^MTEa|^MLTR|^MaLR|^ETn|^ORR|^MMETn|^MMVL|^MTD|^MTC|^MTA|^MTB|^MTE|^URR|^MMERVK|^MurERV|^ERVB|^MLT|^MT2|^MER") ~ "ERV",
        str_detect(family, "^L1|^LINE|^L1Md|^Lx|^L2|^L3|^IMPB")                                                                                                                 ~ "LINE",
        str_detect(family, "^B1|^B2|^B4|^ID|^Alu|^SINE|^B3|^Dip|^MIR|^RSINE|^PB1")                                                                                              ~ "SINE",
        str_detect(family, "^MurSat|^GSAT|^MSAT|^MMSAT|^BGLII|^RCHARR")                                                                                                         ~ "Satellite",
        str_detect(family, "^OldhAT|^hAT|^Mariner|^Tc|^piggyBac|^Charlie|^Harbinger")                                                                                           ~ "DNA",
        TRUE ~ "Other"
      ),
      direction = ifelse(log2FoldChange > 0, "up", "down"),
      contrast  = ifelse(str_detect(comparison, "KO"), "SIRT6-KO", "SIRT6-OE/K3Q")
    )
}

all     <- classify_te(all)
sig_all <- classify_te(sig_all)

# Descriptive part
family_updown <- sig_all %>%
  filter(contrast == "SIRT6-KO",
         !is.na(log2FoldChange), !is.na(padj)) %>%
  group_by(family, class) %>%
  summarise(
    n_up        = sum(direction == "up"),
    n_down      = sum(direction == "down"),
    n_exp_up    = n_distinct(experiment_id[direction == "up"]),
    n_exp_down  = n_distinct(experiment_id[direction == "down"]),
    mean_lfc_up = mean(log2FoldChange[direction == "up"], na.rm = TRUE),
    mean_lfc_dn = mean(log2FoldChange[direction == "down"], na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(
    n_total  = n_up + n_down,
    pct_up   = n_up / n_total,
    net_loci = n_up - n_down
  ) %>%
  filter(n_total >= 5) %>%
  slice_max(n_total, n = 30)

# Plot
COLORS <- c(
  ERV       = "#4C9F70",
  LINE      = "#E07A5F",
  SINE      = "#9B6FB0",
  Satellite = "#3D6B99",
  DNA       = "#E8A857",
  Other     = "#B5B5B5"
)

class_order <- c("ERV", "LINE", "SINE", "Satellite", "DNA", "Other")

family_updown %>%
  mutate(class = factor(class, levels = class_order)) %>%
  ggplot(aes(fct_reorder(family, pct_up), pct_up, fill = class)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.3) +
  geom_text(
    aes(
      label = scales::percent(pct_up, accuracy = 1),
      hjust = ifelse(pct_up >= 0.5, -0.15, 1.15)
    ),
    size = 2.5, color = "black", fontface = "bold"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 0.4, alpha = 0.7) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.05), expand = c(0, 0),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_fill_manual(values = COLORS, breaks = class_order, name = NULL) +
  labs(x = NULL, y = "% upregulated over all DE TEs") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_line(color = "grey90", linewidth = 0.3),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position       = "top",
    legend.text          = element_text(size = 15, color = "black"),
    legend.key.size       = unit(0.6, "cm"),
    axis.text.y          = element_text(size = 11, color = "black"),
    axis.text.x          = element_text(size = 14, color = "black"),
    axis.title.x          = element_text(size = 14, color = "black", margin = margin(t = 8)),
    plot.margin          = margin(12, 40, 10, 10)
  )


es_mouse  <- escalc(measure = "PLO", xi = n_up, ni = n_total, data = family_updown)
res_class <- rma(yi, vi, mods = ~ class, data = es_mouse)
res_class
predict(res_class, transf = transf.ilogit)  

family_updown %>% count(class)
family_updown %>% filter(class == "Satellite") %>% select(family, n_up, n_down, n_total, pct_up)

forest(res_class, slab = paste(family_updown$family, family_updown$class))




library(dplyr); library(forcats); library(ggplot2); library(scales)

make_forest_plot <- function(es, res, class_colors, class_order,
                             title = NULL,
                             x_lab = "Log-odds (up vs down), SIRT6-KO") {
  
  z <- qnorm(0.975)
  
  df <- es %>%
    mutate(
      ci.lb  = yi - z * sqrt(vi),
      ci.ub  = yi + z * sqrt(vi),
      class  = factor(class, levels = class_order),
      weight = 1 / vi,
      sig    = ci.lb > 0 | ci.ub < 0
    ) %>%
    arrange(class, yi) %>%
    mutate(row_id = row_number())
  
  n <- nrow(df)
  
  stripes <- df %>%
    group_by(class) %>%
    summarise(ymin = min(row_id) - 0.5, ymax = max(row_id) + 0.5, .groups = "drop")
  
  
  line_y      <- 0.45
  zone_bottom <- -0.9
  pooled_row  <- (line_y + zone_bottom) / 2
  half_height <- 0.4
  
  pooled <- data.frame(
    x = c(res$ci.lb, res$b[1], res$ci.ub, res$b[1]),
    y = pooled_row + c(0, half_height, 0, -half_height)
  )
  
  
  label_df <- df %>%
    mutate(label = sprintf("%6.2f [%6.2f, %6.2f]", yi, ci.lb, ci.ub))
  
  sig_df <- df %>% filter(sig)
  
  max_abs  <- max(abs(c(df$ci.lb, df$ci.ub, res$ci.lb, res$ci.ub)))
  x_breaks <- pretty(c(-max_abs, max_abs))
  x_lim_lo <- min(x_breaks) * 1.05
  
  family_x <- x_lim_lo - max_abs * 0.55
  label_x  <- max_abs * 1.25
  star_x   <- label_x + max_abs * 1.6
  xlim_hi  <- star_x + max_abs * 0.4
  
  model_text <- sprintf(
    "Random-effects model (REML): k = %d   I\u00b2 = %.1f%%   \u03c4\u00b2 = %.4f\nQ(%d) = %.1f, p %s   |   Pooled: %.2f [%.2f, %.2f]",
    res$k, res$I2, res$tau2, res$k - 1, res$QE,
    ifelse(res$QEp < 0.0001, "< .0001", paste0("= ", round(res$QEp, 4))),
    res$b[1], res$ci.lb, res$ci.ub
  )
  
  ggplot() +
    geom_rect(data = stripes,
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = class),
              alpha = 0.18) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey30", linewidth = 0.5) +
    geom_segment(data = df, aes(y = row_id, yend = row_id, x = ci.lb, xend = ci.ub, color = class),
                 linewidth = 0.4) +
    geom_point(data = df, aes(y = row_id, x = yi, size = weight, color = class), shape = 15) +
    geom_polygon(data = pooled, aes(x = x, y = y), fill = "firebrick", color = "black", linewidth = 0.3) +
    geom_text(data = df, aes(y = row_id, x = family_x, label = family),
              hjust = 0, size = 3.8, fontface = "bold", color = "black") +
    geom_text(data = label_df, aes(y = row_id, x = label_x, label = label),
              hjust = 0, size = 3.3, fontface = "bold", color = "black", family = "mono") +
    geom_text(data = sig_df, aes(y = row_id, x = star_x), label = "*",
              size = 7, fontface = "bold", color = "black") +
    annotate("text", x = label_x, y = n + 1.2,
             label = "Estimate [95% CI]", hjust = 0, fontface = "bold", size = 3.6, color = "black") +
    annotate("text", x = family_x, y = n + 1.2,
             label = "Family", hjust = 0, fontface = "bold", size = 3.6, color = "black") +
    annotate("text", x = star_x, y = n + 1.2,
             label = "Sig.", hjust = 0, fontface = "bold", size = 3.6, color = "black") +
    annotate("segment", x = family_x, xend = xlim_hi, y = line_y, yend = line_y,
             linewidth = 0.6, color = "black") +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values = class_colors) +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    coord_cartesian(
      xlim = c(family_x, xlim_hi),
      ylim = c(zone_bottom, n + 1.5),
      clip = "off"
    ) +
    labs(title = title, x = x_lab, y = NULL, caption = model_text) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid            = element_blank(),
      panel.border          = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.title             = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
      axis.text.y            = element_blank(),
      axis.ticks.y           = element_blank(),
      axis.text.x            = element_text(color = "black", size = 14, face = "bold"),
      axis.title.x           = element_text(color = "black", size = 13, face = "bold", margin = margin(t = 8)),
      plot.caption            = element_text(hjust = 0, face = "bold", size = 9.5,
                                             color = "black", margin = margin(t = 10),
                                             lineheight = 1.2),
      plot.caption.position   = "plot",
      legend.position         = "none",
      plot.margin = margin(5, 160, 5, 20)
    )
}






# Mouse
COLORS_MOUSE <- c(ERV = "#4C9F70", LINE = "#E07A5F", SINE = "#9B6FB0", Satellite = "#3D6B99")
class_order_mouse <- c("ERV", "LINE", "SINE", "Satellite")

res_mouse_plain <- rma(yi, vi, data = es_mouse)

p_forest_mouse <- make_forest_plot(
  es_mouse, res_mouse_plain, COLORS_MOUSE, class_order_mouse,
)

p_forest_mouse
n_rows <- nrow(es_mouse)

row_height_in <- 0.22
extra_in      <- 1.8
fig_height    <- n_rows * row_height_in + extra_in

ggsave("forest_mouse.pdf", p_forest_mouse,
       width = 11, height = fig_height, device = cairo_pdf)




plot_forest_te <- function(es, class_colors, class_order, title = NULL,
                           x_lab = "Log-odds (up vs down), SIRT6-KO",
                           xlim = NULL, alim = NULL,
                           text_cex = 1) {
  
  te_data <- es %>%
    mutate(class = factor(class, levels = class_order)) %>%
    arrange(class, yi)
  
  res <- rma(
    yi = yi, vi = vi, data = te_data,
    method = "REML",
    control = list(iter.max = 1000, rel.tol = 1e-8)
  )
  
  family_width <- max(nchar(te_data$family), na.rm = TRUE)
  labels <- sprintf(paste0("%-", family_width, "s"), te_data$family)
  header <- sprintf(paste0("%-", family_width, "s"), "Family")
  
  class_bg <- setNames(
    adjustcolor(unname(class_colors[class_order]), alpha.f = 0.25),
    class_order
  )
  
  if (is.null(alim)) {
    max_abs <- max(abs(c(te_data$yi - 1.96 * sqrt(te_data$vi),
                         te_data$yi + 1.96 * sqrt(te_data$vi))))
    alim <- c(-max_abs * 1.1, max_abs * 1.1)
  }
  if (is.null(xlim)) xlim <- c(alim[1] * 3.2, alim[2] * 1.3)
  
  x_breaks <- pretty(alim)
  x_breaks <- x_breaks[abs(x_breaks) <= max(abs(alim)) * 1.001]
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
  
  par(mar = c(4, 4, 5, 1))
  par(family = "mono", font.main = 1)
  
  forest(
    res, slab = labels, xlab = x_lab, main = title,
    cex.main = 1.6 * text_cex, refline = 0,
    cex = 1 * text_cex, cex.lab = 1 * text_cex, cex.axis = 1 * text_cex,
    header = header, mlab = "", alim = alim, xlim = xlim, at = x_breaks,
    annotate = FALSE
  )
  
  k <- nrow(te_data)
  usr <- par("usr")
  
  for (i in 1:k) {
    row_y <- k - i + 1
    rect(xleft = usr[1], ybottom = row_y - 0.5,
         xright = usr[2], ytop = row_y + 0.5,
         col = class_bg[as.character(te_data$class[i])], border = NA)
  }
  
  row_diamond <- -1
  half_height <- 0.55
  polygon(
    x = c(res$ci.lb, res$b[1], res$ci.ub, res$b[1]),
    y = row_diamond + c(0, half_height, 0, -half_height),
    col = "firebrick", border = "black", lwd = 1.3
  )
  
  tau2 <- res$tau2
  I2   <- res$I2
  base_y <- usr[3] + 0.3
  
  text(x = usr[1], y = base_y + 0.9, pos = 4, cex = 1 * text_cex, "Random-effects model:")
  text(x = usr[1], y = base_y, pos = 4, cex = 1 * text_cex,
       bquote(I^2 ~ "=" ~ .(formatC(I2, digits = 1, format = "f")) * "%," ~
                tau^2 ~ "=" ~ .(formatC(tau2, digits = 4, format = "f")) * "," ~
                "k =" ~ .(res$k)))
  
  par(mar = c(4, 0, 5, 1), family = "mono")
  plot.new()
  legend("left", legend = class_order, fill = class_bg[class_order],
         border = NA, bty = "n", cex = 1.1 * text_cex, title = "TE class")
  
  layout(1)
  invisible(res)
}

# Human
COLORS_HUMAN      <- c(ERV = "blue1", LINE = "goldenrod2")
class_order_human <- c("ERV", "LINE")

plot_forest_te(es_human, COLORS_HUMAN, class_order_human)

pdf("forest_human.pdf", width = 8, height = 7, family = "mono")
plot_forest_te(es_human, COLORS_HUMAN, class_order_human)
dev.off()

# Mouse
COLORS_MOUSE      <- c(ERV = "blue1", LINE = "goldenrod2", SINE = "darkmagenta", Satellite = "slategrey")
class_order_mouse <- c("ERV", "LINE", "SINE", "Satellite")

pdf("forest_mouse.pdf", width = 8, height = 7, family = "mono")
plot_forest_te(es_mouse, COLORS_MOUSE, class_order_mouse, 
               text_cex = 0.75)
dev.off()

## Rattus norvegicus ----

# Data
sig_rat <- read.csv('/tank/projects/ekashuk/TE/sig_deseq2_by_experiment_rattus_norvegicus.csv')

samples_rat <- as.data.frame(nanoparquet::read_parquet(
  "/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/rattus_norvegicus/samples.parquet"))
s2e_rat     <- as.data.frame(nanoparquet::read_parquet(
  "/tank/projects/public_data/Sirt6_datasets/Expression/SIRT6_db/metadata/rattus_norvegicus/samples_to_experiment.parquet"))

meta_rat <- samples_rat %>%
  left_join(s2e_rat, by = "sample_id") %>%
  mutate(genotype      = factor(genotype, levels = c("WT","SIRT6-KO","SIRT6-OE/K3Q")),
         cell_type     = factor(cell_type),
         experiment_id = factor(experiment_id))

# Classification; rat IDs: Family_chrN_start_end, take everything before _chr
sig_rat <- sig_rat %>%
  mutate(
    family = {
      fam <- str_extract(transcript, "^(.+?)(?=_chr)")
      ifelse(is.na(fam), str_extract(transcript, "^[^_]+"), fam)
    },
    class = case_when(
      str_detect(family, "^IAP|^RLTR|^RMER|^MMERGLN|^MERVL|^MERV|^MTEa|^MLTR|^MaLR|^ETn|^ORR|^MMETn|^MMVL|^MTD|^MTC|^MTA|^MTB|^MTE|^URR|^MMERVK|^MurERV|^ERVB|^MLT|^MT2|^LTR|^RNLTR|^RnERV|^THE|^MER|^MYSERV|^ZP3AR|^ERVL|^ERV3|^RNERVK|^RAL_Rn|^EUTREP") ~ "ERV",
      str_detect(family, "^L1|^LINE|^Lx|^L2|^L3|^IMPB|^L1_Rat|^L1_Rn|^MamRTE")                                                                                                                                                                                    ~ "LINE",
      str_detect(family, "^B1|^B2|^B3|^B4|^ID|^Alu|^SINE|^MIR|^RSINE|^PB1|^BC1|^AmnSINE|^RNHAL")                                                                                                                                                                  ~ "SINE",
      str_detect(family, "^RatSat|^MurSat|^GSAT|^MSAT|^MMSAT|^BGLII|^RCHARR|^RNSAT|^RatRep")                                                                                                                                                                       ~ "Satellite",
      str_detect(family, "^Charlie|^Harbinger|^hAT|^OldhAT|^X2a_DNA|^Mariner|^Tc|^piggyBac|^TcMar|^DNA|^Tigger|^BLACKJACK")                                                                                                                                         ~ "DNA",
      str_detect(family, "^\\(|^GA-rich|^AT-rich|^Simple|^Low|^A-rich|^G-rich|^T-rich|^C-rich|^tRNA|^U2|^snRNA|^rRNA") ~ "Simple_repeat",
      TRUE ~ "Other"
    ),
    direction = ifelse(log2FoldChange > 0, "up", "down"),
    contrast  = ifelse(str_detect(comparison, "KO"), "SIRT6-KO", "SIRT6-OE/K3Q")
  )

# Check Other
sig_rat %>% filter(class == "Other") %>%
  count(family, sort = TRUE) %>% head(10) %>% as.data.frame()

# Barplot by class
sig_rat %>%
  filter(!class %in% c("Simple_repeat", "Other")) %>%
  count(contrast, class, direction) %>%
  ggplot(aes(class, n, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~ contrast, scales = "free_y") +
  scale_fill_manual(values = c(up = "#e74c3c", down = "#3498db")) +
  labs(title = "Rat TE classes (GSE276440)", x = NULL, y = "Count") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Top families: lollipop
top_rat <- sig_rat %>%
  filter(!class %in% c("Simple_repeat", "Other"),
         contrast == "SIRT6-KO") %>%
  group_by(family, class, direction) %>%
  summarise(mean_lfc = mean(log2FoldChange),
            n_loci   = n(), .groups = "drop") %>%
  filter(n_loci >= 3) %>%
  group_by(direction) %>%
  slice_max(abs(mean_lfc), n = 20) %>%
  ungroup()

top_rat %>%
  ggplot(aes(fct_reorder(family, mean_lfc), mean_lfc, color = class)) +
  geom_segment(aes(xend = family, y = 0, yend = mean_lfc), linewidth = 0.8) +
  geom_point(aes(size = n_loci)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = COLORS) +
  labs(title    = "Top TE families — SIRT6-KO rat (GSE276440)",
       subtitle = "≥3 significant loci",
       x = NULL, y = "Mean log2FoldChange", size = "n loci") +
  theme_custom()

# Forest plot with CI
sig_rat %>%
  filter(!class %in% c("Simple_repeat", "Other"),
         contrast == "SIRT6-KO") %>%
  group_by(family, class) %>%
  summarise(
    mean_lfc = mean(log2FoldChange),
    se       = sqrt(mean(lfcSE^2)),
    n_loci   = n(), .groups = "drop"
  ) %>%
  filter(n_loci >= 3) %>%
  mutate(ci_lb = mean_lfc - 1.96 * se,
         ci_ub = mean_lfc + 1.96 * se) %>%
  group_by(sign(mean_lfc)) %>%
  slice_max(abs(mean_lfc), n = 15) %>%
  ungroup() %>%
  ggplot(aes(fct_reorder(family, mean_lfc), mean_lfc, color = class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = ci_lb, ymax = ci_ub), width = 0.3, linewidth = 0.7) +
  geom_point(aes(size = n_loci)) +
  coord_flip() +
  scale_color_manual(values = COLORS) +
  scale_size_continuous(range = c(3, 8), name = "n loci") +
  labs(title    = "TE families — SIRT6-KO rat (GSE276440)",
       subtitle = "Mean log2FC ± 95% CI | ≥3 significant loci",
       x = NULL, y = "Mean log2FoldChange") +
  theme_custom()

sig_rat %>%
  filter(!class %in% c("Simple_repeat", "Other"),
         contrast == "SIRT6-KO") %>%
  group_by(family, class) %>%
  summarise(mean_lfc = mean(log2FoldChange), n_loci = n(), .groups = "drop") %>%
  filter(n_loci >= 3) %>%
  arrange(desc(mean_lfc)) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  as.data.frame() %>% print()

