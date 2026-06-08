##################################################################################################
######################################## Forest plot for DEGs ####################################
##################################################################################################

source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

library(dplyr)
library(stringr)
library(ape)

meta_ko <- read.csv("~/SIRT6_db/meta_analysis/meta_ko.csv")

# 1. Phylogenetic matrix

# Upload Newick tree with divergence times
phylo_tree <- read.tree("~/SIRT6_db/meta_analysis/six_species_tree.tre")
# Clean up tip labels for display
phylo_tree$tip.label <- gsub("_", " ", phylo_tree$tip.label)
phylo_cor <- vcv(phylo_tree, corr = TRUE)

# 2. Create a function for constructing a forest plot
plot_forest_gene <- function(gene_name, data, phylo_cor) {
  
  gene_data <- data %>%
    filter(human_gene_symbol == gene_name) %>%
    arrange(organism, log2FoldChange)
  
  res <- rma.mv(
    yi = log2FoldChange,
    V  = lfcSE^2,
    random = ~ 1 | organism/experiment_id,
    R = list(organism = phylo_cor),
    data = gene_data,
    method = "REML",
    control = list(iter.max = 1000, rel.tol = 1e-8)
  )
  
  species_width <- max(nchar(gene_data$organism_short), na.rm = TRUE)
  exp_width     <- max(nchar(gene_data$experiment_id), na.rm = TRUE)
  
  labels <- sprintf(
    paste0("%-", species_width, "s | %-", exp_width, "s | %s"),
    gene_data$organism_short,
    gene_data$experiment_id,
    gene_data$biological_system
  )
  
  header <- sprintf(
    paste0("%-", species_width, "s | %-", exp_width, "s | %s"),
    "Species", "Experiment", "System"
  )
  
  # Background colors per organism (semi-transparent)
  organism_bg <- c(
    "Homo sapiens"             = adjustcolor("pink2", alpha.f = 0.25),
    "Mus musculus"             = adjustcolor("skyblue2", alpha.f = 0.25),
    "Rattus norvegicus"        = adjustcolor("lightgreen", alpha.f = 0.25),
    "Sus scrofa"               = adjustcolor("mediumpurple1", alpha.f = 0.25),
    "Macaca fascicularis"      = adjustcolor("coral2", alpha.f = 0.25),
    "Drosophila melanogaster"  = adjustcolor("#A65628", alpha.f = 0.25)
  )
  
  par(mar = c(4, 4, 5, 2))
  par(family = "mono")
  
  # Forest plot without shading
  forest(
    res,
    slab = labels,
    xlab = "log2 Fold Change (SIRT6 KO vs WT)",
    main = gene_name,
    cex.main = 3,
    refline = 0,
    cex = 1.9,
    cex.lab = 1.6,
    cex.axis = 1.5,
    header = header,
    mlab = "",
    alim = c(-4, 4),
    xlim = c(-19, 12)
  )
  
  # Draw colored background rectangles per organism
  k <- nrow(gene_data)
  usr <- par("usr")
  
  for (i in 1:k) {
    row_y <- k - i + 1
    rect(
      xleft   = usr[1],
      ybottom = row_y - 0.5,
      xright  = usr[2],
      ytop    = row_y + 0.5,
      col     = organism_bg[gene_data$organism[i]],
      border  = NA
    )
  }
  
  # Add meta-analysis diamond
  addpoly(res, row = -1, mlab = "", col = "firebrick", border = "firebrick")
  
  # Heterogeneity stats
  tau2_total <- sum(res$sigma2)
  tau2_phylo <- res$sigma2[1]
  tau2_experiment <- res$sigma2[2]
  I2 <- 100 * tau2_total / (tau2_total + mean(gene_data$lfcSE^2))
  
  line_height <- 1.1
  base_y <- usr[3] + 0.3
  
  text(x = usr[1], y = base_y + line_height, pos = 4, cex = 1.7,
       "Multilevel model:")
  
  text(x = usr[1], y = base_y, pos = 4, cex = 1.7,
       bquote(I^2 ~ "=" ~ .(formatC(I2, digits = 1, format = "f")) * "%," ~
                tau[phylo]^2 ~ "=" ~ .(formatC(tau2_phylo, digits = 3, format = "e")) * "," ~
                tau[exp]^2 ~ "=" ~ .(formatC(tau2_experiment, digits = 3, format = "e"))))
}

# 3. Save 
out_dir <- "~/SIRT6_db/meta_analysis/"

# APLN
png(
  filename = paste0(out_dir, "forest_APLN.png"),
  width = 10000,
  height = 7500,
  res = 600
)

par(family = "mono")

plot_forest_gene("APLN", meta_ko, phylo_cor)

dev.off()

# ZNF518A
png(
  filename = paste0(out_dir, "forest_ZNF518A.png"),
  width = 11000,
  height = 9000,
  res = 600
)

par(family = "mono")

plot_forest_gene("ZNF518A", meta_ko, phylo_cor)

dev.off()
