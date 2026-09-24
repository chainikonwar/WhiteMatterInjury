###############################################################################
# Figure 2: CMR-based EWAS volcano plot
###############################################################################

library(ggplot2)
library(scales)

# Load volcano plotting function
source("R/functions/Volcano_Plot_Function.R")

# ---------------------------------------------------------------------------
# Prepare CMR results
# ---------------------------------------------------------------------------

T1_CMR_WM_rlm <- T1_CMR_WM_rlm_all_T1_T2_site

dB_threshold <- 0.03

pval_threshold <- max(
  T1_CMR_WM_rlm$pvalue[T1_CMR_WM_rlm$FDR < 0.10],
  na.rm = TRUE
)

legend_title <- "DNAme – White Matter Injury Volume"

T1_CMR_WM_rlm <- T1_CMR_WM_rlm[
  c("CpG", "pvalue", "FDR", "DB")
]

colnames(T1_CMR_WM_rlm) <- c(
  "CpG",
  "pvalue",
  "qvalues",
  "deltabeta"
)

pvalue <- T1_CMR_WM_rlm$pvalue
deltabeta <- T1_CMR_WM_rlm$deltabeta


# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

p <- makeVolcano(
  pvalue,
  deltabeta,
  dB_threshold,
  pval_threshold,
  legend_title,
  plot_title = "Co-Methylated Region-based EWAS",
  DB_title = "Delta Beta",
  FDR_title = "BH FDR",
  dot_size = 1.5,
  x_lim = c(-0.15, 0.15),
  text_size = 8
)

p <- p +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.3, "cm")
  )

p


# ---------------------------------------------------------------------------
# Save figure
# ---------------------------------------------------------------------------

ggsave(
  "figures/Figure2_CMR_volcano.tiff",
  plot = p,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)
