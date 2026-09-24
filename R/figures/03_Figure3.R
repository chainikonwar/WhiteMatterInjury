###############################################################################
# Figure 3: Replication analysis
# Scatter plot of effect sizes and forest plot across timepoints
# White Matter Injury Project
###############################################################################

library(ggplot2)
library(patchwork)

# ---------------------------------------------------------------------------
# Scatter plot
# ---------------------------------------------------------------------------

cor_spearman <- cor.test(
  T2_replication$T1_deltabeta,
  T2_replication$deltabeta,
  method = "spearman"
)

r_label <- sprintf(
  "\u03c1 = %.2f, p %s",
  cor_spearman$estimate,
  ifelse(
    cor_spearman$p.value < 0.001,
    "< 0.001",
    sprintf("= %.3f", cor_spearman$p.value)
  )
)

p_scatter <- ggplot(
  data.frame(
    T1 = T2_replication$T1_deltabeta,
    T2 = T2_replication$deltabeta
  ),
  aes(x = T1, y = T2)
) +
  geom_point(size = 2.5, alpha = 0.8, color = "#08519c") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#08519c",
    fill = "#bdd7e7",
    alpha = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.4
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = r_label,
    hjust = 1.1,
    vjust = -0.5,
    size = 3.5,
    fontface = "italic",
    color = "#08519c"
  ) +
  labs(
    x = "\u0394\u03b2 at Neonatal period",
    y = "\u0394\u03b2 at TEA",
    title = "Neonatal Period vs Term-Equivalent Age"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11)
  )

# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------

forest_df <- rbind(
  data.frame(
    ID = T2_replication$ID,
    coef = T2_replication$T1_deltabeta,
    se = T2_replication$T1_stderror * fixed_range,
    Timepoint = "Neonatal period"
  ),
  data.frame(
    ID = T2_replication$ID,
    coef = T2_replication$deltabeta,
    se = T2_replication$stderror * fixed_range,
    Timepoint = "Term-Equivalent Age"
  )
)

forest_df$lower <- forest_df$coef - 1.96 * forest_df$se
forest_df$upper <- forest_df$coef + 1.96 * forest_df$se

cpg_ids <- replication_targets$id[replication_targets$source == "CpG"]
cmr_ids <- replication_targets$id[replication_targets$source == "CMR"]

id_order <- c(rev(cpg_ids), rev(cmr_ids))
id_order_labeled <- ifelse(
  id_order %in% cmr_ids,
  paste0("* ", id_order),
  paste0("\u2020 ", id_order)
)

forest_df$ID_labeled <- ifelse(
  forest_df$ID %in% cmr_ids,
  paste0("* ", forest_df$ID),
  paste0("\u2020 ", forest_df$ID)
)

forest_df$ID_labeled <- factor(
  forest_df$ID_labeled,
  levels = id_order_labeled
)

n_single <- length(cpg_ids)

p_forest <- ggplot(
  forest_df,
  aes(x = coef, y = ID_labeled, color = Timepoint, shape = Timepoint)
) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0.3,
    position = position_dodge(width = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(
    yintercept = n_single + 0.5,
    linetype = "dotted",
    color = "grey70",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c(
      "Neonatal period" = "#08519c",
      "Term-Equivalent Age" = "#69b1db"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Neonatal period" = 16,
      "Term-Equivalent Age" = 15
    )
  ) +
  labs(
    x = "\u0394\u03b2",
    y = NULL,
    caption = "* Co-methylated region (CMR); \u2020 Single CpG",
    color = NULL,
    shape = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    plot.caption = element_text(size = 9, hjust = 0),
    plot.margin = margin(5, 5, 5, 10, "mm")
  )

# ---------------------------------------------------------------------------
# Combine panels
# ---------------------------------------------------------------------------

p_combined <- p_scatter / p_forest +
  plot_layout(heights = c(1.3, 2.5)) +
  plot_annotation(
    title = " ",
    theme = theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  )

p_combined

# ---------------------------------------------------------------------------
# Save figure
# ---------------------------------------------------------------------------

ggsave(
  "figures/Figure3_replication_scatter_forest.tiff",
  plot = p_combined,
  width = 180,
  height = 240,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)
