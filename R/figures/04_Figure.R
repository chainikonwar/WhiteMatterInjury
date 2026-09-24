###############################################################################
# Figure 4: Delta Beta effect sizes correlations
#
# T1 = Neonatal Period
# T2 = Term-Equivalent Age
# T3 = 36 Months
###############################################################################


# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(dplyr)


# ---------------------------------------------------------------------------
# Plot settings
# ---------------------------------------------------------------------------

lims <- c(-0.03, 0.03)
breaks <- seq(-0.03, 0.03, by = 0.01)


# ---------------------------------------------------------------------------
# Function to calculate Spearman correlation labels
# ---------------------------------------------------------------------------

fmt_label <- function(x, y, method = "spearman") {

  keep <- complete.cases(x, y)

  ct <- cor.test(
    x[keep],
    y[keep],
    method = method
  )

  rho <- round(unname(ct$estimate), 2)

  plabel <- if (ct$p.value < 0.001) {
    "p < 0.001"
  } else {
    paste0(
      "p = ",
      formatC(ct$p.value, digits = 2, format = "f")
    )
  }

  paste0("\u03c1 = ", rho, ", ", plabel)
}


# ---------------------------------------------------------------------------
# Correlation labels
# ---------------------------------------------------------------------------

lab_T1T2 <- fmt_label(
  Deltabeta_46_T1_T2_T3$DB_T1,
  Deltabeta_46_T1_T2_T3$DB_T2
)

lab_T1T3 <- fmt_label(
  Deltabeta_46_T1_T2_T3$DB_T1,
  Deltabeta_46_T1_T2_T3$DB_T3
)

lab_T2T3 <- fmt_label(
  Deltabeta_46_T1_T2_T3$DB_T2,
  Deltabeta_46_T1_T2_T3$DB_T3
)


# ---------------------------------------------------------------------------
# Panel A: T1 vs T2
# Neonatal Period vs Term-Equivalent Age
# ---------------------------------------------------------------------------

p1 <- ggplot(
  data.frame(
    T1 = Deltabeta_46_T1_T2_T3$DB_T1,
    T2 = Deltabeta_46_T1_T2_T3$DB_T2
  ),
  aes(x = T1, y = T2)
) +
  geom_point(
    size = 2.5,
    alpha = 0.8,
    color = "#08519c"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#08519c",
    fill = "#bdd7e7",
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.6
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
    label = lab_T1T2,
    hjust = 1.1,
    vjust = -0.5,
    size = 3.5,
    fontface = "italic",
    color = "#08519c"
  ) +
  scale_x_continuous(
    limits = lims,
    breaks = breaks
  ) +
  scale_y_continuous(
    limits = lims,
    breaks = breaks
  ) +
  labs(
    x = "\u0394\u03b2 at Neonatal Period",
    y = "\u0394\u03b2 at Term-Equivalent Age",
    title = "Neonatal Period vs Term-Equivalent Age"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    axis.title = element_text(size = 9)
  )


# ---------------------------------------------------------------------------
# Panel B: T1 vs T3
# Neonatal Period vs 36 Months
# ---------------------------------------------------------------------------

p2 <- ggplot(
  data.frame(
    T1 = Deltabeta_46_T1_T2_T3$DB_T1,
    T3 = Deltabeta_46_T1_T2_T3$DB_T3
  ),
  aes(x = T1, y = T3)
) +
  geom_point(
    size = 2.5,
    alpha = 0.8,
    color = "#69b1db"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#69b1db",
    fill = "#69b1db",
    alpha = 0.6
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
    label = lab_T1T3,
    hjust = 1.1,
    vjust = -0.5,
    size = 3.5,
    fontface = "italic",
    color = "#69b1db"
  ) +
  scale_x_continuous(
    limits = lims,
    breaks = breaks
  ) +
  scale_y_continuous(
    limits = lims,
    breaks = breaks
  ) +
  labs(
    x = "\u0394\u03b2 at Neonatal Period",
    y = "\u0394\u03b2 at 36 Months",
    title = "Neonatal Period vs 36 Months"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    axis.title = element_text(size = 9)
  )


# ---------------------------------------------------------------------------
# Panel C: T2 vs T3
# Term-Equivalent Age vs 36 Months
# ---------------------------------------------------------------------------

p3 <- ggplot(
  data.frame(
    T2 = Deltabeta_46_T1_T2_T3$DB_T2,
    T3 = Deltabeta_46_T1_T2_T3$DB_T3
  ),
  aes(x = T2, y = T3)
) +
  geom_point(
    size = 2.5,
    alpha = 0.8,
    color = "#8ca6b8"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#8ca6b8",
    fill = "#8ca6b8",
    alpha = 0.6
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
    label = lab_T2T3,
    hjust = 1.1,
    vjust = -0.5,
    size = 3.5,
    fontface = "italic",
    color = "#8ca6b8"
  ) +
  scale_x_continuous(
    limits = lims,
    breaks = breaks
  ) +
  scale_y_continuous(
    limits = lims,
    breaks = breaks
  ) +
  labs(
    x = "\u0394\u03b2 at Term-Equivalent Age",
    y = "\u0394\u03b2 at 36 Months",
    title = "Term-Equivalent Age vs 36 Months"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    axis.title = element_text(size = 9)
  )


# ---------------------------------------------------------------------------
# Combine panels
# ---------------------------------------------------------------------------

p_final <- (p1 | p2) /
           (p3 | plot_spacer())

p_final


# ---------------------------------------------------------------------------
# Save Figure 4
# ---------------------------------------------------------------------------

ggsave(
  "figures/Figure4_longitudinal_DeltaBeta_correlations.tiff",
  plot = p_final,
  width = 180,
  height = 160,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)
