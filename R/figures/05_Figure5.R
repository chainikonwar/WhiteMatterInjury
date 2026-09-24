###############################################################################
# Figure 5: culture-positive infection

###############################################################################


# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggtext)


# ---------------------------------------------------------------------------
# Function to create one CpG panel
# ---------------------------------------------------------------------------

make_cpg_panel <- function(cpg_id) {

  d_df <- data.frame(
    DNAme = betas_T1[cpg_id, ]
  )

  b_df <- plot_df[plot_df$CpG == cpg_id, ]

  ann <- annotations[annotations$CpG == cpg_id, ]


  # Density plot
  p_d <- ggplot(d_df, aes(x = DNAme)) +
    geom_density(
      fill = NA,
      color = "#444444",
      linewidth = 0.8
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    labs(
      x = "Beta value",
      y = "Density"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      axis.ticks.y = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.ticks.x = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      plot.margin = margin(2, 5, 2, 5)
    )


  # Boxplot + jitter
  p_b <- ggplot(
    b_df,
    aes(
      x = cult_pos_infection,
      y = DNAme,
      fill = cult_pos_infection
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7,
      width = 0.5
    ) +
    geom_jitter(
      width = 0.15,
      size = 1.2,
      alpha = 0.6,
      color = "black"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    scale_fill_manual(
      values = c(
        "No Infection" = "#e8e8e8",
        "Culture-Positive\nInfection" = "#777777"
      )
    ) +
    geom_text(
      data = ann,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "italic"
    ) +
    labs(
      x = NULL,
      y = "DNA Methylation (\u03b2-value)",
      fill = NULL,
      title = cpg_id
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_textbox_simple(
        face = "bold",
        halign = 0.5,
        size = 11,
        fill = "#d9d9d9",
        box.color = "black",
        linewidth = 0.5,
        padding = margin(4, 4, 4, 4),
        margin = margin(b = 4)
      ),
      axis.text.x = element_text(size = 9)
    )

  p_d / p_b + plot_layout(heights = c(0.5, 2))
}


# ---------------------------------------------------------------------------
# Build all panels
# ---------------------------------------------------------------------------

panels <- lapply(sig_cpgs, make_cpg_panel)

p_final <- wrap_plots(
  panels,
  ncol = 2
)

p_final


# ---------------------------------------------------------------------------
# Save Figure 5
# ---------------------------------------------------------------------------

ggsave(
  "figures/Figure5_CpG_infection_combined.tiff",
  plot = p_final,
  width = 170,
  height = 230,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)
