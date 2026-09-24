###############################################################################
# Volcano Plot Function
#
# Creates a volcano plot using:
#   - P-values
#   - Delta Beta effect sizes
#   - Delta Beta threshold
#   - Significance-derived P-value threshold
###############################################################################


makeVolcano <- function(
  pvalue,
  deltabeta,
  dB_threshold,
  pval_threshold,
  legend_title,
  DB_title,
  FDR_title,
  plot_title,
  x_lim = c(-0.15, 0.15),
  dot_size = 2,
  text_size = 30,
  increase_color = "#004c99",
  decrease_color = "#994d00",
  ylab = "P-value (-log10)"
) {

  # -------------------------------------------------------------------------
  # Create plotting dataframe
  # -------------------------------------------------------------------------

  volcano <- data.frame(
    Pvalue = pvalue,
    Delta_Beta = deltabeta
  )


  # -------------------------------------------------------------------------
  # Identify significant loci passing Delta Beta threshold
  # -------------------------------------------------------------------------

  sta_delbeta <- deltabeta[pvalue <= pval_threshold]

  sta_delbeta <- sta_delbeta[
    abs(sta_delbeta) >= dB_threshold
  ]

  print(
    paste(
      "Increase in DNAm",
      length(sta_delbeta[sta_delbeta >= dB_threshold]),
      sep = ": "
    )
  )

  print(
    paste(
      "Decrease in DNAm",
      length(sta_delbeta[sta_delbeta <= -dB_threshold]),
      sep = ": "
    )
  )


  # -------------------------------------------------------------------------
  # Classify points for plotting
  # -------------------------------------------------------------------------

  volcano$Significance <- sapply(
    seq_len(nrow(volcano)),
    function(x) {

      if (volcano$Pvalue[x] <= pval_threshold) {

        if (abs(volcano$Delta_Beta[x]) >= dB_threshold) {

          if (volcano$Delta_Beta[x] > 0) {
            "Increased DNAm\nWith Delta Beta"
          } else {
            "Decreased DNAm\nWith Delta Beta"
          }

        } else {

          if (volcano$Delta_Beta[x] > 0) {
            "Increased DNAm"
          } else {
            "Decreased DNAm"
          }
        }

      } else {

        "Not Significant"
      }
    }
  )


  # -------------------------------------------------------------------------
  # Thresholds
  # -------------------------------------------------------------------------

  DB_threshold <- data.frame(
    x = c(-dB_threshold, dB_threshold)
  )

  FDR_threshold <- data.frame(
    x = -log10(pval_threshold)
  )


  # -------------------------------------------------------------------------
  # Volcano plot
  # -------------------------------------------------------------------------

  volcano_plot <- ggplot(
    volcano,
    aes(
      x = Delta_Beta,
      y = -log10(Pvalue),
      color = Significance
    )
  ) +
    geom_point(
      shape = 19,
      size = dot_size
    ) +

    scale_color_manual(
      values = c(
        "Decreased DNAm\nWith Delta Beta" = decrease_color,
        "Increased DNAm\nWith Delta Beta" = increase_color,
        "Decreased DNAm" = scales::muted(
          decrease_color,
          l = 80,
          c = 30
        ),
        "Increased DNAm" = scales::muted(
          increase_color,
          l = 70,
          c = 40
        ),
        "Not Significant" = "grey"
      ),
      name = legend_title
    ) +

    geom_vline(
      aes(
        xintercept = x,
        linetype = DB_title
      ),
      data = DB_threshold
    ) +

    geom_hline(
      aes(
        yintercept = x,
        linetype = FDR_title
      ),
      data = FDR_threshold
    ) +

    scale_linetype_discrete(
      name = "Thresholds"
    ) +

    labs(
      x = "Delta Beta",
      y = ylab,
      title = plot_title
    ) +

    coord_cartesian(
      xlim = x_lim
    ) +

    theme_test(
      base_size = text_size
    ) +

    theme(
      legend.text = element_text(size = 20),
      axis.text = element_text(
        size = 20,
        color = "black"
      ),
      axis.title = element_text(size = 25),
      legend.title = element_text(size = 25),
      plot.title = element_text(
        hjust = 0.5
      )
    )


  # -------------------------------------------------------------------------
  # Return plot
  # -------------------------------------------------------------------------

  volcano_plot
}
