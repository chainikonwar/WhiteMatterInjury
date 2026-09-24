###############################################################################
# Figure 1: Gestational age at birth by white matter injury group
###############################################################################

# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)


# ---------------------------------------------------------------------------
# Prepare metadata
# ---------------------------------------------------------------------------

meta <- SampleInfo_T1_T2_matched_n190 %>%
  mutate(
    WMI_3group = case_when(
      wm_injury == 0 ~ "None",
      wm_injury == 1 ~ "Mild",
      wm_injury %in% c(2, 3) ~ "Moderate_and_Severe"
    ),
    WMI_3group = factor(
      WMI_3group,
      levels = c("None", "Mild", "Moderate_and_Severe"),
      ordered = TRUE
    )
  )

meta_T1 <- subset(meta, Timepoint == "preterm")


# ---------------------------------------------------------------------------
# Figure 1
# ---------------------------------------------------------------------------

p1 <- ggplot(
  meta_T1,
  aes(
    x = WMI_3group,
    y = GA_birth,
    fill = WMI_3group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.5,
    linewidth = 0.7
  ) +
  geom_jitter(
    size = 2,
    shape = 21,
    alpha = 0.7,
    width = 0.15,
    color = "black",
    fill = "grey90"
  ) +
  scale_fill_manual(
    values = c(
      "None" = "#D3D3D3",
      "Mild" = "#808080",
      "Moderate_and_Severe" = "#636161"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "None" = "None",
      "Mild" = "Mild",
      "Moderate_and_Severe" = "Moderate &\nSevere"
    )
  ) +
  annotate(
    "text",
    x = 2,
    y = max(meta_T1$GA_birth, na.rm = TRUE) + 0.5,
    label = "H = 3.32, p = 0.19",
    size = 5,
    hjust = 0.5
  ) +
  labs(
    x = "White matter injury",
    y = "Gestational Age at Birth (weeks)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, vjust = -0.3),
    axis.title.y = element_text(size = 16, vjust = 2),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7)
  )

p1


# ---------------------------------------------------------------------------
# Save figure
# ---------------------------------------------------------------------------

ggsave(
  "Figure1_GA_WMI.tiff",
  plot = p1,
  width = 5,
  height = 5,
  dpi = 300,
  compression = "lzw"
)
