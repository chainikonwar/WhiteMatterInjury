###############################################################################
# CpG-level robust linear regression (same code used for CMR analysis)
# White Matter Injury Project
#
# Outcome: DNA methylation beta value
# Exposure: White matter injury volume (wmi_mm3)
# T1: refers to the neonatal period timepoint
#
# Model:
# DNAm ~ wmi_mm3 + Sex + PMA + GA_birth + Epi + Site
#
# Robust regression: robustbase::lmrob()
# Setting: KS2014
#
# Multiple-testing correction:
# Bacon correction followed by Benjamini-Hochberg FDR
###############################################################################


# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------

library(robustbase)
library(R.utils)
library(bacon)
library(ggplot2)


# ---------------------------------------------------------------------------
# Robust linear regression
# ---------------------------------------------------------------------------

T1_CpG_WM_lmrob_all <- as.data.frame(
  t(
    sapply(1:nrow(CpGs_T1), function(x) {

      CpG <- rownames(CpGs_T1)[x]

      if (x %% 5000 == 0) {
        cat("Row", x, "of", nrow(CpGs_T1), "\n")
      }

      fit <- tryCatch(
        withTimeout(
          lmrob(
            CpGs_T1[x, ] ~ wmi_mm3 + Sex + PMA + GA_birth + Epi + Site,
            data = T1_all,
            setting = "KS2014"
          ),
          timeout = 10,
          onTimeout = "error"
        ),
        error = function(e) {
          message("ERROR row ", x, ": ", conditionMessage(e))
          NULL
        }
      )

      if (is.null(fit) || !isTRUE(fit$converged)) {
        return(c(CpG, NA, NA, NA, NA))
      }

      s <- summary(fit)$coefficients

      c(
        CpG,
        s["wmi_mm3", "Pr(>|t|)"],
        s["wmi_mm3", "Estimate"],
        s["wmi_mm3", "t value"],
        s["wmi_mm3", "Std. Error"]
      )
    })
  )
)


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

colnames(T1_CpG_WM_lmrob_all) <-
  c("CpG", "pvalue", "coef", "tval", "stderror")

rownames(T1_CpG_WM_lmrob_all) <-
  T1_CpG_WM_lmrob_all$CpG

T1_CpG_WM_lmrob_all$pvalue <-
  as.numeric(as.character(T1_CpG_WM_lmrob_all$pvalue))

T1_CpG_WM_lmrob_all$coef <-
  as.numeric(as.character(T1_CpG_WM_lmrob_all$coef))

T1_CpG_WM_lmrob_all$tval <-
  as.numeric(as.character(T1_CpG_WM_lmrob_all$tval))

T1_CpG_WM_lmrob_all$stderror <-
  as.numeric(as.character(T1_CpG_WM_lmrob_all$stderror))


# ---------------------------------------------------------------------------
# Effect size (Delta beta)
# ---------------------------------------------------------------------------

quant_range <- function(y, q_lo, q_hi) {
  quantile(y, q_hi, na.rm = TRUE) -
    quantile(y, q_lo, na.rm = TRUE)
}

T1_CpG_WM_lmrob_all$DB <-
  T1_CpG_WM_lmrob_all$coef *
  quant_range(T1_all$wmi_mm3, 0.05, 0.95)


# ---------------------------------------------------------------------------
# Remove non-converged models
# ---------------------------------------------------------------------------

sum(is.na(T1_CpG_WM_lmrob_all$pvalue))

T1_CpG_WM_lmrob_clean <-
  T1_CpG_WM_lmrob_all[
    !is.na(T1_CpG_WM_lmrob_all$pvalue),
  ]


# ---------------------------------------------------------------------------
# Bacon correction
# inflation(bc) gives you the empirical-null-based inflation estimate and bias(bc) tells you if your test statistics are systematically shifted
# ---------------------------------------------------------------------------

z <- T1_CpG_WM_lmrob_clean$tval

bc <- bacon(z)

inflation(bc)
bias(bc)

T1_CpG_WM_lmrob_clean$pvalue_corrected <-
  as.numeric(pval(bc, corrected = TRUE))

T1_CpG_WM_lmrob_clean$FDR_corrected <-
  p.adjust(
    T1_CpG_WM_lmrob_clean$pvalue_corrected,
    method = "BH"
  )


# ---------------------------------------------------------------------------
# Prepare results for plotting
# ---------------------------------------------------------------------------

T1_CpG_WM_lmrob_plot <-
  T1_CpG_WM_lmrob_clean[
    c(
      "CpG",
      "pvalue_corrected",
      "FDR_corrected",
      "DB"
    )
  ]

colnames(T1_CpG_WM_lmrob_plot) <-
  c(
    "CpG",
    "pvalue",
    "qvalues",
    "deltabeta"
  )


# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

dB_threshold <- 0.03 (|Δβ| of >0.03, exceeding the root mean squared error (RMSE) estimated from technical replicates following preprocessing (RMSE = 0.026)

pval_threshold <-
  max(
    T1_CpG_WM_lmrob_plot$pvalue[
      T1_CpG_WM_lmrob_plot$qvalues < 0.05
    ],
    na.rm = TRUE
  )

legend_title <- "WMI Volume (lmrob)"

pvalue <- T1_CpG_WM_lmrob_plot$pvalue
deltabeta <- T1_CpG_WM_lmrob_plot$deltabeta


# Volcano plot function
source("R/functions/Volcano_Plot_Function.R")


p <- makeVolcano(
  pvalue,
  deltabeta,
  dB_threshold,
  pval_threshold,
  legend_title,
  plot_title = "",
  DB_title = "Delta Beta",
  FDR_title = "BH FDR",
  dot_size = 1.5,
  x_lim = c(-0.12, 0.12),
  text_size = 20
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


ggsave(
  "figures/volcano_WMI_DNAme_lmrob_CpG.tiff",
  plot = p,
  width = 180,
  height = 140,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)


# ---------------------------------------------------------------------------
# Significant CpGs
# ---------------------------------------------------------------------------

T1_CpG_WM_lmrob_finalhits <-
  T1_CpG_WM_lmrob_plot[
    !is.na(T1_CpG_WM_lmrob_plot$qvalues) &
      T1_CpG_WM_lmrob_plot$qvalues < 0.05 &
      abs(T1_CpG_WM_lmrob_plot$deltabeta) >= 0.03,
  ]

dim(T1_CpG_WM_lmrob_finalhits)
