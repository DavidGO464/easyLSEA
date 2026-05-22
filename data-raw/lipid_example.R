# data-raw/lipid_example.R
# Synthetic lipidomics dataset for easyLSEA examples and tests.
# Designed with known enrichment patterns to validate the package.

set.seed(2026L)

n <- 200L
classes    <- c("PC", "PE", "TG", "SM", "Cer", "LPC", "CE", "DG")
class_prob <- c(0.28, 0.22, 0.18, 0.10, 0.08, 0.06, 0.05, 0.03)
chains     <- c("34:1", "36:2", "38:4", "16:0", "18:1",
                "18:0", "20:4", "22:6", "16:1", "18:2")

cls <- sample(classes, n, replace = TRUE, prob = class_prob)
chn <- sample(chains,  n, replace = TRUE)

lipid_example <- data.frame(
  LipidName  = paste0(cls, " ", chn),
  LipidClass = cls,
  logFC      = c(
    rnorm(56,  mean =  0.9, sd = 1.1),  # PC/PE: enriched in case
    rnorm(36,  mean = -0.7, sd = 1.0),  # TG: depleted in case
    rnorm(108, mean =  0.0, sd = 1.0)   # others: no pattern
  ),
  P.Value    = runif(n, min = 1e-5, max = 1),
  stringsAsFactors = FALSE
)

lipid_example$adj.P.Val <- p.adjust(lipid_example$P.Value, method = "BH")
lipid_example$sig       <- as.integer(
  lipid_example$adj.P.Val < 0.05 &
    abs(lipid_example$logFC) > log2(1.25)
)

usethis::use_data(lipid_example, overwrite = TRUE, compress = "xz")

