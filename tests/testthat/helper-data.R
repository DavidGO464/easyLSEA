# tests/testthat/helper-data.R

make_minimal_data <- function(seed = 42L) {
  set.seed(seed)
  names_pc <- paste0("PC ", c("34:1","36:2","38:4","16:0","20:4",
                              "18:1","22:6","34:2","36:4","38:2"))
  names_pe <- paste0("PE ", c("34:1","36:2","38:4","16:0","20:4",
                              "18:1","34:2","36:4"))
  names_tg <- paste0("TG(", c("54:3","52:2","56:4","50:1","54:6",
                              "52:3","54:1"), ")")
  names_sm <- paste0("SM d18:1/", c("16:0","18:0","20:0","22:0","24:0"))

  lip_names <- c(names_pc, names_pe, names_tg, names_sm)
  cls <- c(rep("PC", 10L), rep("PE", 8L), rep("TG", 7L), rep("SM", 5L))
  n <- length(lip_names)

  data.frame(
    LipidName        = lip_names,
    LipidClass       = cls,
    Shorthand        = lip_names,
    Confidence_rank  = "A",
    logFC            = c(rnorm(10L, 0.8), rnorm(8L, -0.5), rnorm(12L, 0)),
    P.Value          = runif(n, 1e-4, 1),
    stringsAsFactors = FALSE
  )
}

MINIMAL_DATA <- make_minimal_data()
ANNOTATED    <- suppressMessages(annotate_lipids(MINIMAL_DATA, verbose = FALSE))
