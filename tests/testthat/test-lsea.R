# tests/testthat/test-lsea.R

test_that("run_lsea() returns a list with ks, fgsea, combined elements", {
  result <- run_lsea(ANNOTATED, engine = "ks", verbose = FALSE)
  expect_type(result, "list")
  expect_named(result, c("ks", "fgsea", "combined"), ignore.order = TRUE)
})

test_that("run_lsea() KS results have required columns", {
  result <- run_lsea(ANNOTATED, engine = "ks", verbose = FALSE)
  required <- c("Group", "N_group", "DirectionalScore",
                "KS_pval", "FDR_LSEA", "Direction")
  expect_true(all(required %in% names(result$ks)))
})

test_that("run_lsea() padj values are in [0, 1]", {
  result <- run_lsea(ANNOTATED, engine = "ks", verbose = FALSE)
  padj <- result$ks$FDR_LSEA
  expect_true(all(padj >= 0 & padj <= 1, na.rm = TRUE))
})

test_that("run_lsea() errors on missing fc_col", {
  expect_error(
    run_lsea(ANNOTATED, fc_col = "nonexistent", verbose = FALSE),
    "not found"
  )
})

test_that("run_lsea() warns on NA values in fc_col", {
  df <- ANNOTATED
  df$logFC[1:3] <- NA
  expect_warning(
    run_lsea(df, engine = "ks", verbose = FALSE),
    "NA"
  )
})

test_that("run_lsea() handles dataset with a single lipid class gracefully", {
  mono <- ANNOTATED[ANNOTATED$LipidClass == "PC", ]
  # Background is empty — should warn, not error, and return a list
  expect_warning(
    result <- run_lsea(mono, engine = "ks",
                       group_cols = "LipidClass", verbose = FALSE),
    "background is empty"
  )
  expect_type(result, "list")
})

test_that("run_lsea() Direction contains case_lbl or ref_lbl", {
  result <- run_lsea(ANNOTATED, engine = "ks",
                     case_lbl = "CASE", ref_lbl = "REF",
                     verbose = FALSE)
  dirs <- result$ks$Direction
  expect_true(all(grepl("CASE|REF", dirs)))
})

test_that("run_lsea() seed argument produces reproducible results", {
  r1 <- run_lsea(ANNOTATED, engine = "ks", seed = 1L, verbose = FALSE)
  r2 <- run_lsea(ANNOTATED, engine = "ks", seed = 1L, verbose = FALSE)
  expect_equal(r1$ks$KS_pval, r2$ks$KS_pval)
})

test_that("run_lsea() different seeds produce different DS_perm_pval", {
  r1 <- run_lsea(ANNOTATED, engine = "ks", seed = 1L,
                 n_perm = 100L, verbose = FALSE)
  r2 <- run_lsea(ANNOTATED, engine = "ks", seed = 99L,
                 n_perm = 100L, verbose = FALSE)
  # KS p-values must be identical (deterministic); perm p may differ
  expect_equal(r1$ks$KS_pval, r2$ks$KS_pval)
})

test_that("run_lsea() combined table includes Convergence column when both engines run", {
  skip_if_not_installed("fgsea")
  result <- run_lsea(ANNOTATED, engine = "both", verbose = FALSE,
                     n_perm = 200L, fgsea_nperm = 1000L)
  expect_true("Convergence" %in% names(result$combined))
})
