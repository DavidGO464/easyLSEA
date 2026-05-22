# tests/testthat/test-easyLSEA.R

test_that("easyLSEA() returns an easyLSEA_result object", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     plots = FALSE, verbose = FALSE)
  expect_s3_class(result, "easyLSEA_result")
})

test_that("easyLSEA_result has the five expected slots", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     plots = FALSE, verbose = FALSE)
  expect_named(result, c("meta", "lsea", "chains", "plots", "input"),
               ignore.order = TRUE)
})

test_that("easyLSEA() $meta captures call metadata correctly", {
  result <- easyLSEA(MINIMAL_DATA,
                     case_lbl = "A", ref_lbl = "B",
                     engine = "ks", plots = FALSE, verbose = FALSE)
  expect_equal(result$meta$case_lbl, "A")
  expect_equal(result$meta$ref_lbl,  "B")
  expect_equal(result$meta$engine,   "ks")
  expect_equal(result$meta$n_lipids, nrow(MINIMAL_DATA))
})

test_that("easyLSEA() $input$data is annotated (has LipidClass column)", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     plots = FALSE, verbose = FALSE)
  expect_true("LipidClass" %in% names(result$input$data))
})

test_that("easyLSEA() runs chain analysis by default", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     plots = FALSE, verbose = FALSE)
  expect_false(is.null(result$chains))
  expect_named(result$chains, c("parsed", "summary"))
})

test_that("easyLSEA() skips chains when run_chains = FALSE", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     run_chains = FALSE, plots = FALSE, verbose = FALSE)
  expect_null(result$chains)
})

test_that("easyLSEA() output = 'separate' returns list with lsea and chains", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     output = "separate", plots = FALSE, verbose = FALSE)
  expect_type(result, "list")
  expect_named(result, c("lsea", "chains"), ignore.order = TRUE)
})

test_that("print.easyLSEA_result() prints without error", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     plots = FALSE, verbose = FALSE)
  expect_output(print(result), "easyLSEA result")
})

test_that("summary.easyLSEA_result() runs without error", {
  result <- easyLSEA(MINIMAL_DATA, engine = "ks",
                     plots = FALSE, verbose = FALSE)
  expect_output(summary(result), "easyLSEA summary")
})

test_that("easyLSEA() is reproducible with same seed", {
  r1 <- easyLSEA(MINIMAL_DATA, engine = "ks", seed = 7L,
                 plots = FALSE, verbose = FALSE)
  r2 <- easyLSEA(MINIMAL_DATA, engine = "ks", seed = 7L,
                 plots = FALSE, verbose = FALSE)
  expect_equal(r1$lsea$ks$KS_pval, r2$lsea$ks$KS_pval)
})

test_that("easyLSEA() works end-to-end with lipid_example dataset", {
  data("lipid_example", package = "easyLSEA")
  result <- suppressWarnings(suppressMessages(easyLSEA(
    data      = lipid_example,
    case_lbl  = "NASH",
    ref_lbl   = "Control",
    engine    = "ks",
    n_perm    = 200L,
    plots     = FALSE,
    verbose   = FALSE
  )))
  expect_s3_class(result, "easyLSEA_result")
  expect_true(nrow(result$lsea$ks) > 0L)
})
