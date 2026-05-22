# tests/testthat/test-chains.R

# -- Low-level parser tests ---------------------------------------------------

test_that(".parse_one() extracts chain length and saturation", {
  result <- easyLSEA:::.parse_one("18:2")
  expect_equal(result[["cl"]], 18L)
  expect_equal(result[["cs"]], 2L)
})

test_that(".parse_one() strips modification suffixes before parsing", {
  # ;O2 = hydroxylation modifier in Liebisch notation
  result <- easyLSEA:::.parse_one("18:1;O2")
  expect_equal(result[["cl"]], 18L)
  expect_equal(result[["cs"]], 1L)
})

test_that(".parse_one() returns NA for non-parseable input", {
  result <- easyLSEA:::.parse_one("unknown")
  expect_true(is.na(result[["cl"]]))
  expect_true(is.na(result[["cs"]]))
})

test_that(".parse_sn2() extracts sn-2 chain from PC notation", {
  result <- easyLSEA:::.parse_sn2("PC 18:0/20:4")
  expect_false(is.null(result))
  expect_equal(result$sn2_cl, 20L)
  expect_equal(result$sn2_cs, 4L)
  expect_equal(result$chain_type, "sn2")
})

test_that(".parse_sn2() accepts underscore separator (Liebisch notation)", {
  result <- easyLSEA:::.parse_sn2("PE 18:0_18:2")
  expect_false(is.null(result))
  expect_equal(result$sn2_cl, 18L)
})

test_that(".parse_sn2() returns NULL for total notation (no separator)", {
  result <- easyLSEA:::.parse_sn2("PC 36:2")
  expect_null(result)
})

test_that(".parse_nacyl() extracts N-acyl chain from SM notation", {
  result <- easyLSEA:::.parse_nacyl("SM d18:1/16:0")
  expect_false(is.null(result))
  expect_equal(result$nacyl_cl, 16L)
  expect_equal(result$nacyl_cs, 0L)
  expect_equal(result$chain_type, "nacyl")
})

test_that(".parse_nacyl() accepts underscore separator for SM", {
  result <- easyLSEA:::.parse_nacyl("SM d18:1_16:0")
  expect_false(is.null(result))
  expect_equal(result$nacyl_cl, 16L)
})

test_that(".parse_pi() returns excluded status for total notation", {
  result <- easyLSEA:::.parse_pi("PI 38:4")
  expect_equal(result$status, "excluded_total_notation_unresolved_PI")
  expect_null(result$data)
})

test_that(".parse_pi() parses resolved PI with two rows", {
  result <- easyLSEA:::.parse_pi("PI 18:0/20:4")
  expect_equal(result$status, "parsed_PI_resolved")
  expect_equal(nrow(result$data), 2L)
  expect_equal(result$data$chain_type, c("long_format", "long_format"))
})

test_that(".parse_long() returns NULL for TG total notation", {
  expect_null(easyLSEA:::.parse_long("TAG 54:3"))
})

test_that(".parse_long() returns NULL for ether TG", {
  expect_null(easyLSEA:::.parse_long("TG O-54:3"))
})

test_that(".parse_long() returns 3 rows for resolved TG", {
  result <- easyLSEA:::.parse_long("TG(16:0/18:1/18:2)")
  expect_equal(nrow(result), 3L)
  expect_true(all(result$chain_type == "long_format"))
})

# -- parse_lipid_chains() integration tests -----------------------------------

test_that("parse_lipid_chains() returns a list with parsed and summary", {
  result <- parse_lipid_chains(ANNOTATED)
  expect_type(result, "list")
  expect_named(result, c("parsed", "summary"))
})

test_that("parse_lipid_chains() parsed data.frame has chain columns", {
  # Requires resolved-chain notation (with "/" separator) to produce parsed rows
  df_resolved <- data.frame(
    LipidName  = c("PC 18:0/20:4", "PE 16:0/18:2", "SM d18:1/16:0",
                   "PC 18:1/18:2", "PE 18:0/20:4"),
    LipidClass = c("PC", "PE", "SM", "PC", "PE"),
    logFC      = c(1.2, -0.5, 0.3, 0.8, -0.9),
    stringsAsFactors = FALSE
  )
  annotated_r <- annotate_lipids(df_resolved, verbose = FALSE)
  result <- parse_lipid_chains(annotated_r)
  expect_true(nrow(result$parsed) > 0L,
              label = "parsed must have rows for resolved-chain lipids")
  expect_true("analysis_chain_cl" %in% names(result$parsed))
  expect_true("analysis_chain_cs" %in% names(result$parsed))
  expect_true("chain_type"        %in% names(result$parsed))
})

test_that("parse_lipid_chains() summary has a status column", {
  result <- parse_lipid_chains(ANNOTATED)
  expect_true("status" %in% names(result$summary))
  valid_statuses <- c("parsed", "excluded_class", "excluded_sm_unresolved",
                      "excluded_unresolved", "excluded_rank_P_or_NA",
                      "excluded_total_notation_unresolved_PI",
                      "excluded_unresolved_PI", "parsed_PI_resolved")
  expect_true(all(result$summary$status %in% valid_statuses))
})

test_that("parse_lipid_chains() parsed output has no negative chain lengths", {
  result <- parse_lipid_chains(ANNOTATED)
  cl <- result$parsed$analysis_chain_cl
  expect_true(all(cl >= 0L, na.rm = TRUE))
})

test_that("parse_lipid_chains() errors on missing class column", {
  df <- ANNOTATED
  df$LipidClass <- NULL
  expect_error(
    parse_lipid_chains(df, class_col = "LipidClass"),
    "not found"
  )
})
