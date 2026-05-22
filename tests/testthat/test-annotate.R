# tests/testthat/test-annotate.R

test_that("annotate_lipids() returns a data.frame with 5 annotation columns", {
  out <- annotate_lipids(MINIMAL_DATA, verbose = FALSE)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("LipidClass", "LipidClass_Full",
                    "LipidCategory_LMAPS", "LipidCategory_functional",
                    "LipidCategory") %in% names(out)))
  expect_equal(nrow(out), nrow(MINIMAL_DATA))
})

test_that("annotate_lipids() preserves all original columns", {
  out <- annotate_lipids(MINIMAL_DATA, verbose = FALSE)
  expect_true(all(names(MINIMAL_DATA) %in% names(out)))
})

test_that("annotate_lipids() errors on non-data.frame input", {
  expect_error(annotate_lipids(list(a = 1)), "must be a data.frame")
})

test_that("annotate_lipids() errors when lipid_col not found", {
  expect_error(
    annotate_lipids(MINIMAL_DATA, lipid_col = "nonexistent"),
    "not found"
  )
})

test_that(".assign_lipid_class() correctly classifies glycerophospholipids", {
  cases <- list(
    "PC 36:2"    = "PC",
    "PE 34:1"    = "PE",
    "PI 38:4"    = "PI",
    "PS 36:1"    = "PS",
    "LPC 18:0"   = "LPC",
    "LPE 16:0"   = "LPE"
  )
  for (name in names(cases)) {
    result <- easyLSEA:::.assign_lipid_class(name)
    expect_equal(result$LipidClass, cases[[name]],
                 label = paste("LipidClass for", name))
    expect_equal(result$LipidCategory_LMAPS, "Glycerophospholipids",
                 label = paste("Category for", name))
  }
})

test_that(".assign_lipid_class() correctly classifies sphingolipids", {
  cases <- list(
    "SM d18:1/16:0"   = "SM",
    "Cer(d18:1/24:0)" = "Cer",
    "HexCer 34:1"     = "HexCer",
    "Hex2Cer 36:1"    = "Hex2Cer"
  )
  for (name in names(cases)) {
    result <- easyLSEA:::.assign_lipid_class(name)
    expect_equal(result$LipidClass, cases[[name]],
                 label = paste("LipidClass for", name))
    expect_equal(result$LipidCategory_LMAPS, "Sphingolipids",
                 label = paste("Category for", name))
  }
})

test_that(".assign_lipid_class() correctly classifies glycerolipids", {
  cases <- list(
    "TG(54:3)"   = "TG",
    "DG 36:2"    = "DG",
    "TAG 52:2"   = "TG",   # alias
    "DAG 34:1"   = "DG"    # alias
  )
  for (name in names(cases)) {
    result <- easyLSEA:::.assign_lipid_class(name)
    expect_equal(result$LipidClass, cases[[name]],
                 label = paste("LipidClass for", name))
  }
})

test_that(".assign_lipid_class() correctly classifies ether lipids", {
  expect_equal(easyLSEA:::.assign_lipid_class("PC O-36:2")$LipidClass, "PC O")
  expect_equal(easyLSEA:::.assign_lipid_class("PE O-34:1")$LipidClass, "PE O")
  expect_equal(easyLSEA:::.assign_lipid_class("TG O-54:3")$LipidClass, "TG O")
})

test_that(".assign_lipid_class() returns Unknown for unrecognised names", {
  result <- easyLSEA:::.assign_lipid_class("Xanthophyll 42:0")
  expect_equal(result$LipidClass, "Unknown")
})

test_that(".assign_lipid_class() handles oxylipins and bile acids", {
  expect_equal(
    easyLSEA:::.assign_lipid_class("12-hydroxy-eicosatetraenoic acid")$LipidClass,
    "Oxylipin"
  )
  expect_equal(
    easyLSEA:::.assign_lipid_class("cholic acid")$LipidClass,
    "BA"
  )
})

test_that(".assign_lipid_class() handles acylcarnitines", {
  expect_equal(
    easyLSEA:::.assign_lipid_class("CAR 16:0")$LipidClass,
    "CAR"
  )
  expect_equal(
    easyLSEA:::.assign_lipid_class("palmitoylcarnitine")$LipidClass,
    "CAR"
  )
})

test_that("annotation summary prints unclassified count correctly", {
  df <- data.frame(
    LipidName = c("PC 36:2", "Unknown_lipid_XYZ"),
    logFC     = c(1, -1),
    stringsAsFactors = FALSE
  )
  expect_message(
    annotate_lipids(df, verbose = TRUE),
    "unclassified"
  )
})
