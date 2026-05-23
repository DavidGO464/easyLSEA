# R/utils.R
# Internal utility functions -- none exported.
# Sources:
#   .assign_lipid_class() / .get_lipid_category() adapted from
#   generate_lipids_full_v2_adapted.R (hierarchical prefix matching
#   validated against 394 VLDL features).

# Suppress R CMD check NOTEs for NSE column names used in dplyr/ggplot2
utils::globalVariables(c(
  # LSEA columns
  "Group", "N_group", "N_sig", "N_notsig",
  "N_sig_bg", "N_notsig_bg", "fisher_p",
  "Expected", "Fold_enr", "FDR_ORA",
  "Median_logFC", "N_up_case", "N_up_ref", "Enriched_in",
  "KS_pval", "FDR_LSEA", "DirectionalScore", "Direction",
  "NES", "FDR_fgsea", "LeadingEdge", "N_leading",
  # chain analysis columns
  "analysis_chain_cl", "analysis_chain_cs", "chain_type",
  "mean_logFC", "n_lip_cell", "n_sig_cell", "cl_f", "cs_f",
  # ggplot2 aes variables
  "n", "direction",
  # input data columns
  ".", "logFC", "AveExpr", "adj.P.Val", "sig",
  "LipidClass", "LipidCategory_LMAPS", "LipidCategory_functional",
  "LipidName", "Shorthand", "Confidence_rank",
  "pi_value", "t_stat", "weight"
))

# -- Entry-list constructors (reduce repetition) -------------------------------

#' @noRd
.gpl_entry <- function(cls, full) {
  list(LipidClass               = cls,
       LipidClass_Full          = full,
       LipidCategory_LMAPS      = "Glycerophospholipids",
       LipidCategory_functional = "Glycerophospholipids")
}

#' @noRd
.spl_entry <- function(cls, full) {
  list(LipidClass               = cls,
       LipidClass_Full          = full,
       LipidCategory_LMAPS      = "Sphingolipids",
       LipidCategory_functional = "Sphingolipids")
}

#' @noRd
.gl_entry <- function(cls, full) {
  list(LipidClass               = cls,
       LipidClass_Full          = full,
       LipidCategory_LMAPS      = "Glycerolipids",
       LipidCategory_functional = "Glycerolipids")
}

#' @noRd
.fa_entry <- function(cls, full) {
  list(LipidClass               = cls,
       LipidClass_Full          = full,
       LipidCategory_LMAPS      = "Fatty Acyls",
       LipidCategory_functional = "Fatty Acyls")
}

# -- Main annotation function --------------------------------------------------

#' Assign lipid class from a lipid name string
#'
#' Hierarchical prefix-matching against LIPID MAPS nomenclature.
#' Priority order matters: more specific prefixes are checked first.
#' Validated against 394 VLDL features covering all classes in a
#' real-world untargeted lipidomics dataset.
#'
#' @param name Character(1). Lipid name in shorthand or common notation.
#' @return Named list with four elements: \code{LipidClass},
#'   \code{LipidClass_Full}, \code{LipidCategory_LMAPS},
#'   \code{LipidCategory_functional}.
#'
#' @noRd
.assign_lipid_class <- function(name) {

  n <- trimws(name)

  # -- Bile acids (common names checked first -- no standard prefix) -----------
  if (grepl(
    paste0("cholic acid|chenodeoxycholic|ursodeoxycholic|lithocholic|",
           "taurocheno|glycochol|glycocheno|glycoursodeo|taurochenodeo"),
    n, ignore.case = TRUE, perl = TRUE
  )) {
    return(list(LipidClass               = "BA",
                LipidClass_Full          = "Bile acid",
                LipidCategory_LMAPS      = "Sterol Lipids",
                LipidCategory_functional = "Bile Acids"))
  }

  # -- Sterol (Cholesterol sulfate and related) -------------------------------
  if (grepl("^Cholesterol", n, ignore.case = TRUE)) {
    return(list(LipidClass               = "ST",
                LipidClass_Full          = "Sterol",
                LipidCategory_LMAPS      = "Sterol Lipids",
                LipidCategory_functional = "Sterol Lipids"))
  }

  # -- Oxylipins (hydroxy/oxo/keto FAs and eicosanoids) ----------------------
  # intToUtf8(177L) produces the +/- sign (U+00B1) at runtime without
  # embedding a non-ASCII character in the source file.
  oxylipin_pat <- paste0(
    "^([(][", intToUtf8(177L), "][)]",
    "|[0-9]+[(][A-Z][)]-|[0-9]+-hydroxy|[0-9]+-oxo|",
    "[0-9]+-keto|alpha-hydroxy|[0-9]+R-hydroxy|[0-9]+,|2E,4E|",
    "3-hydroxy|2-hydroxy|Hexacosanedioic|Prostaglandin|HpODE|",
    "isoprostane|resolvin|leukotriene|thromboxane)"
  )
  if (grepl(oxylipin_pat, n, ignore.case = TRUE, perl = TRUE)) {
    return(list(LipidClass               = "Oxylipin",
                LipidClass_Full          = "Oxylipin / Hydroxy-FA",
                LipidCategory_LMAPS      = "Fatty Acyls",
                LipidCategory_functional = "Oxylipins"))
  }

  # -- Ether phospholipids (MUST precede plain class checks) -----------------
  # PS(O-...) parenthesis format = plain PS (ether notation within parens)
  if (grepl("^PS[(]O-", n))         return(.gpl_entry("PS",   "PS"))
  if (grepl("^PS O-",   n))         return(.gpl_entry("PS O", "Ether-PS"))
  if (grepl("^PC O-|^PC P-", n))    return(.gpl_entry("PC O", "Ether-PC"))
  if (grepl("^PE O-|^PE P-", n))    return(.gpl_entry("PE O", "Ether-PE"))
  if (grepl("^PG O-", n))           return(.gpl_entry("PG O", "Ether-PG"))
  if (grepl("^MG O-", n))           return(.gl_entry("MG O",  "Ether-MG"))

  # -- PE-Ceramide (MUST precede plain PE) -----------------------------------
  if (grepl("^PE-Cer ", n))
    return(list(LipidClass               = "PE-Cer",
                LipidClass_Full          = "PE-Ceramide",
                LipidCategory_LMAPS      = "Glycerophospholipids",
                LipidCategory_functional = "Glycerophospholipids"))

  # -- Glycerophospholipids ---------------------------------------------------
  gp_map <- list(
    "^PC[ (]|^PC$"   = list("PC",    "PC"),
    "^PE[ (]|^PE$"   = list("PE",    "PE"),
    "^PI[ (]|^PI$"   = list("PI",    "PI"),
    "^PS[ (]|^PS[(]" = list("PS",    "PS"),
    "^PG[ (]|^PG$"   = list("PG",    "PG"),
    "^PA[ (]|^PA$"   = list("PA",    "PA"),
    "^LPC[ (]"       = list("LPC",   "LPC"),
    "^LPE[ (]"       = list("LPE",   "LPE"),
    "^LPI[ (]"       = list("LPI",   "LPI"),
    "^LPG[ (]"       = list("LPG",   "LPG"),
    "^LPA[ (]"       = list("LPA",   "LPA"),
    "^LPS[ (]"       = list("LPS",   "LPS"),
    "^LNAPE "        = list("LNAPE", "LNAPE")
  )
  for (pat in names(gp_map)) {
    if (grepl(pat, n)) {
      v <- gp_map[[pat]]
      return(.gpl_entry(v[[1L]], v[[2L]]))
    }
  }

  # -- Sphingolipids ----------------------------------------------------------
  if (grepl("^SHexCer ",           n)) return(.spl_entry("SHexCer", "Sulfatide"))
  if (grepl("^HexCer ",            n)) return(.spl_entry("HexCer",  "HexCer"))
  if (grepl("^GlcCer[(]|^GlcCer ", n)) return(.spl_entry("HexCer", "GlcCer"))
  if (grepl("^Hex2Cer ",           n)) return(.spl_entry("Hex2Cer", "Dihexosylceramide"))
  if (grepl("^Hex3Cer ",           n)) return(.spl_entry("Hex3Cer", "Trihexosylceramide"))
  if (grepl("^CerPE ",             n)) return(.spl_entry("CerPE",   "Ceramide-PE"))
  if (grepl("^Cer[ (]",            n)) return(.spl_entry("Cer",     "Ceramide"))
  if (grepl("^SM[ (]|^SM$",        n)) return(.spl_entry("SM",      "SM"))

  # -- Glycolipids (Saccharolipids per LIPID MAPS) ---------------------------
  if (grepl("^DGDG ", n))
    return(list(LipidClass               = "DGDG",
                LipidClass_Full          = "DGDG",
                LipidCategory_LMAPS      = "Saccharolipids",
                LipidCategory_functional = "Glycolipids"))
  if (grepl("^MGDG ", n))
    return(list(LipidClass               = "MGDG",
                LipidClass_Full          = "MGDG",
                LipidCategory_LMAPS      = "Saccharolipids",
                LipidCategory_functional = "Glycolipids"))

  # -- Ether glycerolipids (before plain TG/DG) ------------------------------
  if (grepl("^TG O-",       n))     return(.gl_entry("TG O", "Ether-TG"))
  if (grepl("^DG O-|^DG[(]P-", n)) return(.gl_entry("DG O", "Ether-DG"))

  # -- Glycerolipids ----------------------------------------------------------
  # DAG/TAG/MAG are aliases used by some lipidomics tools
  if (grepl("^MG[ (]|^MG$|^MAG[ (]|^MAG [0-9]", n)) return(.gl_entry("MG", "MG"))
  if (grepl("^DG[ (]|^DG$|^DAG[ (]|^DAG [0-9]", n)) return(.gl_entry("DG", "DG"))
  if (grepl("^TG[ (]|^TG$|^TAG[ (]|^TAG [0-9]", n)) return(.gl_entry("TG", "TG"))

  # -- Fatty Acyls -----------------------------------------------------------
  if (grepl("^FAHFA ",                             n)) return(.fa_entry("FAHFA", "FAHFA"))
  if (grepl("^FFA[ (]|^FFA[0-9]|[(]FFA",          n)) return(.fa_entry("FFA",   "FFA"))
  if (grepl("^FA[ (]|^FA [0-9]",                  n)) return(.fa_entry("FA",    "FA"))
  # Common fatty acid names (palmitic, oleic, etc.)
  named_fa_pat <- paste0(
    "palmitic acid|palmitoleic acid|stearic acid|oleic acid|",
    "linoleic acid|linolenic acid|myristic acid|lauric acid|",
    "lignoceric acid|arachidic acid|nervonic acid|behenic acid"
  )
  if (grepl(named_fa_pat, n, ignore.case = TRUE, perl = TRUE))
    return(.fa_entry("FFA", "FFA"))

  # -- Acylcarnitines --------------------------------------------------------
  if (grepl("^CAR [0-9]|^CAR$|carnitine", n, ignore.case = TRUE))
    return(list(LipidClass               = "CAR",
                LipidClass_Full          = "Acylcarnitine",
                LipidCategory_LMAPS      = "Fatty Acyls",
                LipidCategory_functional = "Acylcarnitines"))

  # -- Free and esterified cholesterol ---------------------------------------
  if (grepl("^FC ", n))
    return(list(LipidClass               = "FC",
                LipidClass_Full          = "Free Cholesterol",
                LipidCategory_LMAPS      = "Sterol Lipids",
                LipidCategory_functional = "Sterol Lipids"))
  if (grepl("^CE ", n))
    return(list(LipidClass               = "CE",
                LipidClass_Full          = "Cholesterol Ester",
                LipidCategory_LMAPS      = "Sterol Lipids",
                LipidCategory_functional = "Sterol Lipids"))

  # -- Coenzyme Q (Prenol Lipids) --------------------------------------------
  if (grepl("^CoQ", n))
    return(list(LipidClass               = "CoQ",
                LipidClass_Full          = "Coenzyme Q",
                LipidCategory_LMAPS      = "Prenol Lipids",
                LipidCategory_functional = "Prenol Lipids"))

  # -- Fallback ---------------------------------------------------------------
  list(LipidClass               = "Unknown",
       LipidClass_Full          = "Unclassified",
       LipidCategory_LMAPS      = "Unclassified",
       LipidCategory_functional = "Unclassified")
}

#' Derive functional category from LIPID MAPS category and lipid class
#'
#' Used when annotating with \code{annotate_lipid()} (full parser) to map
#' LIPID MAPS structural categories to the functional categories used by
#' easyLSEA. Oxylipins and Bile Acids are separated from Fatty Acyls;
#' Acylcarnitines are separated from Fatty Acyls; Saccharolipids become
#' Glycolipids.
#'
#' @param lm_category_name Character(1). LIPID MAPS category name from
#'   \code{annotate_lipid(detail = "full")}.
#' @param cls Character(1). Lipid class abbreviation (e.g. "HETE", "BA",
#'   "CAR").
#' @return Character(1). Functional category label.
#'
#' @noRd
.get_lipid_category_functional <- function(lm_category_name, cls) {
  oxylipin_classes <- c(
    "HETE", "HEPE", "EET", "DHET", "DiHETE", "OxoETE", "HpETE", "HHTrE",
    "EpETE", "HODE", "HOTrE", "EpOME", "HDoHE", "Resolvin", "Maresin",
    "Protectin", "LXA4", "LXB4", "LTB4", "LTC4", "LTD4",
    "PGE2", "PGD2", "PGF2a", "PGI2", "PGB2", "15d-PGJ2", "TXB2",
    "Oxylipin"
  )
  if (!is.na(cls) && cls %in% oxylipin_classes) return("Oxylipins")
  if (!is.na(cls) && cls == "BA")               return("Bile Acids")
  if (!is.na(cls) && cls == "CAR")              return("Acylcarnitines")
  if (!is.na(lm_category_name) &&
      lm_category_name == "Saccharolipids")     return("Glycolipids")
  if (is.na(lm_category_name))                  return("Unclassified")
  lm_category_name
}

#' Derive LipidCategory from LMAPS and functional category
#'
#' Oxylipins and Bile Acids are treated as standalone categories rather than
#' being nested under Fatty Acyls (which is the LIPID MAPS placement).
#' Saccharolipids are mapped to the more familiar term Glycolipids.
#'
#' @noRd
.get_lipid_category <- function(lmaps_cat, func_cat) {
  if (func_cat %in% c("Oxylipins", "Bile Acids")) return(func_cat)
  if (lmaps_cat == "Saccharolipids")              return("Glycolipids")
  lmaps_cat
}

#' Annotate a character vector of lipid names
#'
#' Applies \code{.assign_lipid_class()} and \code{.get_lipid_category()} to
#' each element and returns a \code{data.frame} with annotation columns
#' appended. Used internally by \code{annotate_lipids()} when
#' \code{method = "internal"}.
#'
#' @param names Character vector of lipid names.
#' @return \code{data.frame} with columns \code{LipidClass},
#'   \code{LipidClass_Full}, \code{LipidCategory_LMAPS},
#'   \code{LipidCategory_functional}, and \code{LipidCategory}.
#'
#' @noRd
.annotate_internal <- function(names) {
  annot <- lapply(names, .assign_lipid_class)

  out <- data.frame(
    LipidClass               = vapply(annot, `[[`, character(1L), "LipidClass"),
    LipidClass_Full          = vapply(annot, `[[`, character(1L), "LipidClass_Full"),
    LipidCategory_LMAPS      = vapply(annot, `[[`, character(1L), "LipidCategory_LMAPS"),
    LipidCategory_functional = vapply(annot, `[[`, character(1L), "LipidCategory_functional"),
    stringsAsFactors         = FALSE
  )

  out$LipidCategory <- mapply(
    .get_lipid_category,
    out$LipidCategory_LMAPS,
    out$LipidCategory_functional,
    SIMPLIFY = TRUE, USE.NAMES = FALSE
  )

  out
}

#' Validate the main input data.frame
#'
#' Called at the top of \code{easyLSEA()} and modular functions.
#' Stops with an informative message on hard failures; warns on soft ones.
#'
#' @noRd
.validate_input <- function(data, lipid_col, rank_col) {

  if (!is.data.frame(data))
    stop("'data' must be a data.frame, not ",
         class(data)[[1L]], ".", call. = FALSE)

  missing_cols <- setdiff(c(lipid_col, rank_col), names(data))
  if (length(missing_cols) > 0L)
    stop("Column(s) not found in 'data': ",
         paste(sQuote(missing_cols), collapse = ", "), ".\n",
         "Available columns: ",
         paste(sQuote(names(data)), collapse = ", "),
         call. = FALSE)

  if (!is.numeric(data[[rank_col]]))
    stop("Column ", sQuote(rank_col), " must be numeric, not ",
         class(data[[rank_col]])[[1L]], ".", call. = FALSE)

  n_na <- sum(is.na(data[[rank_col]]))
  if (n_na > 0L)
    warning(n_na, " NA value(s) in '", rank_col,
            "' will be removed before ranking.", call. = FALSE)

  n_dup <- sum(duplicated(data[[lipid_col]]))
  if (n_dup > 0L)
    warning(n_dup, " duplicate lipid name(s) in '", lipid_col,
            "'. Only the first occurrence will be used.", call. = FALSE)

  invisible(NULL)
}
