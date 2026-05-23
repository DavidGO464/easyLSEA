# =============================================================================
# annotate_lipid.R
#
# Universal lipid name parser with canonical LIPID MAPS classification
# and optional shorthand notation output.
#
# References:
#   Liebisch G et al. J Lipid Res. 2020;61(12):1539-1555. PMID: 33037133
#   Conroy MJ et al. Nucleic Acids Res. 2024;52(D1):D1677-D1682. PMID: 37855672
#   Fahy E et al. J Lipid Res. 2009;50 Suppl:S9-14. PMID: 19098281
#
# Usage:
#   df <- annotate_lipid(lipid_names)                          # compact (default)
#   df <- annotate_lipid(lipid_names, detail = "standard")
#   df <- annotate_lipid(lipid_names, detail = "full")
#   df <- annotate_lipid(lipid_names, shorthand = TRUE)
#   df <- annotate_lipid(lipid_names, detail = "full", shorthand = TRUE)
# =============================================================================


# =============================================================================
# OUTPUT COLUMNS BY DETAIL LEVEL
# =============================================================================

.LM_COLS <- list(

  compact = c(
    "Molecule", "Class", "is_lyso", "is_ambig", "not_matched",
    "chain1", "chain2", "chain3", "chain4",
    "l_1", "s_1", "l_2", "s_2", "l_3", "s_3", "l_4", "s_4",
    "total_cl", "total_cs", "n_chains"
  ),

  standard = c(
    "Molecule", "Class", "lm_category", "lm_class_id", "annotation_level",
    "is_lyso", "is_ambig", "is_ether", "is_plasmalogen", "is_istd",
    "not_matched", "sphingoid_prefix",
    "chain1", "chain2", "chain3", "chain4",
    "l_1", "s_1", "l_2", "s_2", "l_3", "s_3", "l_4", "s_4",
    "total_cl", "total_cs", "n_chains"
  ),

  full = c(
    "Molecule", "clean_name", "Class", "class_stub",
    "lm_category", "lm_category_name", "lm_class_id", "lm_class_name",
    "headgroup", "annotation_level",
    "is_lyso", "is_ambig", "is_ether", "is_plasmalogen", "is_istd",
    "not_matched", "ether_link", "sphingoid_prefix", "hydroxyl_base",
    "chain1", "chain2", "chain3", "chain4",
    "l_1", "s_1", "l_2", "s_2", "l_3", "s_3", "l_4", "s_4",
    "total_cl", "total_cs", "n_chains"
  )
)


# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Annotate lipid names with LIPID MAPS classification
#'
#' Parses lipid names in any format used by lipidomics software (LipidSearch,
#' MS-DIAL, LipidView) and returns a structured data frame with LIPID MAPS
#' canonical classification, chain-level metadata, and optional shorthand
#' notation per Liebisch et al. (2020).
#'
#' @param molecules Character vector of lipid names to parse.
#' @param detail Level of detail in the output table:
#'   \describe{
#'     \item{\code{"compact"}}{Class, chains, totals. Equivalent to lipidR. Default.}
#'     \item{\code{"standard"}}{Adds LIPID MAPS category, class ID, and structural flags.}
#'     \item{\code{"full"}}{All columns. Recommended for joining with LMSD.}
#'   }
#' @param shorthand Logical. If \code{TRUE}, adds \code{shorthand_lm} column
#'   with canonical shorthand per Liebisch et al. (2020). Default \code{FALSE}.
#' @param sn_confirmed Logical. If \code{TRUE}, marks chains as sn-confirmed
#'   in shorthand output. Requires MS/MS-directed analysis. Default \code{FALSE}.
#' @param lyso_explicit Logical. If \code{TRUE}, lyso-lipids include the empty
#'   sn-position: \code{LPC(18:1/0:0)} instead of \code{LPC(18:1)}.
#'   Default \code{FALSE}.
#' @param no_match How to handle unparsed names: \code{"warn"} (default),
#'   \code{"remove"}, or \code{"ignore"}.
#' @param sphingoid_default Default sphingoid base prefix for sphingolipids
#'   without explicit prefix. \code{"d"} = dihydroxy (mammalian default).
#'   Use \code{NA} for non-mammalian data.
#'
#' @return A data frame with one row per unique lipid name. Key columns include
#'   \code{Class}, \code{lm_category}, \code{lm_class_id}, \code{annotation_level},
#'   \code{is_ether}, \code{is_plasmalogen}, \code{is_istd}, \code{sphingoid_prefix},
#'   \code{total_cl}, \code{total_cs}, and optionally \code{shorthand_lm}.
#'
#' @references
#' Liebisch G et al. Update on LIPID MAPS classification, nomenclature, and
#' shorthand notation for MS-derived lipid structures.
#' \emph{J Lipid Res.} 2020;61(12):1539-1555. \doi{10.1194/jlr.S120001025}
#'
#' Conroy MJ et al. LIPID MAPS: update to databases and tools for the
#' lipidomics community. \emph{Nucleic Acids Res.} 2024;52(D1):D1677-D1682.
#' \doi{10.1093/nar/gkad896}
#'
#' @examples
#' lipids <- c("PC 16:0/18:1", "PC O-18:1/20:4", "Cer d18:1/16:0",
#'             "TG(16:0/18:1/18:1)", "Lyso PE 18:1(d7)",
#'             "plasmenylPE (16:0/18:1)", "Sa1P d 18:0",
#'             "WE 16:0/18:1", "CoA 16:0",
#'             "15-HETE", "PGE2", "LTB4", "Resolvin D1", "12(13)-EpOME")
#'
#' annotate_lipid(lipids)
#' annotate_lipid(lipids, detail = "standard")
#' annotate_lipid(lipids, detail = "full", shorthand = TRUE)
#'
#' @export
annotate_lipid <- function(molecules,
                           detail            = c("compact", "standard", "full"),
                           shorthand         = FALSE,
                           sn_confirmed      = FALSE,
                           lyso_explicit     = FALSE,
                           no_match          = c("warn", "remove", "ignore"),
                           sphingoid_default = "d") {

  detail   <- match.arg(detail)
  no_match <- match.arg(no_match)
  molecules <- unique(as.character(trimws(molecules)))

  results <- lapply(molecules, .parse_lm_single,
                    sphingoid_default = sphingoid_default)
  df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))

  num_cols <- c("l_1", "s_1", "l_2", "s_2", "l_3", "s_3", "l_4", "s_4",
                "total_cl", "total_cs", "n_chains", "hydroxyl_base")
  for (col in num_cols) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))

  if (any(df$not_matched)) {
    bad <- df$Molecule[df$not_matched]
    if (no_match == "warn") {
      warning("Lipid names not parsed:\n  ", paste(bad, collapse = ", "))
    } else if (no_match == "remove") {
      df <- df[!df$not_matched, ]
    }
  }

  if (shorthand) {
    df$shorthand_lm <- vapply(
      seq_len(nrow(df)),
      function(i) {
        row <- df[i, ]
        if (isTRUE(row$not_matched) || is.na(row$Class)) return(NA_character_)
        .build_shorthand(row,
                         sn_confirmed  = sn_confirmed,
                         lyso_explicit = lyso_explicit)
      },
      character(1)
    )
  }

  cols <- .LM_COLS[[detail]]
  if (shorthand) cols <- c(cols, "shorthand_lm")
  df[, cols, drop = FALSE]
}


# =============================================================================
# LIPID MAPS CLASS LOOKUP TABLE
# Liebisch et al. 2020 + Conroy et al. 2024
# =============================================================================

.LM_CLASS_TABLE <- local({

  cls <- function(abbrev, cat, cat_name, cls_id, cls_name,
                  n_chains, headgroup = NA_character_, sphingoid = FALSE) {
    data.frame(
      abbrev             = abbrev,
      lm_category        = cat,
      lm_category_name   = cat_name,
      lm_class_id        = cls_id,
      lm_class_name      = cls_name,
      n_chains_expected  = n_chains,
      headgroup          = headgroup,
      sphingoid_base_req = sphingoid,
      stringsAsFactors   = FALSE
    )
  }

  rbind(
    # ── Fatty Acyls ────────────────────────────────────────────────────────
    cls("FA",       "FA", "Fatty Acyls", "FA0101", "Fatty acids",                           1),
    cls("FAHFA",    "FA", "Fatty Acyls", "FA0701", "Fatty acid estolides",                  2),
    cls("FAL",      "FA", "Fatty Acyls", "FA0601", "Fatty aldehydes",                       1),
    cls("FOH",      "FA", "Fatty Acyls", "FA0501", "Fatty alcohols",                        1),
    cls("CAR",      "FA", "Fatty Acyls", "FA0701", "Acylcarnitines",                        1),
    cls("NAE",      "FA", "Fatty Acyls", "FA0801", "N-Acylethanolamines",                   1),
    # Wax esters & acyl-CoA
    cls("WE",       "FA", "Fatty Acyls", "FA0901", "Wax esters",                            2),
    cls("CoA",      "FA", "Fatty Acyls", "FA0702", "Acyl coenzyme As",                      1),
    # Eicosanoids — prostanoids
    cls("PGE2",     "FA", "Fatty Acyls", "FA0301", "Prostaglandins",                        0),
    cls("PGD2",     "FA", "Fatty Acyls", "FA0301", "Prostaglandins",                        0),
    cls("PGF2a",    "FA", "Fatty Acyls", "FA0301", "Prostaglandins",                        0),
    cls("PGI2",     "FA", "Fatty Acyls", "FA0301", "Prostaglandins",                        0),
    cls("PGB2",     "FA", "Fatty Acyls", "FA0301", "Prostaglandins",                        0),
    cls("15d-PGJ2", "FA", "Fatty Acyls", "FA0301", "Prostaglandins",                        0),
    cls("TXB1",     "FA", "Fatty Acyls", "FA0302", "Thromboxanes",                          0),
    cls("TXB2",     "FA", "Fatty Acyls", "FA0302", "Thromboxanes",                          0),
    cls("TXB3",     "FA", "Fatty Acyls", "FA0302", "Thromboxanes",                          0),
    # Eicosanoids — leukotrienes
    cls("LTB4",     "FA", "Fatty Acyls", "FA0303", "Leukotrienes",                          0),
    cls("LTC4",     "FA", "Fatty Acyls", "FA0303", "Leukotrienes",                          0),
    cls("LTD4",     "FA", "Fatty Acyls", "FA0303", "Leukotrienes",                          0),
    # Eicosanoids — HETEs / HEPEs / EETs / oxo
    cls("HETE",     "FA", "Fatty Acyls", "FA0304", "Hydroxyeicosatetraenoic acids",         0),
    cls("HEPE",     "FA", "Fatty Acyls", "FA0304", "Hydroxyeicosapentaenoic acids",         0),
    cls("EET",      "FA", "Fatty Acyls", "FA0304", "Epoxyeicosatrienoic acids",             0),
    cls("DHET",     "FA", "Fatty Acyls", "FA0304", "Dihydroxyeicosatrienoic acids",         0),
    cls("DiHETE",   "FA", "Fatty Acyls", "FA0304", "Dihydroxyeicosatetraenoic acids",       0),
    cls("OxoETE",   "FA", "Fatty Acyls", "FA0304", "Oxoeicosatetraenoic acids",             0),
    cls("HpETE",    "FA", "Fatty Acyls", "FA0304", "Hydroperoxyeicosatetraenoic acids",     0),
    cls("HHTrE",    "FA", "Fatty Acyls", "FA0304", "Hydroxyheptadecatrienoic acids",        0),
    cls("EpETE",    "FA", "Fatty Acyls", "FA0304", "Epoxyeicosatetraenoic acids",           0),
    # Octadecanoids
    cls("HODE",     "FA", "Fatty Acyls", "FA0305", "Hydroxyoctadecadienoic acids",          0),
    cls("HOTrE",    "FA", "Fatty Acyls", "FA0305", "Hydroxyoctadecatrienoic acids",         0),
    cls("EpOME",    "FA", "Fatty Acyls", "FA0305", "Epoxyoctadecenoic acids",               0),
    # Docosanoids / specialized pro-resolving mediators
    cls("HDoHE",    "FA", "Fatty Acyls", "FA0306", "Hydroxydocosahexaenoic acids",          0),
    cls("Resolvin", "FA", "Fatty Acyls", "FA0306", "Resolvins",                             0),
    cls("Maresin",  "FA", "Fatty Acyls", "FA0306", "Maresins",                              0),
    cls("Protectin","FA", "Fatty Acyls", "FA0306", "Protectins",                            0),
    cls("LXA4",     "FA", "Fatty Acyls", "FA0306", "Lipoxins",                              0),
    cls("LXB4",     "FA", "Fatty Acyls", "FA0306", "Lipoxins",                              0),

    # ── Glycerolipids ──────────────────────────────────────────────────────
    cls("MG",       "GL", "Glycerolipids", "GL0101", "Monoradylglycerols",                  1),
    cls("DG",       "GL", "Glycerolipids", "GL0201", "Diradylglycerols",                    2),
    cls("TG",       "GL", "Glycerolipids", "GL0301", "Triradylglycerols",                   3),
    cls("MGDG",     "GL", "Glycerolipids", "GL0501", "Monogalactosyldiacylglycerols",       2),
    cls("DGDG",     "GL", "Glycerolipids", "GL0502", "Digalactosyldiacylglycerols",         2),
    cls("SQDG",     "GL", "Glycerolipids", "GL0503", "Sulfoquinovosyldiacylglycerols",      2),
    cls("SQMG",     "GL", "Glycerolipids", "GL0504", "Sulfoquinovosylmonoacylglycerols",    1),

    # ── Glycerophospholipids ───────────────────────────────────────────────
    cls("PA",           "GP", "Glycerophospholipids", "GP1001", "Phosphatidic acids",                2, "phosphate"),
    cls("LPA",          "GP", "Glycerophospholipids", "GP1002", "Lysophosphatidic acids",            1, "phosphate"),
    cls("PC",           "GP", "Glycerophospholipids", "GP0101", "Phosphatidylcholines",              2, "choline"),
    cls("LPC",          "GP", "Glycerophospholipids", "GP0105", "Lysophosphatidylcholines",          1, "choline"),
    cls("PE",           "GP", "Glycerophospholipids", "GP0201", "Phosphatidylethanolamines",         2, "ethanolamine"),
    cls("LPE",          "GP", "Glycerophospholipids", "GP0205", "Lysophosphatidylethanolamines",     1, "ethanolamine"),
    cls("PG",           "GP", "Glycerophospholipids", "GP0401", "Phosphatidylglycerols",             2, "glycerol"),
    cls("LPG",          "GP", "Glycerophospholipids", "GP0405", "Lysophosphatidylglycerols",         1, "glycerol"),
    cls("PI",           "GP", "Glycerophospholipids", "GP0601", "Phosphatidylinositols",             2, "inositol"),
    cls("LPI",          "GP", "Glycerophospholipids", "GP0605", "Lysophosphatidylinositols",         1, "inositol"),
    cls("PS",           "GP", "Glycerophospholipids", "GP0301", "Phosphatidylserines",               2, "serine"),
    cls("LPS",          "GP", "Glycerophospholipids", "GP0305", "Lysophosphatidylserines",           1, "serine"),
    cls("PIP",          "GP", "Glycerophospholipids", "GP0701", "Phosphoinositolphosphates",         2, "inositolP"),
    cls("PIP2",         "GP", "Glycerophospholipids", "GP0702", "Bisphosphoinositolphosphates",      2, "inositolPP"),
    cls("PIP3",         "GP", "Glycerophospholipids", "GP0703", "Trisphosphoinositolphosphates",     2, "inositolPPP"),
    cls("CL",           "GP", "Glycerophospholipids", "GP1201", "Cardiolipins",                      4, "glycerophosphate"),
    cls("MLCL",         "GP", "Glycerophospholipids", "GP1203", "Monolysocardiolipins",              3, "glycerophosphate"),
    cls("DLCL",         "GP", "Glycerophospholipids", "GP1202", "Dilysocardiolipins",                2, "glycerophosphate"),
    cls("BMP",          "GP", "Glycerophospholipids", "GP0401", "Bis(monoacylglycero)phosphates",    2, "glycerol"),
    cls("LBPA",         "GP", "Glycerophospholipids", "GP0401", "Lysobiophosphatidic acids",         2, "glycerol"),
    cls("PEtOH",        "GP", "Glycerophospholipids", "GP0210", "Phosphatidylethanol",               2, "ethanol"),
    cls("PThr",         "GP", "Glycerophospholipids", "GP0901", "Phosphatidylthreonines",            2, "threonine"),
    # FIX 1: Trivial plasmalogen names — marked as plasmalogen in parser
    cls("plasmenylPC",  "GP", "Glycerophospholipids", "GP0101", "Phosphatidylcholines",              2, "choline"),
    cls("plasmenylPE",  "GP", "Glycerophospholipids", "GP0201", "Phosphatidylethanolamines",         2, "ethanolamine"),

    # ── Sphingolipids ──────────────────────────────────────────────────────
    cls("Cer",      "SP", "Sphingolipids", "SP0201", "Ceramides",                            2, NA_character_, TRUE),
    cls("CerP",     "SP", "Sphingolipids", "SP0204", "Ceramide-1-phosphates",                2, "phosphate",   TRUE),
    cls("SM",       "SP", "Sphingolipids", "SP0301", "Sphingomyelins",                       2, "phosphocholine", TRUE),
    cls("LSM",      "SP", "Sphingolipids", "SP0305", "Lysosphingomyelins",                   1, "phosphocholine", TRUE),
    cls("HexCer",   "SP", "Sphingolipids", "SP0501", "Hexosylceramides",                     2, "hexose",      TRUE),
    cls("LHexCer",  "SP", "Sphingolipids", "SP0505", "Lysomonohexosylceramides",             1, "hexose",      TRUE),
    cls("Hex2Cer",  "SP", "Sphingolipids", "SP0601", "Dihexosylceramides",                   2, "hexose2",     TRUE),
    cls("GlcCer",   "SP", "Sphingolipids", "SP0501", "Glucosylceramides",                    2, "glucose",     TRUE),
    cls("GalCer",   "SP", "Sphingolipids", "SP0501", "Galactosylceramides",                  2, "galactose",   TRUE),
    cls("LacCer",   "SP", "Sphingolipids", "SP0601", "Lactosylceramides",                    2, "lactose",     TRUE),
    cls("SHexCer",  "SP", "Sphingolipids", "SP0501", "Sulfohexosylceramides",                2, "sulfohexose", TRUE),
    cls("GM3",      "SP", "Sphingolipids", "SP0603", "Ganglioside GM3",                      2, "ganglioside", TRUE),
    cls("GM2",      "SP", "Sphingolipids", "SP0604", "Ganglioside GM2",                      2, "ganglioside", TRUE),
    cls("GM1",      "SP", "Sphingolipids", "SP0605", "Ganglioside GM1",                      2, "ganglioside", TRUE),
    cls("GD1",      "SP", "Sphingolipids", "SP0606", "Ganglioside GD1",                      2, "ganglioside", TRUE),
    cls("GD3",      "SP", "Sphingolipids", "SP0607", "Ganglioside GD3",                      2, "ganglioside", TRUE),
    cls("GT3",      "SP", "Sphingolipids", "SP0608", "Ganglioside GT3",                      2, "ganglioside", TRUE),
    cls("Sph",      "SP", "Sphingolipids", "SP0101", "Sphingoid bases",                      1, NA_character_, TRUE),
    cls("S1P",      "SP", "Sphingolipids", "SP0105", "Sphingoid base-1-phosphates",          1, "phosphate",   TRUE),
    # FIX 2: Trivial sphingoid base names from legacy software
    cls("Sa",       "SP", "Sphingolipids", "SP0101", "Sphingoid bases",                      1, NA_character_, TRUE),
    cls("So",       "SP", "Sphingolipids", "SP0101", "Sphingoid bases",                      1, NA_character_, TRUE),
    cls("Sa1P",     "SP", "Sphingolipids", "SP0105", "Sphingoid base-1-phosphates",          1, "phosphate",   TRUE),
    cls("So1P",     "SP", "Sphingolipids", "SP0105", "Sphingoid base-1-phosphates",          1, "phosphate",   TRUE),
    cls("PE-Cer",   "SP", "Sphingolipids", "SP0201", "PE-ceramides",                         2, "ethanolamine", TRUE),
    cls("PI-Cer",   "SP", "Sphingolipids", "SP0201", "PI-ceramides",                         2, "inositol",    TRUE),
    cls("PA-Cer",   "SP", "Sphingolipids", "SP0201", "PA-ceramides",                         2, "phosphate",   TRUE),
    cls("IPC",      "SP", "Sphingolipids", "SP0201", "Inositolphosphoceramides",             2, "inositol",    TRUE),
    cls("deoxyCer", "SP", "Sphingolipids", "SP0202", "1-Deoxy ceramides",                    2, NA_character_, TRUE),

    # ── Sterol Lipids ──────────────────────────────────────────────────────
    cls("Chol",     "ST", "Sterol Lipids", "ST0101", "Cholesterol",                          0),
    cls("CE",       "ST", "Sterol Lipids", "ST0102", "Cholesterol esters",                   1, "cholesterol"),
    cls("FC",       "ST", "Sterol Lipids", "ST0101", "Free cholesterol",                     0),
    cls("BA",       "ST", "Sterol Lipids", "ST0401", "Bile acids",                           0),

    # ── Prenol Lipids ──────────────────────────────────────────────────────
    cls("CoQ",      "PR", "Prenol Lipids", "PR0603", "Ubiquinones",                          0)
  )
})


# =============================================================================
# LIPID MEDIATOR TRIVIAL NAME MAP
# Handles positional prefixes: "15-HETE", "5(6)-EET", "12(13)-EpOME",
# "11,12-DHET", "Resolvin D1", "Maresin 1", "LTB4", "PGE2", "TXB2", etc.
# Checked BEFORE general parsing in .parse_lm_single.
# =============================================================================

.MEDIATOR_MAP <- list(
  # HETEs — e.g. "15-HETE", "5(S)-HETE", "12-HETE"
  list(pat = "^(\\d+[\\(\\-][^-]*-?)?HETE$",      cls = "HETE",     id = "FA0304", nm = "Hydroxyeicosatetraenoic acids"),
  # HEPEs — e.g. "12-HEPE", "15-HEPE"
  list(pat = "^(\\d+-)?HEPE$",                     cls = "HEPE",     id = "FA0304", nm = "Hydroxyeicosapentaenoic acids"),
  # EETs — e.g. "5(6)-EET", "11(12)-EET", "14(15)-EET"
  list(pat = "^\\d+\\(?\\d*\\)?-EET$",             cls = "EET",      id = "FA0304", nm = "Epoxyeicosatrienoic acids"),
  # DHETs — e.g. "11,12-DHET", "14,15-DHET", "8,9-DHET"
  list(pat = "^\\d+,\\d+-DHET$",                   cls = "DHET",     id = "FA0304", nm = "Dihydroxyeicosatrienoic acids"),
  # DiHETEs — e.g. "5,12-DiHETE", "5,6-DiHETE"
  list(pat = "^\\d+,\\d+-DiHETE$",                 cls = "DiHETE",   id = "FA0304", nm = "Dihydroxyeicosatetraenoic acids"),
  # OxoETEs — e.g. "5-OxoETE", "12-OxoETE"
  list(pat = "^\\d+-OxoETE$",                      cls = "OxoETE",   id = "FA0304", nm = "Oxoeicosatetraenoic acids"),
  # HpETE — e.g. "5-HpETE"
  list(pat = "^\\d+-HpETE$",                       cls = "HpETE",    id = "FA0304", nm = "Hydroperoxyeicosatetraenoic acids"),
  # HHTrE — e.g. "12-HHTrE"
  list(pat = "^\\d+-HHTrE$",                       cls = "HHTrE",    id = "FA0304", nm = "Hydroxyheptadecatrienoic acids"),
  # EpETEs — e.g. "14(15)-EpETE"
  list(pat = "^\\d+\\(?\\d*\\)?-EpETE$",           cls = "EpETE",    id = "FA0304", nm = "Epoxyeicosatetraenoic acids"),
  # HODEs — e.g. "9-HODE", "13-HODE"
  list(pat = "^\\d+-HODE$",                        cls = "HODE",     id = "FA0305", nm = "Hydroxyoctadecadienoic acids"),
  # HOTrEs — e.g. "9-HOTrE", "13-HOTrE"
  list(pat = "^\\d+-HOTrE$",                       cls = "HOTrE",    id = "FA0305", nm = "Hydroxyoctadecatrienoic acids"),
  # EpOMEs — e.g. "9(10)-EpOME", "12(13)-EpOME"
  list(pat = "^\\d+\\(?\\d*\\)?-EpOME$",           cls = "EpOME",    id = "FA0305", nm = "Epoxyoctadecenoic acids"),
  # HDoHEs — e.g. "10-HDoHE", "11-HDoHE"
  list(pat = "^\\d+-HDoHE$",                       cls = "HDoHE",    id = "FA0306", nm = "Hydroxydocosahexaenoic acids"),
  # Resolvins — e.g. "Resolvin D1", "Resolvin D2", "Resolvin E1"
  list(pat = "^Resolvin\\s+[DEA-Z]\\d+$",          cls = "Resolvin", id = "FA0306", nm = "Resolvins"),
  # Maresins — e.g. "Maresin 1", "Maresin 2"
  list(pat = "^Maresin\\s+\\d+$",                  cls = "Maresin",  id = "FA0306", nm = "Maresins"),
  # Protectins — e.g. "Protectin D1"
  list(pat = "^Protectin\\s+D\\d+$",               cls = "Protectin",id = "FA0306", nm = "Protectins"),
  # Lipoxins — e.g. "LXA4", "LXB4", "5,6,15-LXA4"
  list(pat = "^(\\d+,\\d+,\\d+-)?LX[AB]\\d+$",    cls = "LXA4",     id = "FA0306", nm = "Lipoxins"),
  # Leukotrienes — e.g. "LTB4", "LTC4", "LTD4"
  list(pat = "^LTB\\d+$",                          cls = "LTB4",     id = "FA0303", nm = "Leukotrienes"),
  list(pat = "^LT[CD]\\d+$",                       cls = "LTC4",     id = "FA0303", nm = "Leukotrienes"),
  # Prostaglandins — e.g. "PGE2", "PGD2", "PGF2alpha", "PGB2"
  list(pat = "^PGE\\d+$",                          cls = "PGE2",     id = "FA0301", nm = "Prostaglandins"),
  list(pat = "^PGD\\d+$",                          cls = "PGD2",     id = "FA0301", nm = "Prostaglandins"),
  list(pat = "^PGF\\d+[a-z]*$",                    cls = "PGF2a",    id = "FA0301", nm = "Prostaglandins"),
  list(pat = "^PGI\\d+$",                          cls = "PGI2",     id = "FA0301", nm = "Prostaglandins"),
  list(pat = "^PGB\\d+$",                          cls = "PGB2",     id = "FA0301", nm = "Prostaglandins"),
  list(pat = "^15d-PGJ\\d+$",                      cls = "15d-PGJ2", id = "FA0301", nm = "Prostaglandins"),
  # Thromboxanes — e.g. "TXB1", "TXB2", "TXB3"
  list(pat = "^TXB\\d+$",                          cls = "TXB2",     id = "FA0302", nm = "Thromboxanes"),
  # Tetranor metabolites — e.g. "Tetranor-12-HETE"
  list(pat = "^[Tt]etranor-\\d+-HETE$",            cls = "HETE",     id = "FA0304", nm = "Hydroxyeicosatetraenoic acids")
)


# =============================================================================
# INTERNAL PARSER — single lipid
# =============================================================================

.parse_lm_single <- function(mol, sphingoid_default = "d") {

  original <- mol
  x        <- trimws(mol)

  out <- list(
    Molecule = original, clean_name = NA_character_, not_matched = FALSE,
    annotation_level = NA_character_, is_ambig = FALSE, is_lyso = FALSE,
    is_istd = FALSE, is_ether = FALSE, is_plasmalogen = FALSE,
    lm_category = NA_character_, lm_category_name = NA_character_,
    lm_class_id = NA_character_, lm_class_name = NA_character_,
    Class = NA_character_, class_stub = NA_character_, headgroup = NA_character_,
    ether_link = NA_character_, sphingoid_prefix = NA_character_,
    hydroxyl_base = NA_real_,
    chain1 = NA_character_, chain2 = NA_character_,
    chain3 = NA_character_, chain4 = NA_character_,
    l_1 = NA_real_, s_1 = NA_real_, l_2 = NA_real_, s_2 = NA_real_,
    l_3 = NA_real_, s_3 = NA_real_, l_4 = NA_real_, s_4 = NA_real_,
    total_cl = NA_real_, total_cs = NA_real_, n_chains = NA_real_
  )

  # ISTD detection: (d7) or -d7 = deuterium label; d18:1 is NOT an ISTD (sphingoid prefix)
  if (grepl("\\(d\\d+\\)|-d\\d+(?!:)|\\(\\d+,\\d+-d\\d+\\)", x, perl = TRUE)) {
    out$is_istd <- TRUE
  }

  x <- gsub("[\\[\\]]", "", x, perl = TRUE)
  x <- sub("\\s+NEG$", "", x)
  x <- sub("\\s+ID\\d+$", "", x)
  x <- sub("15-MHDA", "17:0", x)

  # Special case: 15-MHDA trivial name converts to FA 17:0
  if (trimws(x) == "17:0") {
    out$clean_name <- "17:0"; out$Class <- "FA"; out$class_stub <- "FA"
    out$lm_category <- "FA"; out$lm_category_name <- "Fatty Acyls"
    out$lm_class_id <- "FA0101"; out$lm_class_name <- "Fatty acids"
    out$annotation_level <- "species"; out$is_ambig <- TRUE
    out$chain1 <- "17:0"; out$l_1 <- 17; out$s_1 <- 0
    out$total_cl <- 17; out$total_cs <- 0; out$n_chains <- 1
    return(out)
  }

  # ── Lipid mediator trivial name check ─────────────────────────────────────
  # Must run BEFORE normalization to preserve positional notation (e.g. "15-HETE")
  for (med in .MEDIATOR_MAP) {
    if (grepl(med$pat, x, perl = TRUE)) {
      out$Class            <- med$cls
      out$class_stub       <- med$cls
      out$lm_category      <- "FA"
      out$lm_category_name <- "Fatty Acyls"
      out$lm_class_id      <- med$id
      out$lm_class_name    <- med$nm
      out$annotation_level <- "species"
      out$n_chains         <- 0
      out$total_cl         <- 0
      out$total_cs         <- 0
      out$clean_name       <- x
      return(out)
    }
  }

  # Normalize: "Lyso PE" -> "LPE", number-first -> class-first,
  # C-prefix on chains, spaces between chains -> slash
  x <- gsub("(?i)lyso\\s*([A-Z])", "L\\1", x, perl = TRUE)
  # Fix: "Sa1P d 18:0" -> "Sa1P d18:0" (collapse space between sphingoid prefix and chain)
  x <- gsub("\\s+([dmt])\\s+(\\d+:\\d+)", " \\1\\2", x, perl = TRUE)
  # C-prefix on chains, spaces between chains -> slash
  x <- gsub("(?i)lyso\\s*([A-Z])", "L\\1", x, perl = TRUE)
  if (grepl("^\\d+:\\d+", x))
    x <- sub("^([^[:blank:]]+)[[:blank:]]+(.+)$", "\\2 \\1", x)
  x <- gsub("(?<=[(/\\s,]|^)C(\\d+:\\d+)", "\\1", x, perl = TRUE)
  chain_sp <- "(\\d+:\\d+[^[:space:]/]*)[[:space:]]+((?:[OP]-)?[dmt]?\\d+:\\d+)"
  x <- gsub(chain_sp, "\\1/\\2", x)
  x <- gsub(chain_sp, "\\1/\\2", x)
  x <- gsub("(\\d+:\\d+)-(?=[dmt]?\\d+:\\d+)", "\\1/", x, perl = TRUE)
  out$clean_name <- x

  # Split class and chains
  m_p <- regmatches(x, regexec(
    "^([A-Za-z][A-Za-z0-9-]*)\\s*\\((.+)\\)\\s*$", x, perl = TRUE))[[1]]
  if (length(m_p) == 3) {
    class_raw <- m_p[2]; chains_raw <- m_p[3]
  } else {
    m_s <- regmatches(x, regexec(
      "^([A-Za-z][A-Za-z0-9-]*)\\s+(.+)$", x, perl = TRUE))[[1]]
    if (length(m_s) == 3) {
      class_raw <- m_s[2]; chains_raw <- m_s[3]
    } else {
      class_raw <- x; chains_raw <- ""
    }
  }
  class_raw <- trimws(class_raw)

  # LIPID MAPS classification lookup
  tbl <- .LM_CLASS_TABLE
  idx <- match(class_raw, tbl$abbrev)
  if (!is.na(idx)) {
    out$lm_category <- tbl$lm_category[idx]; out$lm_category_name <- tbl$lm_category_name[idx]
    out$lm_class_id <- tbl$lm_class_id[idx]; out$lm_class_name <- tbl$lm_class_name[idx]
    out$headgroup <- tbl$headgroup[idx]; out$Class <- class_raw
  } else {
    out$lm_category <- "Unknown"; out$lm_category_name <- "Unknown"; out$Class <- class_raw
  }
  out$class_stub <- class_raw
  out$is_lyso    <- grepl("^L[A-Z]", class_raw) && !class_raw %in% c("LacCer", "LHexCer", "LSM", "LXA4", "LXB4")

  # FIX 3: Trivial plasmalogen names -> mark as plasmalogen (P- ether by definition)
  if (class_raw %in% c("plasmenylPC", "plasmenylPE")) {
    out$is_ether       <- TRUE
    out$is_plasmalogen <- TRUE
    out$ether_link     <- "P"
  }

  # Parse chains
  chains_str <- trimws(chains_raw)
  if (nchar(chains_str) == 0) {
    out$annotation_level <- "species"; out$n_chains <- 0
    out$total_cl <- 0; out$total_cs <- 0; return(out)
  }

  # Ether linkage prefix: O- (alkyl) or P- (alkenyl/plasmalogen)
  eth_m <- regmatches(chains_str, regexpr("^([OP])-", chains_str, perl = TRUE))
  if (length(eth_m) > 0 && nchar(eth_m) > 0) {
    out$ether_link     <- substr(eth_m, 1, 1)
    out$is_ether       <- TRUE
    out$is_plasmalogen <- substr(eth_m, 1, 1) == "P"
    chains_str         <- sub("^[OP]-", "", chains_str)
  }

  chain_parts <- trimws(strsplit(chains_str, "/")[[1]])
  chain_parts <- chain_parts[seq_len(min(length(chain_parts), 4))]

  # FIX 4: Remove empty chain positions (0:0) — present in DG/TG from some software
  # e.g. "DG (14:0/16:0/0:0)" -> only count 14:0 and 16:0
  chain_parts <- chain_parts[chain_parts != "0:0"]
  n <- length(chain_parts)

  if (n == 0) {
    out$not_matched <- TRUE
    return(out)
  }

  # Annotation level: species (sum composition) vs molecular_species (chain-resolved)
  out$is_ambig         <- (n == 1 && !out$is_lyso)
  out$annotation_level <- if (out$is_ambig) "species" else "molecular_species"

  parsed <- lapply(chain_parts, .parse_lm_chain)
  empty  <- list(raw = NA_character_, length = NA_real_, unsat = NA_real_,
                 sphingoid_prefix = NA_character_, hydroxyl_n = NA_real_)
  while (length(parsed) < 4) parsed[[length(parsed) + 1]] <- empty

  # Sphingoid base prefix (SP category): d = dihydroxy, m = monohydroxy, t = trihydroxy
  if (!is.na(out$lm_category) && out$lm_category == "SP") {
    sph <- parsed[[1]]$sphingoid_prefix
    out$sphingoid_prefix <- if (!is.na(sph)) sph else sphingoid_default
    out$hydroxyl_base    <- switch(
      if (!is.na(out$sphingoid_prefix)) out$sphingoid_prefix else "",
      "d" = 2, "m" = 1, "t" = 3, NA_real_)
  }

  out$chain1 <- if (n >= 1) chain_parts[1] else NA_character_
  out$chain2 <- if (n >= 2) chain_parts[2] else NA_character_
  out$chain3 <- if (n >= 3) chain_parts[3] else NA_character_
  out$chain4 <- if (n >= 4) chain_parts[4] else NA_character_
  out$l_1 <- parsed[[1]]$length; out$s_1 <- parsed[[1]]$unsat
  out$l_2 <- parsed[[2]]$length; out$s_2 <- parsed[[2]]$unsat
  out$l_3 <- parsed[[3]]$length; out$s_3 <- parsed[[3]]$unsat
  out$l_4 <- parsed[[4]]$length; out$s_4 <- parsed[[4]]$unsat
  out$n_chains <- n
  out$total_cl <- sum(parsed[[1]]$length, parsed[[2]]$length,
                      parsed[[3]]$length, parsed[[4]]$length, na.rm = TRUE)
  out$total_cs <- sum(parsed[[1]]$unsat, parsed[[2]]$unsat,
                      parsed[[3]]$unsat, parsed[[4]]$unsat, na.rm = TRUE)
  if (is.na(out$l_1) && out$n_chains > 0) out$not_matched <- TRUE
  out
}


# =============================================================================
# INDIVIDUAL CHAIN PARSER
# Handles: "18:1", "d18:1", "m18:0", "t18:0", "18:1(9Z)", "18:1(d7)",
#          "22:6(4Z,7Z,...)", "9-O-16:0" (FAHFA ester position)
# Note: "0:0" empty positions are filtered before reaching this function
# =============================================================================

.parse_lm_chain <- function(s) {
  s <- trimws(s)
  # Detect sphingoid base prefix: d = dihydroxy, m = monohydroxy, t = trihydroxy
  sph_m   <- regmatches(s, regexpr("^([dmt])(?=\\d)", s, perl = TRUE))
  sph_pfx <- if (length(sph_m) > 0 && nchar(sph_m) > 0) sph_m else NA_character_
  hydroxy <- switch(if (!is.na(sph_pfx)) sph_pfx else "",
                    "d" = 2, "m" = 1, "t" = 3, NA_real_)
  # Strip prefixes for numeric parsing
  s_num <- sub("^[dmt](?=\\d)", "", s,     perl = TRUE)  # sphingoid prefix
  s_num <- sub("^C(?=\\d)",     "", s_num, perl = TRUE)  # legacy C-prefix
  s_num <- sub("\\(.*\\)",      "", s_num)               # annotations: (9Z), (d7)
  s_num <- sub("^\\d+-[OP]-",   "", s_num)               # FAHFA ester position: 9-O-
  s_num <- trimws(sub("[^0-9:].*$", "", s_num))
  m <- regmatches(s_num, regexec("^(\\d+):(\\d+)", s_num))[[1]]
  list(raw = s,
       length = if (length(m) == 3) as.numeric(m[2]) else NA_real_,
       unsat  = if (length(m) == 3) as.numeric(m[3]) else NA_real_,
       sphingoid_prefix = sph_pfx, hydroxyl_n = hydroxy)
}


# =============================================================================
# CANONICAL SHORTHAND BUILDER (Liebisch et al. 2020)
# =============================================================================

.build_shorthand <- function(row, sn_confirmed = FALSE, lyso_explicit = FALSE) {
  cls <- row$Class
  if (is.na(cls)) return(NA_character_)
  # Lipids with no chains (Cholesterol, bile acids, mediators, etc.)
  if (!is.na(row$n_chains) && row$n_chains == 0) return(cls)
  chain_strs   <- .format_chains(row)
  ether_prefix <- if (!is.na(row$ether_link)) paste0(row$ether_link, "-") else ""
  # Lyso: LPC(18:1) or LPC(18:1/0:0) if lyso_explicit
  if (isTRUE(row$is_lyso)) {
    inner <- if (lyso_explicit && !is.na(chain_strs[1])) {
      paste0(ether_prefix, chain_strs[1], "/0:0")
    } else { paste0(ether_prefix, chain_strs[1]) }
    return(paste0(cls, "(", inner, ")"))
  }
  # Species level (sum composition): PC(34:1)
  if (isTRUE(row$is_ambig))
    return(paste0(cls, "(", ether_prefix, chain_strs[1], ")"))
  # Molecular species level: PC(16:0/18:1)
  valid <- chain_strs[!is.na(chain_strs)]
  inner <- if (sn_confirmed) {
    paste0(ether_prefix, paste(mapply(function(ch, n) paste0("sn", n, "=", ch),
                                      valid, seq_along(valid)), collapse = "/"))
  } else { paste0(ether_prefix, paste(valid, collapse = "/")) }
  paste0(cls, "(", inner, ")")
}


# =============================================================================
# CHAIN FORMATTER — canonical shorthand per chain
# =============================================================================

.format_chains <- function(row) {
  raw_chains <- c(row$chain1, row$chain2, row$chain3, row$chain4)
  lens   <- c(row$l_1, row$l_2, row$l_3, row$l_4)
  unsats <- c(row$s_1, row$s_2, row$s_3, row$s_4)
  result <- vector("character", 4)
  for (i in seq_len(4)) {
    raw <- raw_chains[i]
    if (is.na(raw)) { result[i] <- NA_character_; next }
    s <- trimws(raw)
    # Preserve canonical position annotations: (9Z), (4Z,7Z,...), (d7) for ISTD
    annot <- ""
    annot_m <- regmatches(s, regexpr("\\([^)]+\\)", s))
    if (length(annot_m) > 0 && nchar(annot_m) > 0 &&
        grepl("^\\(([dmt]\\d+|\\d+[EZ](,\\d+[EZ])*)\\)$", annot_m))
      annot <- annot_m
    # Sphingoid prefix for chain1 in SP (with mammalian default)
    sph_pfx <- ""
    if (!is.na(row$lm_category) && row$lm_category == "SP" && i == 1)
      sph_pfx <- if (!is.na(row$sphingoid_prefix)) row$sphingoid_prefix else ""
    # FA estolides (FAHFA): preserve ester position "9-O-16:0"
    if (!is.na(row$lm_category) && row$lm_category == "FA" &&
        grepl("^\\d+-[OP]-", s)) {
      result[i] <- trimws(sub("^C(?=\\d)", "", sub("\\(.*\\)", "", s), perl = TRUE))
      next
    }
    # Reconstruct clean chain string
    s_num <- sub("^[dmt](?=\\d)", "", s,     perl = TRUE)
    s_num <- sub("^C(?=\\d)",     "", s_num, perl = TRUE)
    s_num <- sub("\\(.*\\)",      "", s_num)
    s_num <- sub("^\\d+-[OP]-",   "", s_num)
    s_num <- trimws(sub("[^0-9:].*$", "", s_num))
    nm <- regmatches(s_num, regexec("^(\\d+):(\\d+)", s_num))[[1]]
    result[i] <- if (length(nm) == 3) paste0(sph_pfx, nm[2], ":", nm[3], annot)
    else if (!is.na(lens[i]) && !is.na(unsats[i])) paste0(sph_pfx, lens[i], ":", unsats[i])
    else NA_character_
  }
  result
}


# =============================================================================
# UTILITIES
# =============================================================================

# Null-coalescing operator: returns b if a is NA or length 0
`%||%` <- function(a, b) if (length(a) > 0 && !is.na(a[1])) a else b
