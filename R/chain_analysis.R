# R/chain_analysis.R
# Chain length and unsaturation analysis for lipidomics data.
# Source: ChainAnalysis_v7_unified_weighted.R
#
# Key changes from the original script:
#   1. Global variables (CLS_SN2, CLS_NACYL, etc.) -> arguments with defaults.
#   2. All parsers renamed to .parse_*() -- internal, not exported.
#   3. dispatch() -> .dispatch_chain(), receives class config as argument.
#   4. parse_dataset() -> parse_lipid_chains(), exported with full roxygen2.
#   5. generate_plots() -> plot_chains(), exported; stores ggplot objects,
#      does not write files directly.
#   6. CONFIG block (INPUT_FILE, CASE_LBL, etc.) eliminated entirely.

# -- Default class configuration -----------------------------------------------

#' Default chain analysis class configuration
#'
#' Returns the default list that maps lipid classes to their parsing strategy.
#' Pass the output of this function as the \code{cls_config} argument of
#' \code{parse_lipid_chains()} to override individual entries.
#'
#' @return Named list with elements \code{sn2}, \code{nacyl}, \code{long},
#'   \code{single}, and \code{excl}.
#'
#' @export
default_chain_config <- function() {
  list(
    # sn-2 chain: PC, PE, PE O
    # sn-1 relatively conserved; sn-2 is the biologically variable position.
    sn2    = c("PC", "PE", "PE O"),
    # N-acyl chain: sphingolipids with resolved "/" or "_" notation.
    # GlcCer, Hex2Cer, Hex3Cer follow the same logic as HexCer: sphingoid
    # base (d18:1) is conserved; N-acyl chain drives species diversity.
    nacyl  = c("SM", "Cer", "HexCer", "GlcCer", "Hex2Cer", "Hex3Cer"),
    # Long format: one row per acyl chain.
    # TG: 3 chains; DG, PS, PG, PA: 2 chains; CL: 4 chains (dedicated parser).
    # PI has its own dedicated parser (.parse_pi) and is not listed here.
    long   = c("TG", "DG", "PS", "PG", "PA"),
    # Single chain: lyso-species, acylcarnitines, free fatty acids, cholesteryl esters.
    single = c("LPC", "LPE", "LPI", "LPG", "LPA", "LPS", "CAR", "FFA", "FA", "CE"),
    # Excluded: ether lipids with ambiguous chain assignment, or no acyl chain.
    excl   = c("PC O", "PS O", "PG O", "TG O", "DG O",
               "MG", "CerPE", "SHexCer", "PE-Cer",
               "FC", "ST", "CoQ", "DGDG", "MGDG", "FAHFA", "Unknown")
  )
}

# -- Low-level parsers (all @noRd, not exported) -------------------------------

#' Parse a single "xx:y" chain token
#'
#' Strips modification suffixes (e.g. ";O2", "(OH)") before matching.
#' Returns a named integer vector \code{c(cl = <length>, cs = <saturation>)},
#' or \code{c(cl = NA_integer_, cs = NA_integer_)} if no match.
#'
#' @noRd
.parse_one <- function(s) {
  s <- sub(";.*$",     "", s, perl = TRUE)   # strip modification (e.g. ";O2")
  s <- sub("\\(.*\\)", "", s, perl = TRUE)   # strip annotation in parens
  s <- trimws(s)
  m <- regmatches(s, regexpr("(\\d+):(\\d+)", s, perl = TRUE))
  if (length(m) == 0L)
    return(c(cl = NA_integer_, cs = NA_integer_))
  c(cl = as.integer(sub("(\\d+):.*",  "\\1", m, perl = TRUE)),
    cs = as.integer(sub(".*:(\\d+)$", "\\1", m, perl = TRUE)))
}

#' Parse sn-2 chain from a resolved "sn1/sn2" or "sn1_sn2" annotation
#'
#' Handles double bond geometry notation such as \code{16:1(7Z)/18:1(11E)}
#' and \code{18:2(9Z.12Z)/16:0} by stripping parenthetical annotations
#' before splitting on the chain separator.
#' Returns a one-row data.frame with columns \code{chain_type}, sn-position
#' chain lengths/saturations, total chain, and analysis chain (sn-2).
#' Returns \code{NULL} if the name does not contain a resolved two-chain
#' annotation.
#'
#' @noRd
.parse_sn2 <- function(name) {
  # Match two chain tokens separated by / or _, allowing optional parenthetical
  # geometry annotations between the chain and the separator:
  # e.g. "16:1(7Z)/18:1(11E)" or "18:2(9Z.12Z)_16:0"
  m <- regmatches(name,
                  regexpr("(\\d+:\\d+(?:\\([^)]*\\))?)[_/](\\d+:\\d+(?:\\([^)]*\\))?)",
                          name, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  parts <- strsplit(m, "[_/](?=\\d)", perl = TRUE)[[1L]]
  if (length(parts) < 2L) return(NULL)
  sn1 <- .parse_one(parts[1L])
  sn2 <- .parse_one(parts[2L])
  if (anyNA(c(sn1, sn2))) return(NULL)
  data.frame(
    chain_type         = "sn2",
    sn1_cl             = sn1["cl"],  sn1_cs             = sn1["cs"],
    sn2_cl             = sn2["cl"],  sn2_cs             = sn2["cs"],
    total_cl           = sn1["cl"] + sn2["cl"],
    total_cs           = sn1["cs"] + sn2["cs"],
    analysis_chain_cl  = sn2["cl"],  analysis_chain_cs  = sn2["cs"],
    row.names          = NULL, stringsAsFactors = FALSE
  )
}

#' Parse N-acyl chain from sphingolipid notation
#'
#' Accepts both "/" (sn-resolved) and "_" (chains identified, sn-position
#' unassigned) separators. For sphingolipids the sn-position of the N-acyl
#' chain is unambiguous, so "_" is biologically equivalent to "/".
#' Returns \code{NULL} if the name has no resolved chain separator.
#'
#' @noRd
.parse_nacyl <- function(name) {
  if (!grepl("[/_]", name, perl = TRUE)) return(NULL)
  m <- regmatches(
    name,
    regexpr("(\\d+:\\d+(?:;[^/_]+)?)\\s*[/_]\\s*(\\d+:\\d+)",
            name, perl = TRUE)
  )
  if (length(m) == 0L) return(NULL)
  parts <- strsplit(m, "\\s*[/_]\\s*")[[1L]]
  base  <- .parse_one(parts[1L])
  nacyl <- .parse_one(parts[2L])
  if (anyNA(c(base, nacyl))) return(NULL)
  data.frame(
    chain_type         = "nacyl",
    base_cl            = base["cl"],  base_cs            = base["cs"],
    nacyl_cl           = nacyl["cl"], nacyl_cs           = nacyl["cs"],
    total_cl           = base["cl"] + nacyl["cl"],
    total_cs           = base["cs"] + nacyl["cs"],
    analysis_chain_cl  = nacyl["cl"], analysis_chain_cs  = nacyl["cs"],
    row.names          = NULL, stringsAsFactors = FALSE
  )
}

#' Parse total chain notation (single "xx:y" token)
#' @noRd
.parse_total <- function(name) {
  m <- regmatches(name, regexpr("(\\d+):(\\d+)", name, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  ch <- .parse_one(m)
  if (anyNA(ch)) return(NULL)
  data.frame(
    chain_type         = "total",
    total_cl           = ch["cl"],  total_cs           = ch["cs"],
    analysis_chain_cl  = ch["cl"],  analysis_chain_cs  = ch["cs"],
    row.names          = NULL, stringsAsFactors = FALSE
  )
}

#' Parse single acyl chain (lyso-species, acylcarnitines)
#' @noRd
.parse_single <- function(name) {
  m <- regmatches(name, regexpr("(\\d+):(\\d+)", name, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  ch <- .parse_one(m)
  if (anyNA(ch)) return(NULL)
  data.frame(
    chain_type         = "single",
    total_cl           = ch["cl"],  total_cs           = ch["cs"],
    analysis_chain_cl  = ch["cl"],  analysis_chain_cs  = ch["cs"],
    row.names          = NULL, stringsAsFactors = FALSE
  )
}

#' Parse TG (and other long-format lipids) into one row per acyl chain
#'
#' Excludes total-notation TG (TAG + space + digit) and ether TG O-.
#' Skips the first \code{cl:cs} token when it matches the total composition
#' (i.e. appears before any parenthesis), and filters out empty acyl positions
#' (0:0) commonly used by some software to pad DG to three positions.
#' Returns a multi-row data.frame (one row per resolved chain), or
#' \code{NULL} if the name is total notation, ether, or has < 2 valid chains.
#'
#' @noRd
.parse_long <- function(name) {
  if (grepl("^TAG\\s+\\d", name, perl = TRUE)) return(NULL)  # total notation
  if (grepl("O-",          name, fixed = TRUE)) return(NULL)  # ether lipid

  # Strip everything before the first "(" to remove the total composition token
  # e.g. "DG 30:0 (14:0/16:0/0:0)" -> "(14:0/16:0/0:0)"
  # If no parenthesis, fall back to all tokens (e.g. "DG 16:0/18:1")
  name_inner <- if (grepl("\\(", name, fixed = TRUE))
    sub("^[^(]*\\((.*)\\).*$", "\\1", name, perl = TRUE)
  else
    name

  chains_raw <- regmatches(name_inner,
                           gregexpr("\\d+:\\d+", name_inner, perl = TRUE))[[1L]]

  # Filter out empty acyl positions (0:0) used as sn-3 placeholder in DG
  chains_raw <- chains_raw[chains_raw != "0:0"]

  if (length(chains_raw) < 2L) return(NULL)

  rows <- Filter(Negate(is.null), lapply(chains_raw, function(ch) {
    p <- .parse_one(ch)
    if (anyNA(p)) return(NULL)
    data.frame(
      chain_type         = "long_format",
      source_chain       = ch,
      analysis_chain_cl  = p["cl"],
      analysis_chain_cs  = p["cs"],
      row.names          = NULL, stringsAsFactors = FALSE
    )
  }))
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

#' Parse PI: dedicated handler for phosphoinositols
#'
#' PI 38:4 (total notation) -> excluded with status
#' \code{"excluded_total_notation_unresolved_PI"}.
#' PI 18:0/20:4 or 18:0_20:4 -> long_format, two rows (one per chain).
#' Both PI acyl chains are biologically variable; only resolved annotations
#' enter chain analysis.
#'
#' @return Named list with elements \code{status} (character) and
#'   \code{data} (data.frame or \code{NULL}).
#'
#' @noRd
.parse_pi <- function(name) {
  has_two_chains <- grepl("(\\d+:\\d+)[_/](\\d+:\\d+)", name, perl = TRUE)
  if (!has_two_chains)
    return(list(status = "excluded_total_notation_unresolved_PI", data = NULL))

  m     <- regmatches(name,
                      regexpr("(\\d+:\\d+)[_/](\\d+:\\d+)", name, perl = TRUE))
  parts <- strsplit(m, "[_/]")[[1L]]
  ch1   <- .parse_one(parts[1L])
  ch2   <- .parse_one(parts[2L])

  if (anyNA(c(ch1, ch2)))
    return(list(status = "excluded_unresolved_PI", data = NULL))

  data <- data.frame(
    chain_type         = "long_format",
    source_chain       = parts,
    analysis_chain_cl  = c(ch1["cl"], ch2["cl"]),
    analysis_chain_cs  = c(ch1["cs"], ch2["cs"]),
    total_cl           = ch1["cl"] + ch2["cl"],
    total_cs           = ch1["cs"] + ch2["cs"],
    row.names          = NULL, stringsAsFactors = FALSE
  )
  list(status = "parsed_PI_resolved", data = data)
}

#' Parse CL (cardiolipin) into one row per acyl chain
#'
#' CL carries four acyl chains (two per phosphatidylglycerol arm).
#' CL 72:8 (total notation, single token) -> excluded with status
#' \code{"excluded_total_notation_unresolved_CL"}.
#' CL 16:1/16:1/18:1/14:0 -> long_format, four rows (one per chain).
#' All four acyl positions are biologically variable; only resolved
#' annotations enter chain analysis.
#'
#' @return Named list with elements \code{status} (character) and
#'   \code{data} (data.frame or \code{NULL}).
#'
#' @noRd
.parse_cl <- function(name) {
  chains_raw <- regmatches(name,
                           gregexpr("\\d+:\\d+", name, perl = TRUE))[[1L]]
  if (length(chains_raw) < 2L)
    return(list(status = "excluded_total_notation_unresolved_CL", data = NULL))

  rows <- Filter(Negate(is.null), lapply(chains_raw, function(ch) {
    p <- .parse_one(ch)
    if (anyNA(p)) return(NULL)
    data.frame(
      chain_type         = "long_format",
      source_chain       = ch,
      analysis_chain_cl  = p["cl"],
      analysis_chain_cs  = p["cs"],
      row.names          = NULL, stringsAsFactors = FALSE
    )
  }))
  if (length(rows) == 0L)
    return(list(status = "excluded_unresolved_CL", data = NULL))

  list(status = "parsed_CL_resolved",
       data   = do.call(rbind, rows))
}

# -- Dispatcher ----------------------------------------------------------------

#' Route a single lipid to the correct parser
#'
#' @param name    Character(1). Primary lipid name.
#' @param shorthand Character(1). Shorthand notation (fallback for sn2/single
#'   when the primary name is a common name without ":").
#' @param cls     Character(1). Already-assigned \code{LipidClass}.
#' @param cfg     Named list from \code{default_chain_config()}.
#'
#' @return Named list: \code{type} (character status/chain_type) and
#'   \code{data} (data.frame or \code{NULL}).
#'
#' @noRd
.dispatch_chain <- function(name, shorthand, cls, cfg) {

  if (cls %in% cfg$excl)
    return(list(type = "excluded", data = NULL))

  if (cls %in% cfg$sn2) {
    parsed <- .parse_sn2(name)
    if (is.null(parsed)) parsed <- .parse_sn2(shorthand)   # common-name fallback
    return(list(type = "sn2", data = parsed))
  }

  if (cls %in% cfg$nacyl)
    return(list(type = "nacyl", data = .parse_nacyl(name)))

  if (cls == "PI") {
    result <- .parse_pi(name)
    if (is.null(result$data))
      return(list(type = result$status, data = NULL))
    return(list(type = "long_format", data = result$data))
  }

  if (cls == "CL") {
    result <- .parse_cl(name)
    if (is.null(result$data))
      return(list(type = result$status, data = NULL))
    return(list(type = "long_format", data = result$data))
  }

  if (cls %in% cfg$long)
    return(list(type = "long_format", data = .parse_long(name)))

  if (cls %in% cfg$single) {
    src <- if (!grepl(":", name, fixed = TRUE)) shorthand else name
    return(list(type = "single", data = .parse_single(src)))
  }

  list(type = "excluded", data = NULL)
}

# -- Main exported function ----------------------------------------------------

#' Parse acyl chain composition from a lipidomics data.frame
#'
#' Applies biology-aware chain parsing to each lipid in \code{data},
#' routing each species to the appropriate parser based on its lipid class:
#' sn-2 (PC, PE, PE O), N-acyl (SM, Cer, HexCer, GlcCer, Hex2Cer, Hex3Cer),
#' long-format (TG, DG, PS, PG, PA, PI, CL),
#' single-chain (LPC, LPE, LPI, LPG, LPA, LPS, CAR, FFA, FA, CE), or excluded.
#'
#' @param data A \code{data.frame} with at least the columns specified by
#'   \code{lipid_col} and \code{class_col}. Typically the output of
#'   \code{annotate_lipids()}.
#' @param lipid_col Character(1). Name of the lipid identifier column.
#'   Default: \code{"LipidName"}.
#' @param class_col Character(1). Name of the lipid class column (must contain
#'   abbreviated class names such as "PC", "TG", "SM"). Default:
#'   \code{"LipidClass"}.
#' @param shorthand_col Character(1) or \code{NULL}. Name of an optional
#'   shorthand column used as fallback for sn-2 and single-chain parsing when
#'   the primary name is a common name. Default: \code{"Shorthand"}.
#' @param rank_col Character(1) or \code{NULL}. Name of a confidence-rank
#'   column. Rows with rank \code{"P"} or \code{NA} are excluded from
#'   analysis. Set to \code{NULL} to skip rank filtering. Default:
#'   \code{"Confidence_rank"}.
#' @param cls_config Named list from \code{\link{default_chain_config}}.
#'   Override individual elements to change class routing.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{\code{parsed}}{Long-format \code{data.frame} with one row per
#'       chain observation. Contains all columns from \code{data} plus chain
#'       fields (\code{analysis_chain_cl}, \code{analysis_chain_cs},
#'       \code{chain_type}, etc.).}
#'     \item{\code{summary}}{Per-lipid parsing log \code{data.frame} with
#'       columns \code{LipidName}, \code{LipidClass},
#'       \code{Confidence_rank}, \code{status}, and \code{chain_type}.}
#'     \item{\code{wide}}{Wide-format \code{data.frame} with one row per
#'       lipid. Columns \code{sn1}, \code{sn2}, \code{sn3}, \code{sn4}
#'       contain individual acyl chain positions (e.g. \code{"18:1"}), and
#'       \code{total_carbons} and \code{total_unsat} give the summed totals.
#'       For cardiolipins (CL), all four sn positions are populated.
#'       For sphingolipids, \code{sn1} = sphingoid base, \code{sn2} =
#'       N-acyl chain.}
#'   }
#'
#' @seealso \code{\link{default_chain_config}}, \code{plot_chains()}
#'
#' @examples
#' \dontrun{
#' data("lipid_example", package = "easyLSEA")
#' annotated <- annotate_lipids(lipid_example)
#' chains <- parse_lipid_chains(annotated)
#' head(chains$parsed)
#' head(chains$summary)
#' head(chains$wide)
#' }
#'
#' @export
parse_lipid_chains <- function(
    data,
    lipid_col     = "LipidName",
    class_col     = "LipidClass",
    shorthand_col = "Shorthand",
    rank_col      = "Confidence_rank",
    cls_config    = default_chain_config()
) {

  # -- Input checks ------------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)
  for (col in c(lipid_col, class_col)) {
    if (!col %in% names(data))
      stop("Column ", sQuote(col), " not found in 'data'.", call. = FALSE)
  }
  if (!is.null(shorthand_col) && !shorthand_col %in% names(data)) {
    message("Column ", sQuote(shorthand_col),
            " not found -- shorthand fallback disabled.")
    shorthand_col <- NULL
  }
  if (!is.null(rank_col) && !rank_col %in% names(data)) {
    message("Column ", sQuote(rank_col),
            " not found -- rank filtering skipped.")
    rank_col <- NULL
  }

  # Columns to carry through to the parsed output
  base_cols <- intersect(
    c(lipid_col, class_col, shorthand_col, rank_col,
      "logFC", "AveExpr", "weight",
      "t", "P.Value", "adj.P.Val", "B", "sig",
      "Case", "Reference", "Comparison",
      "LipidClass_Full", "LipidCategory",
      "LipidCategory_LMAPS", "LipidCategory_functional"),
    names(data)
  )

  all_chain_cols <- c(
    "sn1_cl", "sn1_cs", "sn2_cl", "sn2_cs",
    "base_cl", "base_cs", "nacyl_cl", "nacyl_cs",
    "source_chain", "total_cl", "total_cs"
  )

  parsed_rows  <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(data))) {

    row       <- data[i, , drop = FALSE]
    lip_name  <- as.character(row[[lipid_col]])
    lip_cls   <- as.character(row[[class_col]])
    lip_short <- if (!is.null(shorthand_col)) as.character(row[[shorthand_col]]) else lip_name
    lip_rank  <- if (!is.null(rank_col))      as.character(row[[rank_col]])      else "A"

    # Always exclude rank "P" or NA (lowest confidence)
    if (!is.null(rank_col) && (is.na(lip_rank) || lip_rank == "P")) {
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        LipidName       = lip_name,
        LipidClass      = lip_cls,
        Confidence_rank = lip_rank,
        status          = "excluded_rank_P_or_NA",
        chain_type      = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    result <- .dispatch_chain(lip_name, lip_short, lip_cls, cls_config)

    # Determine parsing status -- propagate class-specific exclusion labels
    status <- if (result$type == "excluded") {
      "excluded_class"
    } else if (is.null(result$data)) {
      if      (lip_cls == "SM")                            "excluded_sm_unresolved"
      else if (startsWith(result$type, "excluded_"))       result$type
      else                                                  "excluded_unresolved"
    } else {
      "parsed"
    }

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      LipidName       = lip_name,
      LipidClass      = lip_cls,
      Confidence_rank = lip_rank,
      status          = status,
      chain_type      = if (status == "parsed") result$type else NA_character_,
      stringsAsFactors = FALSE
    )

    if (status != "parsed") next

    parsed     <- result$data
    n_rows_out <- nrow(parsed)

    # Ensure all chain columns present (fill missing with NA)
    for (col in all_chain_cols) {
      if (!col %in% names(parsed)) parsed[[col]] <- NA
    }

    base_rep           <- data[rep(i, n_rows_out), base_cols, drop = FALSE]
    rownames(base_rep) <- NULL
    rownames(parsed)   <- NULL

    parsed_rows[[length(parsed_rows) + 1L]] <- cbind(
      base_rep, parsed, stringsAsFactors = FALSE
    )
  }

  parsed_df  <- if (length(parsed_rows)  > 0L) do.call(rbind, parsed_rows)  else data.frame()
  summary_df <- if (length(summary_rows) > 0L) do.call(rbind, summary_rows) else data.frame()
  wide_df    <- .build_wide(parsed_df, lipid_col = lipid_col)

  list(
    parsed  = parsed_df,
    summary = summary_df,
    wide    = wide_df
  )
}


# -- Wide summary builder ------------------------------------------------------

#' Build a wide-format chain summary (one row per lipid)
#'
#' Collapses the long-format parsed output into one row per lipid, with
#' individual sn positions as columns (sn1, sn2, sn3, sn4 for CL) and
#' total carbon/unsaturation counts.
#'
#' Chain type mapping to sn columns:
#'   sn2        -> sn1 (from sn1_cl/sn1_cs), sn2 (from sn2_cl/sn2_cs)
#'   nacyl      -> sn_base (sphingoid base), sn_nacyl (N-acyl chain)
#'   total      -> total only (individual positions not resolved)
#'   single     -> sn1 only
#'   long_format-> sn1, sn2, sn3, sn4 by order of appearance
#'
#' @param parsed data.frame. Output of parse_lipid_chains()$parsed.
#' @param lipid_col Character(1). Lipid name column. Default: "LipidName".
#' @return data.frame with one row per lipid.
#'
#' @noRd
.build_wide <- function(parsed, lipid_col = "LipidName") {

  if (is.null(parsed) || nrow(parsed) == 0L) return(data.frame())

  # Base columns to carry through (one value per lipid)
  base_cols <- intersect(
    c(lipid_col, "LipidClass", "Confidence_rank", "chain_type",
      "logFC", "adj.P.Val", "sig"),
    names(parsed)
  )

  # Helper: format one chain as "CL:CS" string, e.g. "18:1"
  .fmt <- function(cl, cs) {
    cl <- suppressWarnings(as.numeric(cl))
    cs <- suppressWarnings(as.numeric(cs))
    ifelse(is.na(cl) | is.na(cs), NA_character_,
           paste0(cl, ":", cs))
  }

  # Split by lipid name and collapse each group
  lips <- split(parsed, parsed[[lipid_col]])

  rows <- lapply(lips, function(df) {
    base  <- df[1L, base_cols, drop = FALSE]
    ctype <- df$chain_type[[1L]]
    n     <- nrow(df)

    # Total carbons and unsaturation
    # For long_format (TG, PI, CL) total_cl may be NA — sum from individual chains
    total_cl <- if ("total_cl" %in% names(df) && !all(is.na(df$total_cl))) {
      df$total_cl[[1L]]
    } else if ("analysis_chain_cl" %in% names(df) && !all(is.na(df$analysis_chain_cl))) {
      sum(as.numeric(df$analysis_chain_cl), na.rm = TRUE)
    } else NA_integer_

    total_cs <- if ("total_cs" %in% names(df) && !all(is.na(df$total_cs))) {
      df$total_cs[[1L]]
    } else if ("analysis_chain_cs" %in% names(df) && !all(is.na(df$analysis_chain_cs))) {
      sum(as.numeric(df$analysis_chain_cs), na.rm = TRUE)
    } else NA_integer_

    # Build sn position columns depending on chain type
    if (ctype == "sn2") {
      sn1   <- .fmt(df$sn1_cl[[1L]], df$sn1_cs[[1L]])
      sn2   <- .fmt(df$sn2_cl[[1L]], df$sn2_cs[[1L]])
      extra <- data.frame(sn1 = sn1, sn2 = sn2,
                          sn3 = NA_character_, sn4 = NA_character_,
                          stringsAsFactors = FALSE)

    } else if (ctype == "nacyl") {
      sn_base  <- .fmt(df$base_cl[[1L]],  df$base_cs[[1L]])
      sn_nacyl <- .fmt(df$nacyl_cl[[1L]], df$nacyl_cs[[1L]])
      extra <- data.frame(sn1 = sn_base, sn2 = sn_nacyl,
                          sn3 = NA_character_, sn4 = NA_character_,
                          stringsAsFactors = FALSE)

    } else if (ctype == "long_format") {
      # Chains ordered by appearance (sn1..sn4)
      chains <- .fmt(df$analysis_chain_cl, df$analysis_chain_cs)
      sn1 <- if (n >= 1L) chains[[1L]] else NA_character_
      sn2 <- if (n >= 2L) chains[[2L]] else NA_character_
      sn3 <- if (n >= 3L) chains[[3L]] else NA_character_
      sn4 <- if (n >= 4L) chains[[4L]] else NA_character_
      extra <- data.frame(sn1 = sn1, sn2 = sn2, sn3 = sn3, sn4 = sn4,
                          stringsAsFactors = FALSE)

    } else {
      # total or single — no individual positions
      extra <- data.frame(sn1 = NA_character_, sn2 = NA_character_,
                          sn3 = NA_character_, sn4 = NA_character_,
                          stringsAsFactors = FALSE)
    }

    cbind(base,
          data.frame(total_carbons = total_cl,
                     total_unsat   = total_cs,
                     stringsAsFactors = FALSE),
          extra,
          row.names = NULL,
          stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# -- Plot helpers (internal) ---------------------------------------------------

#' @noRd
.theme_chain <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#2C3E50"),
      strip.text       = ggplot2::element_text(color = "white", face = "bold"),
      panel.grid       = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1L),
      plot.subtitle    = ggplot2::element_text(size = base_size - 2L, color = "grey40"),
      legend.position  = "right"
    )
}

#' @noRd
.axis_lbl <- function(chain_type, axis = c("cl", "cs")) {
  axis <- match.arg(axis)
  lbl  <- list(
    sn2         = c(cl = "sn-2 chain length (carbons)",
                    cs = "sn-2 unsaturation (double bonds)"),
    nacyl       = c(cl = "N-acyl chain length (carbons)",
                    cs = "N-acyl unsaturation (double bonds)"),
    total       = c(cl = "Total chain length (carbons)",
                    cs = "Total unsaturation (double bonds)"),
    long_format = c(cl = "Acyl chain length (carbons)",
                    cs = "Acyl chain unsaturation (double bonds)"),
    single      = c(cl = "Chain length (carbons)",
                    cs = "Unsaturation (double bonds)")
  )
  if (chain_type %in% names(lbl)) return(lbl[[chain_type]][[axis]])
  if (axis == "cl") "Chain length (carbons)" else "Unsaturation (double bonds)"
}

#' @noRd
.method_lbl <- function(chain_type) {
  switch(chain_type,
    sn2         = "sn-2 chain (sn-1 excluded)",
    nacyl       = "N-acyl chain (sphingoid base excluded)",
    total       = "Total chain (both acyl chains summed)",
    long_format = "All acyl chains -- long format, weighted by abundance",
    single      = "Single acyl chain",
    "unknown"
  )
}

# -- Exported plot function ----------------------------------------------------

#' Generate chain analysis plots
#'
#' Produces tile and trend plots for each lipid class with sufficient
#' chain observations. Returns a named list of \code{\link[ggplot2]{ggplot}}
#' objects; does not write files. Use \code{export_lsea()} to save.
#'
#' @param chains_result Named list returned by \code{\link{parse_lipid_chains}}.
#' @param case_lbl Character(1). Label for the case group. Default:
#'   \code{"Case"}.
#' @param ref_lbl Character(1). Label for the reference group. Default:
#'   \code{"Reference"}.
#' @param fdr_thresh Numeric(1). FDR threshold used to mark significant
#'   lipids in tile plot cell labels. Default: \code{0.05}.
#' @param min_n_tile Integer(1). Minimum chain observations per class to
#'   produce a tile plot. Default: \code{4L}.
#' @param min_n_trend Integer(1). Minimum chain observations per class to
#'   produce trend plots. Default: \code{5L}.
#'
#' @return Named list of \code{ggplot} objects with elements
#'   \code{tile_<CLASS>}, \code{trend_length_<CLASS>},
#'   \code{trend_unsat_<CLASS>}.
#'
#' @seealso \code{\link{parse_lipid_chains}}, \code{export_lsea()}
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#' @importFrom ggplot2 labs theme_bw theme element_rect element_text element_blank
#' @importFrom ggplot2 geom_point geom_smooth geom_hline scale_color_manual
#'
#' @export
plot_chains <- function(
    chains_result,
    case_lbl   = "Case",
    ref_lbl    = "Reference",
    fdr_thresh = 0.05,
    min_n_tile  = 4L,
    min_n_trend = 5L
) {

  df <- chains_result$parsed

  if (is.null(df) || nrow(df) == 0L) {
    message("No parsed chain data available -- returning empty list.")
    return(list())
  }

  # Observation count per class
  n_obs_cls    <- tapply(df[[1L]], df$LipidClass, length)
  classes_tile <- names(which(n_obs_cls >= min_n_tile))
  classes_trend <- names(which(n_obs_cls >= min_n_trend))

  # Determine logFC and sig columns (flexible names)
  fc_col  <- intersect(c("logFC", "log2FC"), names(df))[[1L]]
  sig_col <- intersect(c("sig", "significant"), names(df))
  has_sig <- length(sig_col) > 0L
  if (has_sig) sig_col <- sig_col[[1L]]
  wt_col  <- if ("weight" %in% names(df)) "weight" else NULL

  plots <- list()

  # -- Tile plots --------------------------------------------------------------
  for (cls in classes_tile) {
    df_cls <- df[df$LipidClass == cls, ]
    ctype  <- df_cls$chain_type[[1L]]
    n_lip  <- length(unique(df_cls[[1L]]))

    cell_key <- paste(df_cls$analysis_chain_cl,
                      df_cls$analysis_chain_cs, sep = "|")

    tile <- do.call(rbind, lapply(unique(cell_key), function(k) {
      sub <- df_cls[cell_key == k, ]
      wt  <- if (!is.null(wt_col)) sub[[wt_col]] else rep(1, nrow(sub))
      data.frame(
        cl_f       = sub$analysis_chain_cl[[1L]],
        cs_f       = sub$analysis_chain_cs[[1L]],
        mean_logFC = stats::weighted.mean(sub[[fc_col]], w = wt, na.rm = TRUE),
        n_lip_cell = length(unique(sub[[1L]])),
        n_sig_cell = if (has_sig)
                       sum(sub[[sig_col]][!duplicated(sub[[1L]])], na.rm = TRUE)
                     else 0L
      )
    }))

    tile$cl_f <- factor(tile$cl_f, levels = sort(unique(tile$cl_f)))
    tile$cs_f <- factor(tile$cs_f, levels = sort(unique(tile$cs_f)))

    p <- ggplot2::ggplot(tile, ggplot2::aes(cs_f, cl_f, fill = mean_logFC)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.4) +
      ggplot2::scale_fill_gradient2(
        midpoint = 0, low = "#2166AC", mid = "white", high = "#E63946",
        name = "weighted\nmean logFC"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = paste0(
          n_lip_cell,
          ifelse(n_sig_cell > 0,
                 paste0(" (", n_sig_cell, "*)"), "")
        )),
        size = 3, color = "grey20"
      ) +
      ggplot2::labs(
        title    = paste0(cls, " -- Chain Distribution"),
        subtitle = paste0(
          case_lbl, " vs ", ref_lbl, "  |  n = ", n_lip, " lipids  |  ",
          .method_lbl(ctype), "\n",
          "fill = weighted mean logFC  |  n (n_sig*): lipids per cell",
          if (has_sig) paste0(" (FDR < ", fdr_thresh, ")") else ""
        ),
        x = .axis_lbl(ctype, "cs"),
        y = .axis_lbl(ctype, "cl")
      ) +
      .theme_chain()

    plots[[paste0("tile_", cls)]] <- p
  }

  # -- Trend plots: chain length ------------------------------------------------
  for (cls in classes_trend) {
    df_cls <- df[df$LipidClass == cls &
                   !is.na(df$analysis_chain_cl), ]
    if (nrow(df_cls) < min_n_trend) next
    ctype <- df_cls$chain_type[[1L]]
    wt    <- if (!is.null(wt_col)) df_cls[[wt_col]] else rep(1, nrow(df_cls))

    df_agg <- do.call(rbind, lapply(
      split(df_cls, df_cls$analysis_chain_cl),
      function(sub) {
        wts <- if (!is.null(wt_col)) sub[[wt_col]] else rep(1, nrow(sub))
        data.frame(
          analysis_chain_cl = sub$analysis_chain_cl[[1L]],
          mean_logFC        = stats::weighted.mean(sub[[fc_col]], w = wts,
                                                   na.rm = TRUE),
          n                 = nrow(sub)
        )
      }
    ))

    p <- ggplot2::ggplot(df_agg,
                         ggplot2::aes(x = analysis_chain_cl,
                                      y = mean_logFC)) +
      ggplot2::geom_point(ggplot2::aes(size = n), color = "#2C3E50",
                          alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", se = TRUE,
                           color = "#E63946", fill = "#E63946",
                           alpha = 0.15, linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60", linewidth = 0.5) +
      ggplot2::labs(
        title    = paste0(cls, " -- Trend by Chain Length"),
        subtitle = paste0(case_lbl, " vs ", ref_lbl, "  |  ",
                          .method_lbl(ctype)),
        x        = .axis_lbl(ctype, "cl"),
        y        = "weighted mean log2FC"
      ) +
      .theme_chain()

    plots[[paste0("trend_length_", cls)]] <- p
  }

  # -- Trend plots: unsaturation ------------------------------------------------
  for (cls in classes_trend) {
    df_cls <- df[df$LipidClass == cls &
                   !is.na(df$analysis_chain_cs), ]
    if (nrow(df_cls) < min_n_trend) next
    ctype <- df_cls$chain_type[[1L]]

    df_agg <- do.call(rbind, lapply(
      split(df_cls, df_cls$analysis_chain_cs),
      function(sub) {
        wts <- if (!is.null(wt_col)) sub[[wt_col]] else rep(1, nrow(sub))
        data.frame(
          analysis_chain_cs = sub$analysis_chain_cs[[1L]],
          mean_logFC        = stats::weighted.mean(sub[[fc_col]], w = wts,
                                                   na.rm = TRUE),
          n                 = nrow(sub)
        )
      }
    ))

    p <- ggplot2::ggplot(df_agg,
                         ggplot2::aes(x = analysis_chain_cs,
                                      y = mean_logFC)) +
      ggplot2::geom_point(ggplot2::aes(size = n), color = "#2C3E50",
                          alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", se = TRUE,
                           color = "#2166AC", fill = "#2166AC",
                           alpha = 0.15, linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60", linewidth = 0.5) +
      ggplot2::labs(
        title    = paste0(cls, " -- Trend by Unsaturation"),
        subtitle = paste0(case_lbl, " vs ", ref_lbl, "  |  ",
                          .method_lbl(ctype)),
        x        = .axis_lbl(ctype, "cs"),
        y        = "weighted mean log2FC"
      ) +
      .theme_chain()

    plots[[paste0("trend_unsat_", cls)]] <- p
  }

  message("Generated ", length(plots), " chain analysis plot(s).")
  plots
}
