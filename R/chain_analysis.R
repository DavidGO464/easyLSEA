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

  # Exclude lyso-lipids: sn-2 position is 0:0 (unspecified acyl chain).
  # These are lyso-species (e.g. PC(20:4/0:0)) and should not enter
  # chain analysis as if they had a 0-carbon chain.
  # Also exclude sn-1 = 0:0 for the same reason.
  if (sn2["cl"] == 0L && sn2["cs"] == 0L) return(NULL)
  if (sn1["cl"] == 0L && sn1["cs"] == 0L) return(NULL)

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

  # Exclude unresolved N-acyl: 0:0 means the acyl chain is not specified
  if (nacyl["cl"] == 0L && nacyl["cs"] == 0L) return(NULL)
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
#'   column. Rows with rank \code{"P"} or \code{NA} are always excluded.
#'   Set to \code{NULL} to skip rank filtering entirely. Default:
#'   \code{"Confidence_rank"}.
#' @param min_rank Character(1). Minimum confidence rank to include in
#'   analysis. Ranks are ordered \code{A > B > C > D > E > P}. Setting
#'   \code{min_rank = "B"} includes only ranks A and B, excluding C, D, E
#'   and P. Default: \code{"E"} (include all except P and NA).
#' @param cls_config Named list from \code{\link{default_chain_config}}.
#'   Override individual elements to change class routing.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{parsed}}{Long-format \code{data.frame} with one row per
#'       chain observation. Contains all columns from \code{data} plus chain
#'       fields (\code{analysis_chain_cl}, \code{analysis_chain_cs},
#'       \code{chain_type}, etc.).}
#'     \item{\code{summary}}{Per-lipid parsing log \code{data.frame} with
#'       columns \code{LipidName}, \code{LipidClass},
#'       \code{Confidence_rank}, \code{status}, and \code{chain_type}.}
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
#' }
#'
#' @export
parse_lipid_chains <- function(
    data,
    lipid_col     = "LipidName",
    class_col     = "LipidClass",
    shorthand_col = "Shorthand",
    rank_col      = "Confidence_rank",
    min_rank      = "E",
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

  # Rank order: A is highest confidence, P is provisional (always excluded)
  .rank_order <- c(A = 1L, B = 2L, C = 3L, D = 4L, E = 5L, P = 6L)
  min_rank    <- toupper(min_rank)
  if (!min_rank %in% names(.rank_order))
    stop("'min_rank' must be one of: A, B, C, D, E. Got: ", sQuote(min_rank),
         call. = FALSE)
  min_rank_val <- .rank_order[[min_rank]]

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

    # Exclude rank "P", NA, or below min_rank threshold
    lip_rank_val <- if (!is.null(rank_col) && lip_rank %in% names(.rank_order))
      .rank_order[[lip_rank]] else NA_integer_

    if (!is.null(rank_col) && (is.na(lip_rank) || lip_rank == "P" ||
                                (!is.na(lip_rank_val) && lip_rank_val > min_rank_val))) {
      status_lbl <- if (is.na(lip_rank) || lip_rank == "P")
        "excluded_rank_P_or_NA"
      else
        paste0("excluded_rank_below_", min_rank)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        LipidName       = lip_name,
        LipidClass      = lip_cls,
        Confidence_rank = lip_rank,
        status          = status_lbl,
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
      if      (grepl("0:0", lip_name, fixed = TRUE))       "excluded_lyso_unresolved"
      else if (lip_cls == "SM")                            "excluded_sm_unresolved"
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

#' @noRd
.build_wide <- function(parsed, lipid_col = "LipidName") {

  if (is.null(parsed) || nrow(parsed) == 0L) return(data.frame())

  base_cols <- intersect(
    c(lipid_col, "LipidClass", "Confidence_rank", "chain_type",
      "logFC", "adj.P.Val", "sig"),
    names(parsed)
  )

  .fmt <- function(cl, cs) {
    cl <- suppressWarnings(as.numeric(cl))
    cs <- suppressWarnings(as.numeric(cs))
    ifelse(is.na(cl) | is.na(cs), NA_character_, paste0(cl, ":", cs))
  }

  lips <- split(parsed, parsed[[lipid_col]])

  rows <- lapply(lips, function(df) {
    base  <- df[1L, base_cols, drop = FALSE]
    ctype <- df$chain_type[[1L]]
    n     <- nrow(df)

    # Total carbons and unsaturation
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
      chains <- .fmt(df$analysis_chain_cl, df$analysis_chain_cs)
      sn1 <- if (n >= 1L) chains[[1L]] else NA_character_
      sn2 <- if (n >= 2L) chains[[2L]] else NA_character_
      sn3 <- if (n >= 3L) chains[[3L]] else NA_character_
      sn4 <- if (n >= 4L) chains[[4L]] else NA_character_
      extra <- data.frame(sn1 = sn1, sn2 = sn2, sn3 = sn3, sn4 = sn4,
                          stringsAsFactors = FALSE)
    } else {
      extra <- data.frame(sn1 = NA_character_, sn2 = NA_character_,
                          sn3 = NA_character_, sn4 = NA_character_,
                          stringsAsFactors = FALSE)
    }

    cbind(base,
          data.frame(total_carbons = total_cl, total_unsat = total_cs,
                     stringsAsFactors = FALSE),
          extra, row.names = NULL, stringsAsFactors = FALSE)
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

#' @noRd
.trend_annotation <- function(x, y, weights, test, case_lbl, ref_lbl, fdr_thresh) {
  # Runs the selected statistical test on raw (individual) observations
  # and returns a formatted annotation string for geom_label.
  if (test == "none" || length(x) < 3L) return(NULL)

  ann <- tryCatch({
    if (test == "spearman") {
      res  <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
      rho  <- round(res$estimate, 2L)
      pval <- res$p.value
      pstr <- if (pval < 0.001) "p < 0.001"
              else paste0("p = ", round(pval, 3L))
      paste0("Spearman \u03c1 = ", sprintf("%+.2f", rho), ",  ", pstr)
    } else {
      # Weighted linear regression
      df_lm <- data.frame(x = x, y = y, w = weights)
      fit   <- stats::lm(y ~ x, data = df_lm, weights = w)
      cf    <- summary(fit)$coefficients
      beta  <- round(cf["x", "Estimate"], 3L)
      pval  <- cf["x", "Pr(>|t|)"]
      pstr  <- if (pval < 0.001) "p < 0.001"
               else paste0("p = ", round(pval, 3L))
      paste0("\u03b2 = ", sprintf("%+.3f", beta), " per unit,  ", pstr,
             "  (weighted LM)")
    }
  }, error = function(e) NULL)

  ann
}



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
#' @param fdr_thresh Numeric(1). FDR threshold to colour individual lipid
#'   points in trend plots (red = FDR sig, grey = NS) and to label
#'   significant counts in tile cells. Default: \code{0.05}.
#' @param min_n_tile Integer(1). Minimum chain observations per class to
#'   produce a tile plot. Default: \code{4L}.
#' @param min_n_trend Integer(1). Minimum chain observations per class to
#'   produce trend plots. Default: \code{5L}.
#' @param smooth_method Character(1). Smoothing method for trend plots.
#'   \code{"loess"} (default) fits a local polynomial; \code{"lm"} fits a
#'   global linear model. Use \code{"lm"} for small datasets or when a
#'   monotone trend is expected a priori.
#' @param smooth_span Numeric(1). Span for loess smoothing (only used when
#'   \code{smooth_method = "loess"}). Smaller values (e.g. \code{0.4}) produce
#'   a more flexible curve; larger values (e.g. \code{0.9}) produce a smoother
#'   curve. Default: \code{0.75}. A warning is issued when
#'   \code{smooth_span < 0.5} and fewer than 10 observations are available,
#'   as this combination risks overfitting.
#' @param smooth_weighted Logical(1). If \code{TRUE} (default), the smoothing
#'   curve is weighted by the number of chain observations per x-axis position,
#'   giving more influence to well-represented chain lengths/unsaturations.
#'   Mathematically more appropriate than unweighted loess when observation
#'   counts are unequal across positions.
#' @param smooth_se Logical(1). Whether to display the 95\% confidence
#'   interval ribbon around the smoothing curve. Default: \code{TRUE}.
#' @param show_points Logical(1). Whether to display individual lipid points
#'   in trend plots, coloured by FDR significance. Default: \code{TRUE}.
#'   Set to \code{FALSE} to show only the smoothing curve (cleaner for
#'   classes with many lipids).
#' @param trend_test Character(1). Statistical test to annotate on trend plots.
#'   \code{"spearman"} (default) computes Spearman rank correlation between
#'   chain position (length or unsaturation) and logFC across individual lipids,
#'   reporting \eqn{\rho} and p-value. \code{"lm"} fits a weighted linear
#'   regression (weighted by n observations per position) and reports the slope
#'   \eqn{\beta} and p-value. \code{"none"} shows no statistical annotation.
#'   Note: these tests are computed on the individual lipid observations, not
#'   on the smoothed curve.
#' @param trend_x_step_length Integer(1) or \code{NULL}. Step size for
#'   x-axis tick marks in chain length trend plots. Default: \code{2L}
#'   (every 2 carbons), suitable for the typical range of 8--36 carbons.
#'   Use \code{1L} for fine-grained resolution or \code{4L} for very wide
#'   ranges. When \code{NULL}, ggplot2 chooses breaks automatically.
#' @param trend_x_step_unsat Integer(1) or \code{NULL}. Step size for
#'   x-axis tick marks in unsaturation trend plots. Default: \code{1L}
#'   (every double bond), suitable for the typical range of 0--8.
#'   When \code{NULL}, ggplot2 chooses breaks automatically.
#' @param tile_label Character(1). What to display inside each tile cell:
#'   \code{"both"} (default) shows total and significant lipid counts;
#'   \code{"n"} shows only the total; \code{"sig"} shows only significant;
#'   \code{"none"} shows no text.
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
    case_lbl        = "Case",
    ref_lbl         = "Reference",
    fdr_thresh      = 0.05,
    min_n_tile      = 4L,
    min_n_trend     = 5L,
    smooth_method   = c("loess", "lm"),
    smooth_span     = 0.75,
    smooth_weighted = TRUE,
    smooth_se       = TRUE,
    show_points     = TRUE,
    tile_label      = c("both", "n", "sig", "none"),
    trend_test          = c("spearman", "lm", "none"),
    trend_x_step_length = 2L,
    trend_x_step_unsat  = 1L
) {

  smooth_method <- match.arg(smooth_method)
  tile_label    <- match.arg(tile_label)
  trend_test    <- match.arg(trend_test)

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
    ctype  <- df_cls$chain_type[[1L]]
    n_lip  <- length(unique(df_cls[[1L]]))
    n_obs  <- nrow(df_cls)

    # Overfitting warning: small span + few observations
    if (smooth_method == "loess" && smooth_span < 0.5 && n_obs < 10L)
      warning(cls, " trend_length: smooth_span=", smooth_span,
              " with only ", n_obs, " observations may overfit. ",
              "Consider increasing smooth_span.", call. = FALSE)

    # Individual lipid logFC (one row per chain observation)
    df_pts <- df_cls
    if (has_sig) {
      df_pts$sig_label <- ifelse(df_pts[[sig_col]] == 1L,
                                 paste0("FDR < ", fdr_thresh), "NS")
    } else {
      df_pts$sig_label <- "NS"
    }

    # Per-position weights for the smooth (n chain obs per x value)
    df_pts$n_pos <- ave(df_pts$analysis_chain_cl,
                        df_pts$analysis_chain_cl,
                        FUN = length)

    p <- ggplot2::ggplot(df_pts,
                         ggplot2::aes(x = analysis_chain_cl,
                                      y = .data[[fc_col]])) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60", linewidth = 0.5)

    # Individual points coloured by FDR significance
    if (show_points) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes(color = sig_label),
        size = 1.8, alpha = 0.65
      ) +
      ggplot2::scale_color_manual(
        values = c(
          setNames("#E63946", paste0("FDR < ", fdr_thresh)),
          NS = "grey55"
        ),
        name   = paste0("FDR < ", fdr_thresh),
        labels = c(setNames("FDR sig", paste0("FDR < ", fdr_thresh)),
                   NS = "NS")
      )
    }

    # Smoothing curve — weighted by n observations per position if requested
    smooth_args <- list(
      method    = smooth_method,
      se        = smooth_se,
      color     = "#2C3E50",
      fill      = "grey70",
      alpha     = 0.18,
      linewidth = 0.9
    )
    if (smooth_method == "loess")
      smooth_args$method.args <- list(span = smooth_span)

    if (smooth_weighted) {
      p <- p + do.call(ggplot2::geom_smooth,
                       c(list(ggplot2::aes(weight = n_pos)), smooth_args))
    } else {
      p <- p + do.call(ggplot2::geom_smooth, smooth_args)
    }

    # Statistical annotation
    ann_length <- .trend_annotation(
      x          = df_pts$analysis_chain_cl,
      y          = df_pts[[fc_col]],
      weights    = df_pts$n_pos,
      test       = trend_test,
      case_lbl   = case_lbl,
      ref_lbl    = ref_lbl,
      fdr_thresh = fdr_thresh
    )

    p <- p +
      ggplot2::labs(
        title    = paste0(cls, " \u2014 Chain Length Trend"),
        subtitle = paste0(
          case_lbl, " vs ", ref_lbl,
          "  |  n=", n_lip, " lipids / ", n_obs, " chain obs",
          "  |  method=", smooth_method,
          if (smooth_method == "loess") paste0("  span=", smooth_span) else "",
          "\n", .method_lbl(ctype),
          if (!is.null(ann_length)) paste0("\n", ann_length) else ""
        ),
        x = .axis_lbl(ctype, "cl"),
        y = paste0("logFC (", case_lbl, " / ", ref_lbl, ")")
      ) +
      .theme_chain() +
      ggplot2::theme(legend.position = "bottom") +
      {
        x_vals <- df_pts$analysis_chain_cl
        x_min  <- min(x_vals, na.rm = TRUE)
        x_max  <- max(x_vals, na.rm = TRUE)
        if (!is.null(trend_x_step_length) && is.numeric(trend_x_step_length)) {
          brks <- seq(floor(x_min / trend_x_step_length) * trend_x_step_length,
                      ceiling(x_max / trend_x_step_length) * trend_x_step_length,
                      by = trend_x_step_length)
          ggplot2::scale_x_continuous(
            breaks = brks,
            expand = ggplot2::expansion(mult = c(0.05, 0.05))
          )
        } else {
          ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.05))
          )
        }
      }

    plots[[paste0("trend_length_", cls)]] <- p
  }

  # -- Trend plots: unsaturation ------------------------------------------------
  for (cls in classes_trend) {
    df_cls <- df[df$LipidClass == cls &
                   !is.na(df$analysis_chain_cs), ]
    if (nrow(df_cls) < min_n_trend) next
    ctype  <- df_cls$chain_type[[1L]]
    n_lip  <- length(unique(df_cls[[1L]]))
    n_obs  <- nrow(df_cls)

    # Overfitting warning
    if (smooth_method == "loess" && smooth_span < 0.5 && n_obs < 10L)
      warning(cls, " trend_unsat: smooth_span=", smooth_span,
              " with only ", n_obs, " observations may overfit. ",
              "Consider increasing smooth_span.", call. = FALSE)

    df_pts <- df_cls
    if (has_sig) {
      df_pts$sig_label <- ifelse(df_pts[[sig_col]] == 1L,
                                 paste0("FDR < ", fdr_thresh), "NS")
    } else {
      df_pts$sig_label <- "NS"
    }

    df_pts$n_pos <- ave(df_pts$analysis_chain_cs,
                        df_pts$analysis_chain_cs,
                        FUN = length)

    p <- ggplot2::ggplot(df_pts,
                         ggplot2::aes(x = analysis_chain_cs,
                                      y = .data[[fc_col]])) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60", linewidth = 0.5)

    if (show_points) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes(color = sig_label),
        size = 1.8, alpha = 0.65
      ) +
      ggplot2::scale_color_manual(
        values = c(
          setNames("#E63946", paste0("FDR < ", fdr_thresh)),
          NS = "grey55"
        ),
        name   = paste0("FDR < ", fdr_thresh),
        labels = c(setNames("FDR sig", paste0("FDR < ", fdr_thresh)),
                   NS = "NS")
      )
    }

    smooth_args <- list(
      method    = smooth_method,
      se        = smooth_se,
      color     = "#2C3E50",
      fill      = "grey70",
      alpha     = 0.18,
      linewidth = 0.9
    )
    if (smooth_method == "loess")
      smooth_args$method.args <- list(span = smooth_span)

    if (smooth_weighted) {
      p <- p + do.call(ggplot2::geom_smooth,
                       c(list(ggplot2::aes(weight = n_pos)), smooth_args))
    } else {
      p <- p + do.call(ggplot2::geom_smooth, smooth_args)
    }

    # Statistical annotation
    ann_unsat <- .trend_annotation(
      x          = df_pts$analysis_chain_cs,
      y          = df_pts[[fc_col]],
      weights    = df_pts$n_pos,
      test       = trend_test,
      case_lbl   = case_lbl,
      ref_lbl    = ref_lbl,
      fdr_thresh = fdr_thresh
    )

    p <- p +
      ggplot2::labs(
        title    = paste0(cls, " \u2014 Unsaturation Trend"),
        subtitle = paste0(
          case_lbl, " vs ", ref_lbl,
          "  |  n=", n_lip, " lipids / ", n_obs, " chain obs",
          "  |  method=", smooth_method,
          if (smooth_method == "loess") paste0("  span=", smooth_span) else "",
          "\n", .method_lbl(ctype),
          if (!is.null(ann_unsat)) paste0("\n", ann_unsat) else ""
        ),
        x = .axis_lbl(ctype, "cs"),
        y = paste0("logFC (", case_lbl, " / ", ref_lbl, ")")
      ) +
      .theme_chain() +
      ggplot2::theme(legend.position = "bottom") +
      {
        x_vals <- df_pts$analysis_chain_cs
        x_min  <- min(x_vals, na.rm = TRUE)
        x_max  <- max(x_vals, na.rm = TRUE)
        if (!is.null(trend_x_step_unsat) && is.numeric(trend_x_step_unsat)) {
          brks <- seq(floor(x_min / trend_x_step_unsat) * trend_x_step_unsat,
                      ceiling(x_max / trend_x_step_unsat) * trend_x_step_unsat,
                      by = trend_x_step_unsat)
          ggplot2::scale_x_continuous(
            breaks = brks,
            expand = ggplot2::expansion(mult = c(0.05, 0.05))
          )
        } else {
          ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.05))
          )
        }
      }

    plots[[paste0("trend_unsat_", cls)]] <- p
  }

  message("Generated ", length(plots), " chain analysis plot(s).")
  plots
}
