# R/lsea.R

# Suppress R CMD check NOTEs for ggplot2 aes() column names used
# inside plot_distribution() via tidy evaluation.
utils::globalVariables(c(".data", "Convergence", "label"))
# Lipid Set Enrichment Analysis -- KS and fgsea engines.
# Source: LipidEnrichment_unified_v8.R
#
# Key changes from the original script:
#   1. Global variables (CASE_LBL, REF_LBL, N_SIG, N_TOTAL, FDR_THRESH,
#      MIN_N, N_PERM, FGSEA_EPS, FGSEA_NPERM, FGSEA_AVAILABLE, LIPID_ID_COL,
#      CONTRAST_TAG, OUT_DIR) -> function arguments with defaults.
#   2. set.seed(42) inside the function -> exposed as `seed` argument;
#      uses withr::with_seed() to avoid polluting the user's RNG state.
#   3. run_lsea()  (KS engine in script) -> .run_ks() internal.
#   4. run_lsea_fgsea() -> .run_fgsea_engine() internal.
#   5. New exported run_lsea() orchestrates both engines and returns a
#      single tidy data.frame.
#   6. map_dfr() / tibble() -> do.call(rbind, lapply()) / data.frame()
#      (avoids purrr and tibble as hard dependencies).
#   7. ORA kept as an internal helper; not in the primary run_lsea() output
#      by default (can be requested via include_ora = TRUE).
#   8. File I/O removed -- output is an R object, not files.

# -- Internal helpers ----------------------------------------------------------

#' Build lipid sets from a grouping column
#'
#' Converts a categorical column (e.g. LipidClass) to the named list of
#' character vectors required by fgsea.
#'
#' @noRd
.build_sets <- function(data, group_col, lipid_id_col) {
  grps <- unique(data[[group_col]])
  grps <- grps[!is.na(grps)]
  sets <- lapply(grps, function(g) {
    as.character(data[[lipid_id_col]][data[[group_col]] == g])
  })
  stats::setNames(sets, as.character(grps))
}

#' Detect lipid ID column
#'
#' Tries common column names; falls back to row index if none found.
#'
#' @noRd
.detect_lipid_id_col <- function(data) {
  candidates <- c("LipidName", "Lipid", "lipid", "Molecule", "molecule",
                   "feature", "Feature", "Name", "name")
  col <- intersect(candidates, names(data))[[1L]]
  if (length(col) == 0L || is.na(col)) {
    return(NULL)   # caller will create a synthetic ID
  }
  col
}

#' Ensure unique lipid IDs (required by fgsea)
#'
#' @noRd
.ensure_unique_ids <- function(data, id_col) {
  if (anyDuplicated(data[[id_col]])) {
    n_dup <- sum(duplicated(data[[id_col]]))
    warning(n_dup, " duplicate value(s) in lipid ID column '", id_col,
            "'. Making unique with make.unique() -- verify input data.",
            call. = FALSE)
    data[[id_col]] <- make.unique(as.character(data[[id_col]]))
  }
  data
}

#' Compute pi-value rank metric
#'
#' pi-value = sign(logFC) * -log10(P.Value raw).
#' Combines magnitude and statistical evidence for fgsea ranking.
#' Uses raw p-value (not adjusted) -- pi-value is a rank metric, not a
#' significance threshold. Reference: Xiao et al. (2014) Bioinformatics.
#'
#' @noRd
.compute_pi_value <- function(logfc, pval) {
  sign(logfc) * (-log10(pmax(pval, 1e-300)))
}

# -- KS engine -----------------------------------------------------------------

#' KS-based LSEA for one grouping level
#'
#' DirectionalScore: standardized mean difference (Cohen's d analog).
#' This is NOT the classic GSEA NES from the running sum.
#'
#' Inference hierarchy:
#'   PRIMARY  : KS p-value -> BH FDR (FDR_LSEA)
#'   SECONDARY: DS_perm_pval -- permutation p for DirectionalScore
#'              Evaluates whether the EFFECT SIZE exceeds permutation null.
#'              NOT a robustness check for the KS statistic.
#'              Min resolution = 1/n_perm.
#'
#' ContributingLipids_KS: lipids on the enriched side of the max CDF
#' divergence point (the KS statistic location). Analogous to fgsea
#' leading edge but derived from the KS test.
#'
#' @param data      data.frame with logFC column and a grouping column.
#' @param group_col Character(1). Name of the grouping column (e.g.
#'   "LipidClass").
#' @param fc_col    Character(1). Name of the logFC column.
#' @param lipid_id_col Character(1). Column with unique lipid IDs.
#' @param case_lbl  Character(1). Case group label for Direction column.
#' @param ref_lbl   Character(1). Reference group label.
#' @param min_n     Integer(1). Minimum set size to test.
#' @param n_perm    Integer(1). Permutations for DS_perm_pval.
#' @param seed      Integer(1) or NULL. RNG seed (via withr::with_seed).
#'
#' @return data.frame with one row per group: KS_pval, FDR_LSEA,
#'   DirectionalScore, DS_perm_pval, Direction, ContributingLipids_KS, etc.
#'
#' @noRd
.run_ks <- function(
    data,
    group_col,
    fc_col     = "logFC",
    lipid_id_col,
    case_lbl   = "Case",
    ref_lbl    = "Reference",
    min_n      = 3L,
    n_perm     = 2000L,
    seed       = 42L
) {

  n_total <- nrow(data)

  # Groups with sufficient n
  grp_sizes <- tapply(seq_len(n_total), data[[group_col]], length)
  groups    <- names(grp_sizes)[grp_sizes >= min_n]

  if (length(groups) == 0L) {
    warning("No group has >= ", min_n, " lipids. Returning empty KS results.",
            call. = FALSE)
    return(data.frame())
  }

  run_one <- function(grp) {
    grp_idx <- data[[group_col]] == grp
    grp_fc  <- data[[fc_col]][grp_idx]
    bg_fc   <- data[[fc_col]][!grp_idx]

    # Guard: ks.test() requires at least 1 background observation
    if (length(bg_fc) == 0L) {
      warning("Group '", grp, "': background is empty (no other lipids). ",
              "Ensure data contains more than one lipid class.", call. = FALSE)
      return(NULL)
    }

    # KS test (two-sided)
    ks_p <- stats::ks.test(grp_fc, bg_fc, alternative = "two.sided")$p.value

    # DirectionalScore (Cohen's d analog)
    pool <- c(grp_fc, bg_fc)
    ps   <- stats::sd(pool)
    ds   <- if (ps > 0) (mean(grp_fc) - mean(bg_fc)) / ps else 0

    # DS_perm_pval: permutation p for the DirectionalScore (NOT for the KS)
    obs    <- abs(ds)
    perm_v <- replicate(n_perm, {
      p   <- sample(pool)
      ps2 <- stats::sd(p)
      if (ps2 > 0)
        abs((mean(p[seq_along(grp_fc)]) - mean(p[-seq_along(grp_fc)])) / ps2)
      else 0
    })
    ds_perm_p <- mean(perm_v >= obs)

    # ContributingLipids_KS: lipids on enriched side of max CDF divergence
    contrib_str <- tryCatch({
      all_vals  <- sort(unique(c(grp_fc, bg_fc)))
      ecdf_grp  <- stats::ecdf(grp_fc)
      ecdf_bg   <- stats::ecdf(bg_fc)
      cdf_diff  <- ecdf_grp(all_vals) - ecdf_bg(all_vals)
      ks_thresh <- all_vals[which.max(abs(cdf_diff))]

      lipid_ids <- as.character(data[[lipid_id_col]][grp_idx])
      contrib_ids <- if (ds >= 0) {
        lipid_ids[grp_fc >= ks_thresh]
      } else {
        lipid_ids[grp_fc <= ks_thresh]
      }
      if (length(contrib_ids) == 0L) "--"
      else paste(contrib_ids, collapse = "; ")
    }, error = function(e) "--")

    n_contrib <- if (identical(contrib_str, "--")) 0L
                 else length(strsplit(contrib_str, "; ")[[1L]])

    data.frame(
      Group                 = grp,
      N_group               = length(grp_fc),
      Mean_logFC            = mean(grp_fc),
      DirectionalScore      = ds,
      KS_pval               = ks_p,
      DS_perm_pval          = ds_perm_p,
      Direction             = ifelse(ds > 0,
                                     paste0("up in ", case_lbl),
                                     paste0("up in ", ref_lbl)),
      ContributingLipids_KS = contrib_str,
      N_contributing_KS     = n_contrib,
      stringsAsFactors      = FALSE
    )
  }

  rows <- if (!is.null(seed)) {
    withr::with_seed(seed, lapply(groups, run_one))
  } else {
    lapply(groups, run_one)
  }

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame())

  out <- do.call(rbind, rows)
  out$FDR_LSEA <- stats::p.adjust(out$KS_pval, method = "BH")
  out$Level    <- group_col
  out[order(out$KS_pval), ]
}

# -- fgsea engine --------------------------------------------------------------

#' fgsea enrichment for one grouping level
#'
#' Uses the fast preranked GSEA algorithm (Korotkevich et al.).
#' Key differences vs KS:
#'   KS    -> any distributional shift (location, shape, spread)
#'   fgsea -> extreme-rank clustering (members at top/bottom of ranking)
#'
#' Rank metrics:
#'   "pi_value" (default): sign(logFC) * -log10(P.Value raw).
#'   "logFC": magnitude only (same as lipidr default).
#'   "t_stat": LIMMA t-statistic (gold standard for GSEA in genomics).
#'
#' @param data        data.frame with rank metric column and grouping column.
#' @param group_col   Character(1). Grouping column name.
#' @param lipid_id_col Character(1). Column with unique lipid IDs.
#' @param rank_col    Character(1). Name of the rank metric column.
#' @param min_n       Integer(1). Minimum set size.
#' @param eps         Numeric(1). fgsea eps parameter (0 = reduce p-value
#'   approximation error).
#' @param nperm       Integer(1). Monte Carlo permutations.
#' @param seed        Integer(1) or NULL.
#'
#' @return data.frame with NES, fgsea pval/padj, leading edge, or NULL if
#'   fgsea is not installed or rank metric is unavailable.
#'
#' @noRd
.run_fgsea_engine <- function(
    data,
    group_col,
    lipid_id_col,
    rank_col   = "pi_value",
    min_n      = 3L,
    eps        = 0,
    nperm      = 10000L,
    seed       = 42L
) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    message("fgsea not installed -- skipping fgsea engine.\n",
            "Install with: BiocManager::install('fgsea')")
    return(NULL)
  }

  if (!rank_col %in% names(data)) {
    message("Rank column '", rank_col, "' not found -- skipping fgsea engine.")
    return(NULL)
  }

  # Build named rank vector; break ties with tiny jitter
  rv <- data[[rank_col]]
  if (anyNA(rv)) {
    n_na <- sum(is.na(rv))
    warning(n_na, " NA(s) in rank column '", rank_col,
            "' removed before fgsea.", call. = FALSE)
    keep <- !is.na(rv)
    data <- data[keep, ]
    rv   <- rv[keep]
  }

  # Jitter to break ties (fgsea requires no ties)
  set.seed_local <- if (!is.null(seed)) seed else sample.int(.Machine$integer.max, 1L)
  rv <- rv + withr::with_seed(
    set.seed_local,
    stats::rnorm(length(rv), sd = 1e-6)
  )
  names(rv) <- as.character(data[[lipid_id_col]])

  # Build lipid sets
  sets <- .build_sets(data, group_col, lipid_id_col)
  sets <- sets[vapply(sets, length, integer(1L)) >= min_n]

  if (length(sets) == 0L) {
    message("No set has >= ", min_n, " members -- skipping fgsea.")
    return(NULL)
  }

  fg_res <- withr::with_seed(set.seed_local, {
    fgsea::fgsea(
      pathways  = sets,
      stats     = rv,
      minSize   = min_n,
      maxSize   = 10000L,
      eps       = eps,
      nPermSimple = nperm
    )
  })

  if (is.null(fg_res) || nrow(fg_res) == 0L) return(NULL)

  # Tidy output: one row per set
  out <- data.frame(
    Group         = as.character(fg_res$pathway),
    NES           = fg_res$NES,
    fgsea_pval    = fg_res$pval,
    FDR_fgsea     = fg_res$padj,
    N_leading     = fg_res$size,   # size after minSize/maxSize filtering
    LeadingEdge   = vapply(fg_res$leadingEdge,
                            paste, character(1L), collapse = "; "),
    rank_metric   = rank_col,
    Level         = group_col,
    stringsAsFactors = FALSE
  )
  out[order(out$fgsea_pval), ]
}

# -- Interpretation helper -----------------------------------------------------

#' Generate plain-text interpretation for one lipid set
#'
#' Mirrors .interpret() from LipidEnrichment_unified_v8.R.
#' DS_ZERO_THRESH: |DirectionalScore| below this value = no net shift.
#'
#' @noRd
.interpret_set <- function(
    ora_sig, ks_sig, fgsea_sig,
    fgsea_available,
    direction, ds, fold_enr,
    case_lbl, ref_lbl,
    ds_zero_thresh = 0.10
) {

  ds_zero    <- abs(ds) < ds_zero_thresh
  up_in_case <- grepl(paste0("up in ", case_lbl), direction, fixed = TRUE)

  if (!ks_sig && !fgsea_sig && !ora_sig) {
    return("No significant enrichment detected by any method.")
  }

  if (ks_sig && fgsea_sig && ora_sig) {
    dir_str <- if (up_in_case)
      paste0("enriched in ", case_lbl)
    else
      paste0("depleted in ", case_lbl, " (enriched in ", ref_lbl, ")")
    return(paste0(
      "Convergent signal across all three methods (ORA, KS, fgsea): ",
      dir_str, ". The lipid set shows both over-representation among ",
      "significant features (ORA Fold=", round(fold_enr, 1L), ") and a ",
      "coordinated directional shift of logFC values (KS + fgsea). ",
      "High confidence."
    ))
  }

  if (ks_sig && !fgsea_sig) {
    if (ds_zero) {
      return(paste0(
        "KS significant but DirectionalScore \u2248 0. ",
        "The set shows a distributional difference from background ",
        "(shape or variance shift) without a net directional effect. ",
        "Interpret with caution -- this pattern suggests heterogeneous ",
        "within-set logFC rather than coordinated remodeling."
      ))
    }
    dir_str <- if (up_in_case) case_lbl else ref_lbl
    return(paste0(
      "KS-significant (distributional shift toward ", dir_str, ") ",
      "without fgsea support. Consistent with a distributed moderate effect ",
      "across many set members rather than extreme-rank clustering. ",
      "Consider reporting DS=", round(ds, 2L), " alongside the KS FDR."
    ))
  }

  if (!ks_sig && fgsea_sig) {
    dir_str <- if (up_in_case) case_lbl else ref_lbl
    return(paste0(
      "fgsea-significant (extreme-rank enrichment toward ", dir_str, ") ",
      "without KS support. A few strongly regulated lipids drive the signal; ",
      "the overall set distribution does not differ from background. ",
      "Examine the leading edge for the key lipids."
    ))
  }

  if (ora_sig && !ks_sig && !fgsea_sig) {
    return(paste0(
      "Over-represented among significant lipids (ORA Fold=",
      round(fold_enr, 1L), ") but no directional enrichment detected ",
      "by KS or fgsea. The significant lipids in this class are over-",
      "represented without a consistent direction of change."
    ))
  }

  "Partial enrichment signal -- see individual method results."
}

# -- Main exported function ----------------------------------------------------

#' Lipid Set Enrichment Analysis
#'
#' Runs KS-based LSEA, fgsea, or both for each grouping level in
#' \code{group_cols} and returns a tidy \code{data.frame} with enrichment
#' statistics.
#'
#' @param data A \code{data.frame} with at minimum a lipid identifier column,
#'   a log2 fold-change column, and one or more grouping columns (e.g.
#'   \code{LipidClass}).
#' @param group_cols Character vector. Names of grouping columns to test.
#'   Each column defines one level of analysis (e.g. class, LIPID MAPS
#'   category, functional category). Default:
#'   \code{c("LipidClass", "LipidCategory_LMAPS", "LipidCategory_functional")}.
#' @param fc_col Character(1). Log2 fold-change column. Default:
#'   \code{"logFC"}.
#' @param pval_col Character(1) or \code{NULL}. Raw p-value column used to
#'   compute the pi-value rank metric. If \code{NULL}, logFC is used as the
#'   fgsea rank metric. Default: \code{"P.Value"}.
#' @param lipid_id_col Character(1) or \code{NULL}. Column with unique lipid
#'   identifiers. If \code{NULL}, auto-detected from common column names.
#' @param case_lbl Character(1). Case group label. Default: \code{"Case"}.
#' @param ref_lbl Character(1). Reference group label. Default:
#'   \code{"Reference"}.
#' @param engine Character(1). Enrichment engine: \code{"ks"},
#'   \code{"fgsea"}, or \code{"both"}. Default: \code{"both"}.
#' @param fgsea_rank Character(1). Rank metric for fgsea: \code{"pi_value"},
#'   \code{"logFC"}, or \code{"t_stat"}. Default: \code{"pi_value"}.
#' @param min_n Integer(1). Minimum set size to test. Default: \code{3L}.
#' @param n_perm Integer(1). KS permutations for \code{DS_perm_pval}.
#'   Default: \code{2000L}.
#' @param fgsea_nperm Integer(1). fgsea Monte Carlo permutations.
#'   Default: \code{10000L}.
#' @param fgsea_eps Numeric(1). fgsea epsilon (0 = reduce approximation
#'   error). Default: \code{0}.
#' @param seed Integer(1) or \code{NULL}. RNG seed passed to
#'   \code{withr::with_seed()} -- does not alter the user's global RNG state.
#'   Default: \code{42L}.
#' @param verbose Logical(1). Print progress messages. Default: \code{TRUE}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{ks}}{data.frame of KS results (or \code{NULL} if
#'       \code{engine = "fgsea"}).}
#'     \item{\code{fgsea}}{data.frame of fgsea results (or \code{NULL} if
#'       \code{engine = "ks"} or fgsea is not installed).}
#'     \item{\code{combined}}{data.frame merging both engines by Group and
#'       Level, including a Convergence column.}
#'   }
#'
#' @references
#' Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A
#' (2021). Fast gene set enrichment analysis. \emph{bioRxiv}.
#' \doi{10.1101/060012}
#'
#' Xiao Y, Hsiao TH, Suresh U, Chen HI, Wu X, Wolf SE, Chen Y (2014).
#' A novel significance score for gene selection and ranking.
#' \emph{Bioinformatics}, 30(6), 801--807. \doi{10.1093/bioinformatics/btr671}
#'
#' @seealso \code{annotate_lipids()}, \code{\link{plot_lsea}},
#'   \code{export_lsea()}
#'
#' @examples
#' \dontrun{
#' data("lipid_example", package = "easyLSEA")
#' annotated <- annotate_lipids(lipid_example)
#'
#' result <- run_lsea(
#'   data      = annotated,
#'   fc_col    = "logFC",
#'   engine    = "ks",
#'   case_lbl  = "NASH",
#'   ref_lbl   = "Control"
#' )
#'
#' head(result$ks)
#' }
#'
#' @importFrom stats ks.test p.adjust ecdf sd weighted.mean setNames
#' @importFrom withr with_seed
#'
#' @export
run_lsea <- function(
    data,
    group_cols    = c("LipidClass",
                      "LipidCategory_LMAPS",
                      "LipidCategory_functional"),
    fc_col        = "logFC",
    pval_col      = "P.Value",
    lipid_id_col  = NULL,
    case_lbl      = "Case",
    ref_lbl       = "Reference",
    engine        = c("both", "ks", "fgsea"),
    fgsea_rank    = c("pi_value", "logFC", "t_stat"),
    min_n         = 3L,
    n_perm        = 2000L,
    fgsea_nperm   = 10000L,
    fgsea_eps     = 0,
    seed          = 42L,
    verbose       = TRUE
) {

  engine     <- match.arg(engine)
  fgsea_rank <- match.arg(fgsea_rank)

  # -- Input validation --------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)

  missing_gcols <- setdiff(group_cols, names(data))
  if (length(missing_gcols) > 0L)
    stop("group_col(s) not found: ",
         paste(sQuote(missing_gcols), collapse = ", "), call. = FALSE)

  if (!fc_col %in% names(data))
    stop("fc_col ", sQuote(fc_col), " not found in 'data'.", call. = FALSE)

  # -- Lipid ID column ---------------------------------------------------------
  if (is.null(lipid_id_col)) {
    lipid_id_col <- .detect_lipid_id_col(data)
    if (is.null(lipid_id_col)) {
      data$.lipid_id <- paste0("lipid_", seq_len(nrow(data)))
      lipid_id_col   <- ".lipid_id"
      if (verbose)
        message("No lipid ID column found -- using row index as lipid ID.")
    }
  }
  data <- .ensure_unique_ids(data, lipid_id_col)

  # -- Remove NA in fc_col -----------------------------------------------------
  na_fc <- is.na(data[[fc_col]])
  if (any(na_fc)) {
    warning(sum(na_fc), " NA value(s) in '", fc_col, "' removed.",
            call. = FALSE)
    data <- data[!na_fc, ]
  }

  # -- Compute pi-value if available ------------------------------------------
  if (fgsea_rank == "pi_value") {
    if (!is.null(pval_col) && pval_col %in% names(data)) {
      data$pi_value <- .compute_pi_value(data[[fc_col]], data[[pval_col]])
      if (verbose)
        message("pi-value computed from '", fc_col, "' and '", pval_col, "'.")
    } else {
      if (verbose)
        message("'", pval_col, "' not found -- using logFC as fgsea rank metric.")
      data$pi_value <- data[[fc_col]]
    }
  }

  # -- Run engines across all group levels -------------------------------------
  ks_results    <- list()
  fgsea_results <- list()

  for (gcol in group_cols) {

    if (verbose) message("Level: ", gcol, " ...")

    if (engine %in% c("ks", "both")) {
      ks_results[[gcol]] <- .run_ks(
        data         = data,
        group_col    = gcol,
        fc_col       = fc_col,
        lipid_id_col = lipid_id_col,
        case_lbl     = case_lbl,
        ref_lbl      = ref_lbl,
        min_n        = min_n,
        n_perm       = n_perm,
        seed         = seed
      )
    }

    if (engine %in% c("fgsea", "both")) {
      rank_col_use <- switch(fgsea_rank,
        pi_value = "pi_value",
        logFC    = fc_col,
        t_stat   = if ("t" %in% names(data)) "t"
                   else if ("t_stat" %in% names(data)) "t_stat"
                   else { message("t-statistic column not found."); fc_col }
      )
      fgsea_results[[gcol]] <- .run_fgsea_engine(
        data         = data,
        group_col    = gcol,
        lipid_id_col = lipid_id_col,
        rank_col     = rank_col_use,
        min_n        = min_n,
        eps          = fgsea_eps,
        nperm        = fgsea_nperm,
        seed         = seed
      )
    }
  }

  # -- Bind results ------------------------------------------------------------
  ks_df    <- if (length(ks_results) > 0L)
                do.call(rbind, Filter(Negate(is.null), ks_results))
              else NULL

  fgsea_df <- if (length(fgsea_results) > 0L)
                do.call(rbind, Filter(Negate(is.null), fgsea_results))
              else NULL

  # -- Combined table (KS + fgsea + Convergence) ----------------------------
  combined <- NULL
  if (!is.null(ks_df) && !is.null(fgsea_df)) {
    combined <- merge(
      ks_df,
      fgsea_df[, c("Group", "Level", "NES", "FDR_fgsea",
                    "fgsea_pval", "LeadingEdge", "N_leading")],
      by      = c("Group", "Level"),
      all     = TRUE
    )
    # Convergence: both KS and fgsea significant at FDR < 0.05
    fdr_thresh <- 0.05   # default; could be exposed as argument in v2
    combined$Convergence <- ifelse(
      !is.na(combined$FDR_LSEA) & combined$FDR_LSEA < fdr_thresh &
        !is.na(combined$FDR_fgsea) & combined$FDR_fgsea < fdr_thresh,
      "KS+fgsea",
      ifelse(!is.na(combined$FDR_LSEA) & combined$FDR_LSEA < fdr_thresh,
             "KS only",
             ifelse(!is.na(combined$FDR_fgsea) & combined$FDR_fgsea < fdr_thresh,
                    "fgsea only", "NS"))
    )
  } else if (!is.null(ks_df)) {
    combined <- ks_df
  } else if (!is.null(fgsea_df)) {
    combined <- fgsea_df
  }

  if (verbose) {
    n_ks_sig    <- if (!is.null(ks_df))
                     sum(ks_df$FDR_LSEA < 0.05, na.rm = TRUE) else 0L
    n_fgsea_sig <- if (!is.null(fgsea_df))
                     sum(fgsea_df$FDR_fgsea < 0.05, na.rm = TRUE) else 0L
    message("Done. KS sig sets (FDR<0.05): ", n_ks_sig,
            " | fgsea sig sets: ", n_fgsea_sig)
  }

  list(
    ks       = ks_df,
    fgsea    = fgsea_df,
    combined = combined
  )
}

# -- Plot stubs ----------------------------------------------------------------

#' Generate LSEA enrichment plots
#'
#' Produces bubble, barplot, and running sum plots from a \code{run_lsea()}
#' result. Returns a named list of \code{\link[ggplot2]{ggplot}} objects.
#'
#' @param lsea_result Named list returned by \code{\link{run_lsea}}.
#' @param which Character vector. Which plots to generate:
#'   \code{"bubble_ks"}, \code{"bubble_fgsea"}, \code{"bubble_combined"},
#'   \code{"barplot"}, \code{"running_sum"}. Default: all.
#' @param fdr_thresh Numeric(1). Significance threshold for highlighting.
#'   Default: \code{0.05}.
#' @param case_lbl Character(1). Case label for plot annotations.
#' @param ref_lbl Character(1). Reference label for plot annotations.
#' @param bubble_label Character vector. Which statistics to display next to
#'   each bubble. Any subset of \code{"FDR"}, \code{"DS"} (KS plots only),
#'   \code{"NES"} (fgsea plots only), and \code{"n"}. Default: all four.
#'
#' @return Named list of \code{ggplot} objects.
#'
#' @seealso \code{\link{run_lsea}}, \code{export_lsea()}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_col geom_hline
#' @importFrom ggplot2 scale_color_manual scale_size_continuous labs theme_bw theme
#' @importFrom ggplot2 element_text element_blank facet_wrap coord_flip
#'
#' @export
plot_lsea <- function(
    lsea_result,
    which        = c("bubble_ks", "bubble_fgsea", "bubble_combined",
                     "barplot", "running_sum"),
    fdr_thresh   = 0.05,
    case_lbl     = "Case",
    ref_lbl      = "Reference",
    bubble_label = c("FDR", "DS", "NES", "n")
) {

  bubble_label <- match.arg(bubble_label, several.ok = TRUE)

  # Helper: build the per-bubble label string from selected components.
  # KS bubbles can show FDR, DS, n; fgsea bubbles can show FDR, NES, n.
  .bubble_lab <- function(parts, fdr = NULL, ds = NULL, nes = NULL, n = NULL) {
    bits <- character(0)
    if ("FDR" %in% parts && !is.null(fdr))
      bits <- c(bits, paste0("FDR=", formatC(fdr, format = "e", digits = 1)))
    if ("DS"  %in% parts && !is.null(ds))
      bits <- c(bits, paste0("DS=", sprintf("%+.2f", ds)))
    if ("NES" %in% parts && !is.null(nes))
      bits <- c(bits, paste0("NES=", sprintf("%+.2f", nes)))
    if ("n"   %in% parts && !is.null(n))
      bits <- c(bits, paste0("n=", n))
    paste(bits, collapse = "   ")
  }

  plots <- list()

  # -- KS bubble plot — one plot per analysis level ----------------------------
  if ("bubble_ks" %in% which && !is.null(lsea_result$ks)) {
    df_all <- lsea_result$ks
    df_all <- df_all[!is.na(df_all$KS_pval), ]
    df_all$sig       <- df_all$FDR_LSEA < fdr_thresh
    df_all$log_fdr   <- -log10(df_all$FDR_LSEA + 1e-300)
    df_all$direction <- ifelse(df_all$DirectionalScore > 0,
                               paste0("Up in ", case_lbl),
                               paste0("Up in ", ref_lbl))

    levels_ks <- unique(df_all$Level)
    level_labels <- c(
      LipidClass                 = "01_Class",
      LipidCategory_LMAPS        = "02_LMAPS",
      LipidCategory_functional   = "03_Functional"
    )

    for (lvl in levels_ks) {
      df  <- df_all[df_all$Level == lvl, ]
      lbl <- if (lvl %in% names(level_labels)) level_labels[[lvl]] else lvl

      .make_bubble_ks <- function(dat, key_name, sig_only = FALSE) {
        if (nrow(dat) == 0L) return(NULL)
        dat <- dat[order(dat$DirectionalScore), ]
        dat$Group <- factor(dat$Group, levels = dat$Group)
        max_score <- max(abs(dat$DirectionalScore), na.rm = TRUE)
        x_range   <- diff(range(dat$DirectionalScore, na.rm = TRUE))

        # Dynamic axis expansion:
        #   Right — estimate label width in data units.
        #     Label format: "FDR=X.Xe-XX   DS=+X.XX   n=XXX"
        #     ~30 chars at text size 3 (≈2.1 mm/char) = ~63 mm
        #     Convert to data units: 63mm / (plot_width_mm / x_range)
        #     We approximate plot_width as 120mm (safe for 8" wide PDF minus legend).
        #   Left — accommodate the largest bubble radius bleeding past its center.
        #     Bubble size range = c(3,12) pt mapped to N_group.
        #     Largest bubble radius ≈ 12pt ≈ 4mm → in data units: 4/120 * x_range.
        label_chars   <- 32L          # conservative char count for KS label
        char_mm       <- 2.1          # mm per char at geom_text size=3
        plot_width_mm <- 120          # approx panel width in mm
        label_data_units <- (label_chars * char_mm / plot_width_mm) * x_range
        nudge_data_units <- max_score * 0.06
        right_add <- label_data_units + nudge_data_units

        bubble_radius_mm  <- 4        # approx max bubble radius in mm
        left_add  <- max(x_range * 0.15, 0.20)

        p <- ggplot2::ggplot(
          dat,
          ggplot2::aes(x = DirectionalScore, y = Group,
                       size = N_group, color = DirectionalScore)
        ) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                               color = "grey60") +
          ggplot2::geom_point(
            ggplot2::aes(alpha = FDR_LSEA < fdr_thresh),
            show.legend = TRUE
          ) +
          ggplot2::geom_text(
            ggplot2::aes(
              label = mapply(.bubble_lab,
                             MoreArgs = list(parts = bubble_label),
                             fdr = FDR_LSEA, ds = DirectionalScore, n = N_group)
            ),
            nudge_x = nudge_data_units, hjust = 0,
            size = 3.0, fontface = "bold", color = "grey25",
            show.legend = FALSE
          ) +
          ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(add = c(left_add, right_add))
          ) +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::scale_color_gradient2(
            low = "#2166AC", mid = "grey85", high = "#E63946",
            midpoint = 0, name = "Directional\nScore",
            limits = c(-max_score, max_score)
          ) +
          ggplot2::scale_size_continuous(range = c(3, 12), name = "N lipids") +
          ggplot2::scale_alpha_manual(
            values = c("TRUE" = 0.9, "FALSE" = 0.35),
            labels = c("TRUE"  = paste0("FDR < ", fdr_thresh),
                       "FALSE" = paste0("FDR >= ", fdr_thresh)),
            name   = "KS sig."
          ) +
          ggplot2::labs(
            title    = paste0("LSEA (KS) \u00b7 ", key_name,
                              if (sig_only) "  [significant sets only]" else ""),
            subtitle = paste0(case_lbl, " vs ", ref_lbl,
                              " \u00b7 Kolmogorov-Smirnov test \u00b7 BH FDR",
                              if (sig_only) paste0(" \u00b7 FDR < ", fdr_thresh, " only")
                              else paste0(" \u00b7 faded = FDR >= ", fdr_thresh)),
            x = "Directional Score (standardized mean difference)", y = NULL
          ) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(
            plot.title       = ggplot2::element_text(face = "bold", size = 12,
                                                      hjust = 0.5),
            plot.subtitle    = ggplot2::element_text(size = 9, hjust = 0.5,
                                                      color = "grey50"),
            axis.text.y      = ggplot2::element_text(size = 10, face = "bold"),
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin      = ggplot2::margin(t = 8, r = 20, b = 8, l = 20)
          )
        attr(p, "n_sets")   <- nrow(dat)
        attr(p, "sig_only") <- sig_only
        p
      }

      plots[[paste0("bubble_ks_", lbl)]] <- .make_bubble_ks(df, lbl)
      df_sig <- df[df$FDR_LSEA < fdr_thresh, ]
      if (nrow(df_sig) > 0L)
        plots[[paste0("bubble_ks_sig_", lbl)]] <- .make_bubble_ks(
          df_sig, lbl, sig_only = TRUE
        )
    }
  }

  # -- fgsea bubble plot — one plot per analysis level -------------------------
  if ("bubble_fgsea" %in% which && !is.null(lsea_result$fgsea)) {
    df_all <- lsea_result$fgsea
    df_all <- df_all[!is.na(df_all$NES), ]
    df_all$sig       <- df_all$FDR_fgsea < fdr_thresh
    df_all$direction <- ifelse(df_all$NES > 0,
                               paste0("Up in ", case_lbl),
                               paste0("Up in ", ref_lbl))

    levels_fgsea <- unique(df_all$Level)
    level_labels <- c(
      LipidClass                 = "01_Class",
      LipidCategory_LMAPS        = "02_LMAPS",
      LipidCategory_functional   = "03_Functional"
    )

    for (lvl in levels_fgsea) {
      df  <- df_all[df_all$Level == lvl, ]
      lbl <- if (lvl %in% names(level_labels)) level_labels[[lvl]] else lvl

      .make_bubble_fgsea <- function(dat, key_name, sig_only = FALSE) {
        dat <- dat[!is.na(dat$NES), ]
        if (nrow(dat) == 0L) return(NULL)
        dat <- dat[order(dat$NES), ]
        dat$Group <- factor(dat$Group, levels = dat$Group)
        max_nes <- max(abs(dat$NES), na.rm = TRUE)
        x_range <- diff(range(dat$NES, na.rm = TRUE))
        # When all NES values are identical (e.g. single-point sig plot),
        # x_range = 0 which would zero out the expansion. Use max_nes as proxy.
        if (x_range < 1e-6) x_range <- max_nes

        # Dynamic axis expansion (same logic as .make_bubble_ks):
        #   fgsea label: "FDR=X.Xe-XX   NES=+X.XX   n=XXX" ~ 32 chars
        label_chars   <- 32L
        char_mm       <- 2.1
        plot_width_mm <- 120
        label_data_units <- (label_chars * char_mm / plot_width_mm) * x_range
        nudge_data_units <- max_nes * 0.06
        right_add <- label_data_units + nudge_data_units

        bubble_radius_mm <- 4
        left_add <- max(x_range * 0.15, 0.20)

        p <- ggplot2::ggplot(
          dat,
          ggplot2::aes(x = NES, y = Group, size = N_leading, color = NES)
        ) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                               color = "grey60") +
          ggplot2::geom_point(
            ggplot2::aes(alpha = FDR_fgsea < fdr_thresh),
            show.legend = TRUE
          ) +
          ggplot2::geom_text(
            ggplot2::aes(
              label = mapply(.bubble_lab,
                             MoreArgs = list(parts = bubble_label),
                             fdr = FDR_fgsea, nes = NES, n = N_leading)
            ),
            nudge_x = nudge_data_units, hjust = 0,
            size = 3.0, fontface = "bold", color = "grey25",
            show.legend = FALSE
          ) +
          ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(add = c(left_add, right_add))
          ) +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::scale_color_gradient2(
            low = "#2166AC", mid = "grey85", high = "#E63946",
            midpoint = 0, name = "NES",
            limits = c(-max_nes, max_nes)
          ) +
          ggplot2::scale_size_continuous(range = c(3, 12), name = "N lipids") +
          ggplot2::scale_alpha_manual(
            values = c("TRUE" = 0.9, "FALSE" = 0.35),
            labels = c("TRUE"  = paste0("FDR < ", fdr_thresh),
                       "FALSE" = paste0("FDR >= ", fdr_thresh)),
            name   = "fgsea sig."
          ) +
          ggplot2::labs(
            title    = paste0("LSEA (fgsea) \u00b7 ", key_name,
                              if (sig_only) "  [significant sets only]" else ""),
            subtitle = paste0(case_lbl, " vs ", ref_lbl, " \u00b7 fgsea \u00b7 BH FDR",
                              if (sig_only) paste0(" \u00b7 FDR < ", fdr_thresh, " only")
                              else paste0(" \u00b7 faded = FDR >= ", fdr_thresh)),
            x = "NES (Normalized Enrichment Score)", y = NULL
          ) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(
            plot.title       = ggplot2::element_text(face = "bold", size = 12,
                                                      hjust = 0.5),
            plot.subtitle    = ggplot2::element_text(size = 9, hjust = 0.5,
                                                      color = "grey50"),
            axis.text.y      = ggplot2::element_text(size = 10, face = "bold"),
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin      = ggplot2::margin(t = 8, r = 20, b = 8, l = 20)
          )
        attr(p, "n_sets")   <- nrow(dat)
        attr(p, "sig_only") <- sig_only
        p
      }

      plots[[paste0("bubble_fgsea_", lbl)]] <- .make_bubble_fgsea(df, lbl)
      df_sig <- df[df$FDR_fgsea < fdr_thresh, ]
      if (nrow(df_sig) > 0L)
        plots[[paste0("bubble_fgsea_sig_", lbl)]] <- .make_bubble_fgsea(
          df_sig, lbl, sig_only = TRUE
        )
    }
  }

  if (length(plots) == 0L)
    message("No plots generated -- check 'which' and that results are not NULL.")

  plots
}


#' Distribution enrichment boxplot per lipid set
#'
#' Produces a boxplot of logFC distributions for each lipid set, with jittered
#' individual lipid points, FDR/DS/NES labels for significant sets, and
#' red borders for significant sets. When \code{engine = "both"} (KS + fgsea),
#' fill colour encodes convergence (KS only, fgsea only, or KS+fgsea).
#'
#' @param data A \code{data.frame} with at least \code{fc_col} and the grouping
#'   column (e.g. \code{LipidClass}).
#' @param lsea_result A named list as returned by \code{\link{run_lsea}},
#'   with elements \code{ks}, \code{fgsea}, and/or \code{combined}.
#' @param group_col Character(1). Grouping column name
#'   (e.g. \code{"LipidClass"}).
#' @param fc_col Character(1). Column with log fold-change values.
#'   Default: \code{"logFC"}.
#' @param case_lbl Character(1). Label for the case group. Default:
#'   \code{"Case"}.
#' @param ref_lbl Character(1). Label for the reference group. Default:
#'   \code{"Control"}.
#' @param fdr_thresh Numeric(1). FDR threshold for significance.
#'   Default: \code{0.05}.
#' @param min_n Integer(1). Minimum number of lipids per set to include.
#'   Default: \code{3L}.
#' @param sig_only Logical(1). If \code{TRUE}, show only significant sets.
#'   Default: \code{FALSE}.
#' @param label_angle Numeric(1). Angle for FDR labels. \code{0} = horizontal
#'   (default); \code{90} = vertical (useful when many groups).
#'
#' @return A \code{ggplot} object, or \code{NULL} if no groups pass
#'   \code{min_n}.
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_boxplot geom_jitter geom_label scale_color_manual coord_cartesian labs theme_bw theme element_text element_blank unit
#' @importFrom stats median IQR quantile
#' @importFrom grDevices adjustcolor
#'
#' @export
plot_distribution <- function(
    data,
    lsea_result,
    group_col,
    fc_col      = "logFC",
    case_lbl    = "Case",
    ref_lbl     = "Control",
    fdr_thresh  = 0.05,
    min_n       = 3L,
    sig_only    = FALSE,
    label_angle = 0
) {

  # -- Detect available engines -----------------------------------------------
  has_ks    <- !is.null(lsea_result$ks)
  has_fgsea <- !is.null(lsea_result$fgsea)
  has_both  <- has_ks && has_fgsea && !is.null(lsea_result$combined)

  # -- Convergence palette ----------------------------------------------------
  pal_conv <- c(
    "KS+fgsea [strongest]"         = "#E63946",
    "KS only [distributed effect]" = "#0096C7",
    "fgsea only [extreme-driven]"  = "#7B2D8B",
    "Neither"                      = "grey55"
  )

  # -- Build per-group metadata -----------------------------------------------
  if (has_both) {
    sig_meta <- lsea_result$combined
    sig_meta$.is_sig_any <- (!is.na(sig_meta$FDR_LSEA)  & sig_meta$FDR_LSEA  < fdr_thresh) |
                            (!is.na(sig_meta$FDR_fgsea) & sig_meta$FDR_fgsea < fdr_thresh)
    # Map combined Convergence short labels to palette labels
    conv_map <- c(
      "KS+fgsea"   = "KS+fgsea [strongest]",
      "KS only"    = "KS only [distributed effect]",
      "fgsea only" = "fgsea only [extreme-driven]",
      "NS"         = "Neither",
      "Neither"    = "Neither"
    )
    sig_meta$Convergence <- ifelse(
      sig_meta$Convergence %in% names(conv_map),
      conv_map[sig_meta$Convergence],
      "Neither"
    )

  } else if (has_ks) {
    sig_meta <- lsea_result$ks[lsea_result$ks$Level == group_col, ]
    sig_meta$Convergence  <- ifelse(
      !is.na(sig_meta$FDR_LSEA) & sig_meta$FDR_LSEA < fdr_thresh,
      "KS only [distributed effect]", "Neither"
    )
    sig_meta$FDR_fgsea    <- NA_real_
    sig_meta$NES          <- NA_real_
    sig_meta$.is_sig_any  <- !is.na(sig_meta$FDR_LSEA) & sig_meta$FDR_LSEA < fdr_thresh

  } else if (has_fgsea) {
    sig_meta <- lsea_result$fgsea[lsea_result$fgsea$Level == group_col, ]
    sig_meta$Convergence  <- ifelse(
      !is.na(sig_meta$FDR_fgsea) & sig_meta$FDR_fgsea < fdr_thresh,
      "fgsea only [extreme-driven]", "Neither"
    )
    sig_meta$FDR_LSEA         <- NA_real_
    sig_meta$DirectionalScore <- NA_real_
    sig_meta$.is_sig_any      <- !is.na(sig_meta$FDR_fgsea) & sig_meta$FDR_fgsea < fdr_thresh

  } else {
    message("plot_distribution: no LSEA results found.")
    return(NULL)
  }

  # Filter to current level if combined
  if (has_both && "Level" %in% names(sig_meta))
    sig_meta <- sig_meta[sig_meta$Level == group_col, ]

  # Valid groups (min_n)
  valid_groups <- sig_meta$Group[!is.na(sig_meta$N_group) &
                                   sig_meta$N_group >= min_n]

  # -- Prepare plot data ------------------------------------------------------
  plot_df <- data[!is.na(data[[fc_col]]) & !is.na(data[[group_col]]), ]
  plot_df$Group <- as.character(plot_df[[group_col]])
  plot_df <- merge(plot_df,
                   sig_meta[, intersect(names(sig_meta),
                                        c("Group", "Convergence", ".is_sig_any",
                                          "FDR_LSEA", "FDR_fgsea", "NES",
                                          "DirectionalScore", "N_group"))],
                   by = "Group", all.x = TRUE)

  plot_df$Convergence  <- ifelse(is.na(plot_df$Convergence),  "Neither", plot_df$Convergence)
  plot_df$.is_sig_any  <- ifelse(is.na(plot_df$.is_sig_any),  FALSE,     plot_df$.is_sig_any)
  plot_df <- plot_df[plot_df$Group %in% valid_groups, ]

  if (sig_only) plot_df <- plot_df[plot_df$.is_sig_any, ]
  if (nrow(plot_df) == 0L) {
    message("plot_distribution: no groups with n >= ", min_n, " for ", group_col)
    return(NULL)
  }

  # -- Order groups by median logFC -------------------------------------------
  med_order <- tapply(plot_df[[fc_col]], plot_df$Group, stats::median, na.rm = TRUE)
  set_order <- names(sort(med_order))
  plot_df$Group <- factor(plot_df$Group, levels = set_order)

  # -- Y axis limits with label room ------------------------------------------
  y_range      <- range(plot_df[[fc_col]], na.rm = TRUE)
  y_span       <- diff(y_range)
  y_whisker_max <- max(tapply(plot_df[[fc_col]], plot_df$Group, function(x) {
    stats::quantile(x, 0.75, na.rm = TRUE) + 1.5 * stats::IQR(x, na.rm = TRUE)
  }), na.rm = TRUE)
  y_top         <- max(y_range[2], y_whisker_max, na.rm = TRUE)
  label_gap     <- y_span * 0.10
  label_room    <- if (label_angle == 0) y_span * 0.45 else y_span * 0.55
  y_lbl         <- y_top + label_gap
  y_upper_limit <- y_top + label_room
  y_lower_limit <- y_range[1] - y_span * 0.18

  # -- FDR labels for significant sets ----------------------------------------
  fdr_labels <- sig_meta[sig_meta$.is_sig_any & sig_meta$Group %in% valid_groups, ]
  if (nrow(fdr_labels) > 0) {
    fdr_labels$fdr_val <- ifelse(
      fdr_labels$Convergence == "fgsea only [extreme-driven]",
      fdr_labels$FDR_fgsea,
      fdr_labels$FDR_LSEA
    )
    fdr_labels$label <- mapply(function(conv, fdr, ds, nes) {
      fdr_str <- formatC(fdr, format = "e", digits = 1)
      if (conv == "fgsea only [extreme-driven]") {
        sprintf("FDR=%s\nNES=%+.2f", fdr_str, nes)
      } else if (conv == "KS+fgsea [strongest]") {
        sprintf("FDR=%s\nDS=%+.2f\nNES=%+.2f", fdr_str, ds, nes)
      } else {
        sprintf("FDR=%s\nDS=%+.2f", fdr_str, ds)
      }
    }, fdr_labels$Convergence, fdr_labels$fdr_val,
       fdr_labels$DirectionalScore, fdr_labels$NES)

    fdr_labels$Group <- factor(fdr_labels$Group, levels = set_order)
    fdr_labels$y_lbl <- y_lbl
  }

  # -- Build plot -------------------------------------------------------------
  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = Group, y = .data[[fc_col]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey45", linewidth = 0.5)

  # Boxplot layers: non-sig (thin border) then sig (red border)
  avail_conv <- names(pal_conv)
  for (cat in avail_conv) {
    df_ns <- plot_df[plot_df$Convergence == cat & !plot_df$.is_sig_any, ]
    if (nrow(df_ns) > 0)
      p <- p + ggplot2::geom_boxplot(
        data = df_ns, color = pal_conv[[cat]], fill = "white",
        linewidth = 0.38, outlier.shape = NA, width = 0.65
      )
    df_s <- plot_df[plot_df$Convergence == cat & plot_df$.is_sig_any, ]
    if (nrow(df_s) > 0)
      p <- p + ggplot2::geom_boxplot(
        data = df_s, color = "#E63946",
        fill = adjustcolor(pal_conv[[cat]], alpha.f = 0.13),
        linewidth = 1.1, outlier.shape = NA, width = 0.65
      )
  }

  # Jittered points
  p <- p + ggplot2::geom_jitter(
    ggplot2::aes(color = Convergence),
    width = 0.18, size = 0.85, alpha = 0.42
  )

  # FDR labels
  if (nrow(fdr_labels) > 0)
    p <- p + ggplot2::geom_label(
      data          = fdr_labels,
      ggplot2::aes(x = Group, y = y_lbl, label = label),
      angle         = label_angle,
      hjust         = if (label_angle == 90) 0 else 0.5,
      vjust         = 0,
      size          = if (sig_only) {
                         if (label_angle == 0) 3.4 else 3.0
                       } else {
                         if (label_angle == 0) 2.5 else 2.8
                       },
      color         = "#E63946",
      fontface      = "bold",
      fill          = "white",
      label.size    = if (sig_only) 0.15 else 0.12,
      label.padding = ggplot2::unit(if (sig_only) 0.20 else 0.12, "lines"),
      inherit.aes   = FALSE
    )

  # Active palette — drop unused levels
  used_conv <- unique(plot_df$Convergence)
  p <- p +
    ggplot2::scale_color_manual(
      values = pal_conv[names(pal_conv) %in% used_conv],
      drop   = TRUE, name = "Convergence"
    ) +
    ggplot2::coord_cartesian(
      ylim = c(y_lower_limit, y_upper_limit), clip = "off"
    ) +
    ggplot2::labs(
      title    = paste0(case_lbl, " - ", ref_lbl),
      subtitle = paste0(
        "Distribution LSEA \u00b7 ", group_col,
        if (sig_only) " \u00b7 Significant sets only" else "",
        " \u00b7 Red border = significant \u00b7 Fill = detection engine"
      ),
      x       = NULL,
      y       = expression(log[2]*FC),
      caption = paste0(
        "Red border (lw 1.1) = FDR < ", fdr_thresh,
        " \u00b7 Fill: red = KS+fgsea \u00b7 celeste = KS only \u00b7 purple = fgsea only",
        "\nDS = DirectionalScore (standardized mean diff)"
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle      = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey45"),
      plot.caption       = ggplot2::element_text(size = 7, color = "grey55", hjust = 0),
      axis.text.x        = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                  hjust = 1, size = 8,
                                                  face = "bold"),
      axis.text.y        = ggplot2::element_text(size = 9),
      panel.grid.major.y = ggplot2::element_line(color = "grey93"),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_text(face = "bold", size = 9),
      legend.text        = ggplot2::element_text(size = 8),
      plot.margin        = ggplot2::margin(t = 35, r = 12, b = 10, l = 12)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      nrow         = 2,
      override.aes = list(linewidth = 1.1, fill = "white", size = 3.5)
    ))

  attr(p, "n_sets")    <- length(set_order)
  attr(p, "sig_only")  <- sig_only
  p
}
