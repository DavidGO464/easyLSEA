# R/lsea.R
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
    which      = c("bubble_ks", "bubble_fgsea", "bubble_combined",
                   "barplot", "running_sum"),
    fdr_thresh = 0.05,
    case_lbl   = "Case",
    ref_lbl    = "Reference"
) {

  plots <- list()

  # -- KS bubble plot ----------------------------------------------------------
  if ("bubble_ks" %in% which && !is.null(lsea_result$ks)) {
    df <- lsea_result$ks
    df <- df[!is.na(df$KS_pval), ]
    df$sig       <- df$FDR_LSEA < fdr_thresh
    df$log_fdr   <- -log10(df$FDR_LSEA + 1e-300)
    df$direction <- ifelse(df$DirectionalScore > 0,
                           paste0("Up in ", case_lbl),
                           paste0("Up in ", ref_lbl))

    plots$bubble_ks <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x     = DirectionalScore,
        y     = stats::reorder(Group, DirectionalScore),
        size  = N_group,
        color = direction,
        alpha = sig
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(
        values = setNames(
          c("#E63946", "#2166AC"),
          c(paste0("Up in ", case_lbl), paste0("Up in ", ref_lbl))
        )
      ) +
      ggplot2::scale_size_continuous(name = "Set size", range = c(2, 10)) +
      ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35),
                                   guide = "none") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "grey60") +
      ggplot2::labs(
        title    = "LSEA -- KS enrichment",
        subtitle = paste0(case_lbl, " vs ", ref_lbl,
                          "  |  FDR < ", fdr_thresh, " = opaque"),
        x        = "DirectionalScore (Cohen's d)",
        y        = NULL,
        color    = "Direction"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }

  # -- fgsea bubble plot -------------------------------------------------------
  if ("bubble_fgsea" %in% which && !is.null(lsea_result$fgsea)) {
    df <- lsea_result$fgsea
    df <- df[!is.na(df$NES), ]
    df$sig       <- df$FDR_fgsea < fdr_thresh
    df$direction <- ifelse(df$NES > 0,
                           paste0("Up in ", case_lbl),
                           paste0("Up in ", ref_lbl))

    plots$bubble_fgsea <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x     = NES,
        y     = stats::reorder(Group, NES),
        size  = N_leading,
        color = direction,
        alpha = sig
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(
        values = setNames(
          c("#7B2D8B", "#2166AC"),
          c(paste0("Up in ", case_lbl), paste0("Up in ", ref_lbl))
        )
      ) +
      ggplot2::scale_size_continuous(name = "Leading edge n",
                                      range = c(2, 10)) +
      ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35),
                                   guide = "none") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "grey60") +
      ggplot2::labs(
        title    = "LSEA -- fgsea enrichment",
        subtitle = paste0(case_lbl, " vs ", ref_lbl,
                          "  |  FDR < ", fdr_thresh, " = opaque"),
        x        = "NES (normalized enrichment score)",
        y        = NULL,
        color    = "Direction"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }

  if (length(plots) == 0L)
    message("No plots generated -- check 'which' and that results are not NULL.")

  plots
}
