# R/easyLSEA.R
# Main wrapper — orchestrates annotation, LSEA, and chain analysis.

#' Lipid Set Enrichment Analysis — full pipeline
#'
#' One-call interface to the complete easyLSEA workflow:
#' lipid annotation, KS and/or fgsea enrichment across three biological
#' levels (class, LIPID MAPS category, functional category), and fatty
#' acid chain analysis. Returns a structured \code{easyLSEA_result} object
#' that can be plotted and exported.
#'
#' @param data A \code{data.frame} with at least a lipid name column and
#'   a numeric fold-change column. Additional columns (p-values,
#'   confidence ranks, abundance) are used when present.
#' @param lipid_col Character(1). Name of the lipid identifier column.
#'   Default: \code{"LipidName"}.
#' @param fc_col Character(1). Name of the log2 fold-change column.
#'   Default: \code{"logFC"}.
#' @param pval_col Character(1) or \code{NULL}. Name of the raw p-value
#'   column. Used to compute the pi-value rank metric for fgsea.
#'   Default: \code{"P.Value"}.
#' @param case_lbl Character(1). Label for the case group, used in output
#'   tables and plot titles. Default: \code{"Case"}.
#' @param ref_lbl Character(1). Label for the reference group.
#'   Default: \code{"Reference"}.
#' @param engine Character(1). Enrichment engine: \code{"ks"},
#'   \code{"fgsea"}, or \code{"both"}. Default: \code{"both"}.
#' @param annotator Character(1). Lipid annotation method:
#'   \code{"internal"} (built-in regex classifier) or
#'   \code{"lipidAnnotator"} (requires optional package).
#'   Default: \code{"internal"}.
#' @param run_chains Logical(1). Whether to run fatty acid chain analysis
#'   in addition to LSEA. Default: \code{TRUE}.
#' @param min_rank Character(1). Minimum confidence rank for chain analysis.
#'   Ranks are ordered \code{A > B > C > D > E > P}. Lipids with rank lower
#'   than \code{min_rank} (and rank \code{"P"} or \code{NA}) are excluded
#'   from chain parsing. Default: \code{"E"} (include all except P and NA).
#'   Only used when \code{run_chains = TRUE} and a rank column is present.
#' @param group_cols Character vector. Grouping columns to test in LSEA.
#'   If \code{NULL} (default), uses the three standard levels:
#'   \code{LipidClass}, \code{LipidCategory_LMAPS},
#'   \code{LipidCategory_functional}.
#' @param min_n Integer(1). Minimum set size to test. Default: \code{3L}.
#' @param n_perm Integer(1). KS permutations for \code{DS_perm_pval}.
#'   Default: \code{2000L}.
#' @param fgsea_nperm Integer(1). fgsea Monte Carlo permutations.
#'   Default: \code{10000L}.
#' @param plots Logical(1). Whether to generate \pkg{ggplot2} objects.
#'   Set to \code{FALSE} to skip plotting and reduce runtime.
#'   Default: \code{TRUE}.
#' @param bubble_label Character vector. Which statistics to show next to
#'   each bubble in the LSEA bubble plots. Any subset of \code{"FDR"},
#'   \code{"DS"} (KS only), \code{"NES"} (fgsea only), and \code{"n"}.
#'   Use fewer to shorten labels. Default: all four.
#' @param output Character(1). Return format when both modules run:
#'   \code{"combined"} returns a single \code{easyLSEA_result};
#'   \code{"separate"} returns a named list with elements \code{lsea}
#'   and \code{chains}. Default: \code{"combined"}.
#' @param seed Integer(1) or \code{NULL}. RNG seed for reproducibility.
#'   Passed to \code{\link[withr]{with_seed}} — does not alter the
#'   user's global RNG state. Default: \code{42L}.
#' @param verbose Logical(1). Print progress messages. Default: \code{TRUE}.
#'
#' @return An object of class \code{easyLSEA_result}: a named list with
#'   five slots.
#'   \describe{
#'     \item{\code{$meta}}{Named list: call, date, labels, engine, counts.}
#'     \item{\code{$lsea}}{Named list: \code{results} (data.frame with KS
#'       and/or fgsea statistics), \code{combined} (merged table with
#'       Convergence column).}
#'     \item{\code{$chains}}{Named list: \code{parsed} and \code{summary}
#'       from \code{\link{parse_lipid_chains}}, or \code{NULL} if
#'       \code{run_chains = FALSE}.}
#'     \item{\code{$plots}}{Named list of \code{ggplot} objects, or
#'       \code{NULL} if \code{plots = FALSE}.}
#'     \item{\code{$input}}{Named list: \code{data} (annotated input),
#'       \code{group_cols}.}
#'   }
#'   When \code{output = "separate"}, returns
#'   \code{list(lsea = ..., chains = ...)} instead.
#'
#' @seealso
#'   \code{\link{annotate_lipids}} for standalone annotation,
#'   \code{\link{run_lsea}} for the enrichment engine,
#'   \code{\link{parse_lipid_chains}} for chain analysis,
#'   \code{\link{plot_lsea}}, \code{\link{plot_chains}},
#'   \code{export_lsea()} to save results.
#'
#' @examples
#' data("lipid_example", package = "easyLSEA")
#'
#' result <- easyLSEA(
#'   data      = lipid_example,
#'   lipid_col = "LipidName",
#'   fc_col    = "logFC",
#'   case_lbl  = "NASH",
#'   ref_lbl   = "Control",
#'   engine    = "ks",
#'   plots     = FALSE
#' )
#'
#' print(result)
#' head(result$lsea$results)
#'
#' @export
easyLSEA <- function(
    data,
    lipid_col   = "LipidName",
    fc_col      = "logFC",
    pval_col    = "P.Value",
    case_lbl    = "Case",
    ref_lbl     = "Reference",
    engine      = c("both", "ks", "fgsea"),
    annotator   = c("internal", "lipidAnnotator"),
    run_chains  = TRUE,
    min_rank    = "E",
    group_cols  = NULL,
    min_n       = 3L,
    n_perm      = 2000L,
    fgsea_nperm = 10000L,
    plots        = TRUE,
    bubble_label = c("FDR", "DS", "NES", "n"),
    output       = c("combined", "separate"),
    seed         = 42L,
    verbose      = TRUE
) {

  engine    <- match.arg(engine)
  annotator <- match.arg(annotator)
  output    <- match.arg(output)
  bubble_label <- match.arg(bubble_label, several.ok = TRUE)

  call_time <- Sys.time()

  # -- 1. Validate input ------------------------------------------------------
  .validate_input(data, lipid_col, fc_col)

  # -- 2. Annotate lipids -----------------------------------------------------
  if (verbose) message("[1/4] Annotating lipids...")

  annotated <- annotate_lipids(
    data          = data,
    lipid_col     = lipid_col,
    method        = annotator,
    verbose       = verbose
  )

  # -- 3. Set grouping levels -------------------------------------------------
  if (is.null(group_cols)) {
    group_cols <- intersect(
      c("LipidClass", "LipidCategory_LMAPS", "LipidCategory_functional"),
      names(annotated)
    )
  }

  if (length(group_cols) == 0L)
    stop("No grouping columns found after annotation. ",
         "Check that annotation produced LipidClass columns.",
         call. = FALSE)

  # -- 4. LSEA ----------------------------------------------------------------
  if (verbose) message("[2/4] Running LSEA (engine = '", engine, "')...")

  lsea_out <- run_lsea(
    data        = annotated,
    group_cols  = group_cols,
    fc_col      = fc_col,
    pval_col    = pval_col,
    lipid_id_col = lipid_col,
    case_lbl    = case_lbl,
    ref_lbl     = ref_lbl,
    engine      = engine,
    min_n       = min_n,
    n_perm      = n_perm,
    fgsea_nperm = fgsea_nperm,
    seed        = seed,
    verbose     = FALSE
  )

  # -- 5. Chain analysis ------------------------------------------------------
  chains_out <- NULL

  if (run_chains) {
    if (verbose) message("[3/4] Running chain analysis...")
    chains_out <- tryCatch(
      parse_lipid_chains(
        data      = annotated,
        lipid_col = lipid_col,
        class_col = "LipidClass",
        min_rank  = min_rank
      ),
      error = function(e) {
        warning("Chain analysis failed: ", conditionMessage(e),
                "\nReturning NULL for $chains.", call. = FALSE)
        NULL
      }
    )
  }

  # -- 6. Plots ---------------------------------------------------------------
  plot_list <- NULL

  if (plots) {
    if (verbose) message("[4/4] Generating plots...")
    plot_list <- list()

    plot_list$lsea <- tryCatch(
      plot_lsea(
        lsea_result  = lsea_out,
        case_lbl     = case_lbl,
        ref_lbl      = ref_lbl,
        bubble_label = bubble_label
      ),
      error = function(e) {
        warning("LSEA plot generation failed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )

    # Distribution boxplots — one per grouping level
    dist_plots <- list()
    for (gcol in group_cols) {
      lbl <- switch(gcol,
                    LipidClass               = "01_Class",
                    LipidCategory_LMAPS      = "02_LMAPS",
                    LipidCategory_functional = "03_Functional",
                    gcol
      )
      key <- paste0("dist_", lbl)
      dist_plots[[key]] <- tryCatch(
        plot_distribution(
          data        = annotated,
          lsea_result = lsea_out,
          group_col   = gcol,
          fc_col      = fc_col,
          case_lbl    = case_lbl,
          ref_lbl     = ref_lbl
        ),
        error = function(e) {
          warning("Distribution plot failed for ", gcol, ": ",
                  conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      # sig_only variant
      key_sig <- paste0("dist_sig_", lbl)
      dist_plots[[key_sig]] <- tryCatch(
        plot_distribution(
          data        = annotated,
          lsea_result = lsea_out,
          group_col   = gcol,
          fc_col      = fc_col,
          case_lbl    = case_lbl,
          ref_lbl     = ref_lbl,
          sig_only    = TRUE
        ),
        error = function(e) NULL
      )
    }
    plot_list$lsea <- c(plot_list$lsea,
                        Filter(Negate(is.null), dist_plots))

    if (!is.null(chains_out) && !is.null(chains_out$parsed) &&
        nrow(chains_out$parsed) > 0L) {
      plot_list$chains <- tryCatch(
        plot_chains(
          chains_result = chains_out,
          case_lbl      = case_lbl,
          ref_lbl       = ref_lbl
        ),
        error = function(e) {
          warning("Chain plot generation failed: ", conditionMessage(e),
                  call. = FALSE)
          NULL
        }
      )
    }
  }

  # -- 7. Assemble result object ----------------------------------------------
  result <- structure(
    list(
      meta = list(
        call       = match.call(),
        date       = call_time,
        case_lbl   = case_lbl,
        ref_lbl    = ref_lbl,
        engine     = engine,
        annotator  = annotator,
        n_lipids   = nrow(data),
        group_cols = group_cols
      ),
      lsea   = lsea_out,
      chains = chains_out,
      plots  = plot_list,
      input  = list(
        data       = annotated,
        group_cols = group_cols
      )
    ),
    class = "easyLSEA_result"
  )

  if (verbose) message("Done. Use print(), plot(), or export_lsea() to explore results.")

  # -- 8. Return --------------------------------------------------------------
  if (output == "separate") {
    return(list(
      lsea   = result$lsea,
      chains = result$chains
    ))
  }

  result
}

# -- S3 methods ---------------------------------------------------------------

#' Print method for easyLSEA_result
#'
#' @param x An \code{easyLSEA_result} object.
#' @param ... Ignored.
#' @return Invisibly returns the input \code{easyLSEA_result} object
#'   (\code{x}). Called for its side effect of printing a formatted
#'   summary of the enrichment results to the console.
#' @export
print.easyLSEA_result <- function(x, ...) {

  cat("-- easyLSEA result ", strrep("-", 44), "\n", sep = "")
  cat(sprintf("  Comparison : %s vs %s\n", x$meta$case_lbl, x$meta$ref_lbl))
  cat(sprintf("  Date       : %s\n", format(x$meta$date, "%Y-%m-%d %H:%M")))
  cat(sprintf("  Engine     : %s\n", x$meta$engine))
  cat(sprintf("  Lipids     : %d\n", x$meta$n_lipids))
  cat(sprintf("  Levels     : %s\n", paste(x$meta$group_cols, collapse = ", ")))

  if (!is.null(x$lsea$ks)) {
    n_sig_ks <- sum(x$lsea$ks$FDR_LSEA < 0.05, na.rm = TRUE)
    cat("-- LSEA (KS) ", strrep("-", 47), "\n", sep = "")
    cat(sprintf("  Significant sets (FDR < 0.05): %d\n", n_sig_ks))
    if (n_sig_ks > 0L) {
      top <- x$lsea$ks[order(x$lsea$ks$KS_pval), ][1L, ]
      cat(sprintf("  Top hit: %s  (DS = %.2f, FDR = %.3f)\n",
                  top$Group, top$DirectionalScore, top$FDR_LSEA))
    }
  }

  if (!is.null(x$lsea$fgsea)) {
    n_sig_fg <- sum(x$lsea$fgsea$FDR_fgsea < 0.05, na.rm = TRUE)
    cat("-- LSEA (fgsea) ", strrep("-", 44), "\n", sep = "")
    cat(sprintf("  Significant sets (FDR < 0.05): %d\n", n_sig_fg))
    if (n_sig_fg > 0L) {
      top <- x$lsea$fgsea[order(x$lsea$fgsea$fgsea_pval), ][1L, ]
      cat(sprintf("  Top hit: %s  (NES = %.2f, FDR = %.3f)\n",
                  top$Group, top$NES, top$FDR_fgsea))
    }
  }

  if (!is.null(x$chains)) {
    n_parsed <- if (!is.null(x$chains$parsed)) nrow(x$chains$parsed) else 0L
    cat("-- Chain analysis ", strrep("-", 42), "\n", sep = "")
    cat(sprintf("  Chain observations parsed: %d\n", n_parsed))
  }

  cat(strrep("-", 60), "\n", sep = "")
  cat("  Use result$lsea$results, result$chains$parsed to access data.\n")
  cat("  Use export_lsea() to save all outputs.\n")

  invisible(x)
}

#' Summary method for easyLSEA_result
#'
#' @param object An \code{easyLSEA_result} object.
#' @param padj_cutoff Numeric(1). FDR threshold for significant sets.
#'   Default: \code{0.05}.
#' @param ... Ignored.
#' @return Invisibly returns the input \code{easyLSEA_result} object
#'   (\code{object}). Called for its side effect of printing a summary
#'   table of the significant lipid sets to the console.
#' @export
summary.easyLSEA_result <- function(object, padj_cutoff = 0.05, ...) {

  cat("easyLSEA summary\n\n")

  if (!is.null(object$lsea$combined)) {
    df <- object$lsea$combined
    sig <- df[!is.na(df$FDR_LSEA) & df$FDR_LSEA < padj_cutoff, ]
    cat(sprintf("Significant sets (KS FDR < %.2f): %d / %d tested\n",
                padj_cutoff, nrow(sig), nrow(df)))
    if (nrow(sig) > 0L) {
      cat("\n")
      print(sig[order(sig$KS_pval),
                intersect(c("Level", "Group", "N_group",
                            "DirectionalScore", "FDR_LSEA",
                            "NES", "FDR_fgsea", "Convergence"),
                          names(sig))],
            row.names = FALSE)
    }
  }

  invisible(object)
}
