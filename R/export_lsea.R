# R/export_lsea.R
# Export easyLSEA_result objects to CSV, Excel, PDF/PNG, and HTML.

#' Export easyLSEA results to disk
#'
#' Saves the contents of an \code{\link{easyLSEA}} result object to a
#' timestamped output folder. Supported formats: CSV tables, a multi-sheet
#' Excel workbook, PDF or PNG plots, and a standalone HTML report.
#' Any combination of formats can be requested in a single call.
#'
#' @param result An \code{easyLSEA_result} object returned by
#'   \code{\link{easyLSEA}}, or a named list with elements \code{lsea}
#'   and/or \code{chains} (output of \code{output = "separate"}).
#' @param dir Character(1). Base directory where the output folder will be
#'   created. Default: current working directory (\code{"."}).
#' @param prefix Character(1). Prefix for the output folder name. The folder
#'   is named \code{<prefix>_<YYYY-MM-DD>/}. Default: \code{"easyLSEA"}.
#' @param format Character vector. One or more of \code{"csv"},
#'   \code{"excel"}, \code{"pdf"}, \code{"png"}, \code{"html"}.
#'   Default: \code{c("csv", "excel", "pdf")}.
#' @param overwrite Logical(1). If \code{TRUE}, an existing output folder with
#'   the same name is overwritten. Default: \code{FALSE}.
#' @param plot_width Numeric(1). Plot width in inches. Default: \code{8}.
#' @param plot_height Numeric(1). Plot height in inches. Default: \code{6}.
#' @param plot_dpi Integer(1). Resolution for PNG output. Default: \code{300L}.
#' @param verbose Logical(1). Print progress messages. Default: \code{TRUE}.
#'
#' @return Invisibly returns a named character vector of all file paths
#'   created. Useful for programmatic use or verification.
#'
#' @details
#' ## Output folder structure
#' \preformatted{
#' <prefix>_<YYYY-MM-DD>/
#'   tables/
#'     lsea_results_ks.csv
#'     lsea_results_fgsea.csv
#'     lsea_combined.csv
#'     chain_results.csv
#'     chain_parsed.csv
#'   plots/
#'     lsea/
#'       bubble_ks.pdf
#'       bubble_fgsea.pdf
#'     chains/
#'       tile/
#'         tile_PC.pdf
#'         tile_TG.pdf  ...
#'       trend/
#'         trend_length_PC.pdf
#'         trend_unsat_PC.pdf  ...
#'   results.xlsx
#'   report.html
#' }
#'
#' ## Dependencies for optional formats
#' Excel export requires \pkg{openxlsx} (\code{install.packages("openxlsx")}).
#' HTML export requires \pkg{rmarkdown} and \pkg{knitr}.
#'
#' @seealso \code{\link{easyLSEA}}, \code{\link{run_lsea}},
#'   \code{\link{parse_lipid_chains}}
#'
#' @examples
#' data("lipid_example", package = "easyLSEA")
#'
#' result <- suppressWarnings(suppressMessages(easyLSEA(
#'   data      = lipid_example,
#'   engine    = "ks",
#'   n_perm    = 100L,
#'   plots     = FALSE,
#'   verbose   = FALSE
#' )))
#'
#' \dontrun{
#' # Export CSV and PDF to a temp folder
#' paths <- export_lsea(result, dir = tempdir(), format = c("csv", "pdf"))
#' paths
#' }
#'
#' @importFrom utils write.csv
#' @importFrom grDevices pdf png dev.off
#'
#' @export
export_lsea <- function(
    result,
    dir         = ".",
    prefix      = "easyLSEA",
    format      = c("csv", "excel", "pdf"),
    overwrite   = FALSE,
    plot_width  = 8,
    plot_height = 6,
    plot_dpi    = 300L,
    verbose     = TRUE
) {

  format <- match.arg(format,
                      choices   = c("csv", "excel", "pdf", "png", "html"),
                      several.ok = TRUE)

  # -- Validate input ---------------------------------------------------------
  if (!inherits(result, "easyLSEA_result") && !is.list(result))
    stop("'result' must be an easyLSEA_result object or a named list.",
         call. = FALSE)

  # -- Create output folder ---------------------------------------------------
  folder_name <- paste0(prefix, "_", format(Sys.Date(), "%Y-%m-%d"))
  out_dir     <- file.path(dir, folder_name)

  if (dir.exists(out_dir) && !overwrite)
    stop("Output folder already exists: ", out_dir, "\n",
         "Set overwrite = TRUE to replace it.", call. = FALSE)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (verbose) message("Output folder: ", out_dir)

  created_files <- character(0)

  # -- Extract data slots -----------------------------------------------------
  lsea_slot   <- if (inherits(result, "easyLSEA_result")) result$lsea
                 else result$lsea
  chains_slot <- if (inherits(result, "easyLSEA_result")) result$chains
                 else result$chains
  plots_slot  <- if (inherits(result, "easyLSEA_result")) result$plots
                 else NULL

  # -- CSV export -------------------------------------------------------------
  if ("csv" %in% format) {
    tbl_dir <- file.path(out_dir, "tables")
    dir.create(tbl_dir, showWarnings = FALSE)

    tables <- list(
      lsea_results_ks     = if (!is.null(lsea_slot$ks))       lsea_slot$ks,
      lsea_results_fgsea  = if (!is.null(lsea_slot$fgsea))    lsea_slot$fgsea,
      lsea_combined       = if (!is.null(lsea_slot$combined))  lsea_slot$combined,
      chain_results       = if (!is.null(chains_slot$summary)) chains_slot$summary,
      chain_parsed        = if (!is.null(chains_slot$parsed))  chains_slot$parsed
    )
    tables <- Filter(Negate(is.null), tables)

    for (tbl_name in names(tables)) {
      path <- file.path(tbl_dir, paste0(tbl_name, ".csv"))
      write.csv(tables[[tbl_name]], path, row.names = FALSE)
      created_files <- c(created_files, path)
    }

    if (verbose)
      message("  CSV: ", length(tables), " table(s) -> tables/")
  }

  # -- Excel export -----------------------------------------------------------
  if ("excel" %in% format) {

    if (!requireNamespace("openxlsx", quietly = TRUE))
      stop("Package 'openxlsx' is required for Excel export.\n",
           "Install with: install.packages('openxlsx')", call. = FALSE)

    wb <- openxlsx::createWorkbook()

    # Header style
    hdr <- openxlsx::createStyle(
      fontColour = "#FFFFFF", fgFill = "#2C3E50",
      halign = "left", textDecoration = "bold"
    )

    .add_sheet <- function(wb, name, data) {
      if (is.null(data) || nrow(data) == 0L) return(invisible(NULL))
      openxlsx::addWorksheet(wb, name)
      openxlsx::writeData(wb, name, data, headerStyle = hdr)
      openxlsx::freezePane(wb, name, firstRow = TRUE)
      openxlsx::setColWidths(wb, name, cols = seq_len(ncol(data)),
                              widths = "auto")
    }

    .add_sheet(wb, "LSEA_KS",      lsea_slot$ks)
    .add_sheet(wb, "LSEA_fgsea",   lsea_slot$fgsea)
    .add_sheet(wb, "LSEA_combined",lsea_slot$combined)
    .add_sheet(wb, "Chain_summary", chains_slot$summary)
    .add_sheet(wb, "Chain_parsed",  chains_slot$parsed)

    # Metadata sheet
    openxlsx::addWorksheet(wb, "Metadata")
    if (inherits(result, "easyLSEA_result")) {
      meta_df <- data.frame(
        Field = c("Date", "Case", "Reference", "Engine",
                  "Annotator", "N_lipids"),
        Value = c(
          format(result$meta$date, "%Y-%m-%d %H:%M"),
          result$meta$case_lbl,
          result$meta$ref_lbl,
          result$meta$engine,
          result$meta$annotator,
          as.character(result$meta$n_lipids)
        ),
        stringsAsFactors = FALSE
      )
      openxlsx::writeData(wb, "Metadata", meta_df, headerStyle = hdr)
    }

    xl_path <- file.path(out_dir, "results.xlsx")
    openxlsx::saveWorkbook(wb, xl_path, overwrite = TRUE)
    created_files <- c(created_files, xl_path)

    if (verbose) message("  Excel: results.xlsx")
  }

  # -- PDF / PNG plot export --------------------------------------------------
  plot_formats <- intersect(format, c("pdf", "png"))

  if (length(plot_formats) > 0L && !is.null(plots_slot)) {

    plt_dir <- file.path(out_dir, "plots")

    # Subdirectory routing: lsea/ | chains/tile/ | chains/trend/
    .plot_subdir <- function(name) {
      if (grepl("^bubble.*01_Class",      name)) file.path(plt_dir, "lsea", "01_Class")
      else if (grepl("^bubble.*02_LMAPS", name)) file.path(plt_dir, "lsea", "02_LMAPS")
      else if (grepl("^bubble.*03_Func",  name)) file.path(plt_dir, "lsea", "03_Functional")
      else if (grepl("^dist.*01_Class",   name)) file.path(plt_dir, "lsea", "01_Class")
      else if (grepl("^dist.*02_LMAPS",   name)) file.path(plt_dir, "lsea", "02_LMAPS")
      else if (grepl("^dist.*03_Func",    name)) file.path(plt_dir, "lsea", "03_Functional")
      else if (grepl("^bubble",           name)) file.path(plt_dir, "lsea")
      else if (grepl("^dist_",            name)) file.path(plt_dir, "lsea")
      else if (grepl("^tile_",            name)) file.path(plt_dir, "chains", "tile")
      else if (grepl("^trend_",           name)) file.path(plt_dir, "chains", "trend")
      else                                       plt_dir
    }

    # Build named list: lsea plots + chain plots
    lsea_plots   <- if (!is.null(plots_slot$lsea))   plots_slot$lsea   else list()
    chain_plots  <- if (!is.null(plots_slot$chains)) plots_slot$chains else list()
    all_plots    <- Filter(Negate(is.null), c(lsea_plots, chain_plots))

    if (length(all_plots) == 0L && verbose) {
      message("  Plots: none available (run easyLSEA with plots = TRUE)")
    }

    for (dev in plot_formats) {
      n_saved <- 0L
      for (plt_name in names(all_plots)) {
        plt <- all_plots[[plt_name]]
        if (!inherits(plt, "ggplot")) next

        sub <- .plot_subdir(plt_name)
        dir.create(sub, recursive = TRUE, showWarnings = FALSE)

        path <- file.path(sub, paste0(plt_name, ".", dev))

        tryCatch({
          # Dynamic sizing: dist and bubble plots scale with number of groups
          n_sets   <- attr(plt, "n_sets")
          sig_only <- isTRUE(attr(plt, "sig_only"))
          is_bubble <- grepl("^bubble", plt_name)
          w <- if (!is.null(n_sets) && is_bubble) {
            max(8, 1 + n_sets * 0.45)    # bubble: labels are on the right
          } else if (!is.null(n_sets)) {
            max(8, 1 + n_sets * 0.7)     # dist: 2-line labels need wider spacing
          } else plot_width
          h <- if (!is.null(n_sets) && !is_bubble && sig_only) {
            max(7, n_sets * 1.4)         # sig dist: taller so boxes aren't flat
          } else plot_height

          if (dev == "pdf") {
            grDevices::pdf(path, width = w, height = h)
          } else {
            grDevices::png(path, width = w, height = h,
                           units = "in", res = plot_dpi)
          }
          print(plt)
          grDevices::dev.off()
          created_files <- c(created_files, path)
          n_saved <- n_saved + 1L
        }, error = function(e) {
          grDevices::dev.off()
          warning("Could not save plot '", plt_name, "': ",
                  conditionMessage(e), call. = FALSE)
        })
      }
      if (verbose)
        message("  ", toupper(dev), ": ", n_saved,
                " plot(s) -> plots/lsea/ | plots/chains/tile/ | plots/chains/trend/")
    }
  } else if (length(plot_formats) > 0L && is.null(plots_slot) && verbose) {
    message("  Plots: skipped (no plots in result -- ",
            "run easyLSEA with plots = TRUE)")
  }

  # -- HTML report ------------------------------------------------------------
  if ("html" %in% format) {

    for (pkg in c("rmarkdown", "knitr")) {
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg, "' is required for HTML export.\n",
             "Install with: install.packages('", pkg, "')", call. = FALSE)
    }

    html_path <- file.path(out_dir, "report.html")

    # Build the report inline (no external .Rmd template needed)
    rmd_content <- .build_report_rmd(result)
    rmd_tmp     <- tempfile(fileext = ".Rmd")
    writeLines(rmd_content, rmd_tmp)

    tryCatch({
      rmarkdown::render(
        input       = rmd_tmp,
        output_file = html_path,
        quiet       = !verbose,
        envir       = new.env(parent = globalenv())
      )
      created_files <- c(created_files, html_path)
      if (verbose) message("  HTML: report.html")
    }, error = function(e) {
      warning("HTML report generation failed: ", conditionMessage(e),
              call. = FALSE)
    })

    unlink(rmd_tmp)
  }

  # -- Summary ----------------------------------------------------------------
  if (verbose) {
    message(sprintf("Done. %d file(s) saved to: %s",
                    length(created_files), out_dir))
  }

  invisible(stats::setNames(created_files, basename(created_files)))
}

# -- Internal: build inline Rmd content for HTML report -----------------------

#' @noRd
.build_report_rmd <- function(result) {

  meta <- if (inherits(result, "easyLSEA_result")) result$meta else list()

  header <- paste0(
    '---\n',
    'title: "easyLSEA Report"\n',
    'date: "', format(Sys.Date(), "%Y-%m-%d"), '"\n',
    'output:\n',
    '  html_document:\n',
    '    toc: true\n',
    '    toc_float: true\n',
    '    theme: flatly\n',
    '---\n\n'
  )

  intro <- paste0(
    '```{r setup, include=FALSE}\n',
    'knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)\n',
    '```\n\n',
    '## Analysis summary\n\n',
    if (!is.null(meta$case_lbl))
      paste0('**Comparison:** ', meta$case_lbl, ' vs ', meta$ref_lbl, '\n\n')
    else '',
    if (!is.null(meta$date))
      paste0('**Date:** ', format(meta$date, "%Y-%m-%d %H:%M"), '\n\n')
    else '',
    if (!is.null(meta$engine))
      paste0('**Engine:** ', meta$engine, '\n\n')
    else ''
  )

  lsea_section <- '## LSEA results\n\n'
  lsea_slot    <- if (inherits(result, "easyLSEA_result")) result$lsea
                  else result$lsea

  if (!is.null(lsea_slot$ks) && nrow(lsea_slot$ks) > 0L) {
    lsea_section <- paste0(
      lsea_section,
      '### KS enrichment\n\n',
      '```{r}\n',
      'ks_tbl <- result_data$lsea$ks\n',
      'knitr::kable(ks_tbl[order(ks_tbl$KS_pval), ],\n',
      '             digits = 4, row.names = FALSE)\n',
      '```\n\n'
    )
  }

  chain_section <- ''
  chains_slot   <- if (inherits(result, "easyLSEA_result")) result$chains
                   else result$chains

  if (!is.null(chains_slot$summary) && nrow(chains_slot$summary) > 0L) {
    chain_section <- paste0(
      '## Chain analysis\n\n',
      '```{r}\n',
      'knitr::kable(result_data$chains$summary,\n',
      '             digits = 4, row.names = FALSE)\n',
      '```\n\n'
    )
  }

  session_section <- paste0(
    '## Session info\n\n',
    '```{r}\n',
    'sessionInfo()\n',
    '```\n'
  )

  paste0(header, intro, lsea_section, chain_section, session_section)
}
