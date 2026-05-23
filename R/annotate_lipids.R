# R/annotate_lipids.R
# Public lipid annotation interface.
# Dispatches to the internal regex-based classifier (.annotate_internal())
# or to the optional lipidAnnotator package.

#' Annotate lipid names with class and category information
#'
#' Assigns lipid class (e.g. PC, TG, Cer), full class name, LIPID MAPS
#' structural category, and functional category to each lipid in \code{data}.
#' Returns the input data.frame with annotation columns appended, ready for
#' use in \code{\link{run_lsea}} and \code{\link{parse_lipid_chains}}.
#'
#' @param data A \code{data.frame} with at least one column of lipid names.
#' @param lipid_col Character(1). Name of the column containing lipid
#'   identifiers. Default: \code{"LipidName"}.
#' @param shorthand_col Character(1) or \code{NULL}. Optional column with
#'   shorthand notation used as the primary annotation source when present
#'   (more standardised than common names). Falls back to \code{lipid_col}
#'   if \code{NULL} or column not found. Default: \code{"Shorthand"}.
#' @param method Character(1). Annotation method:
#'   \describe{
#'     \item{\code{"internal"}}{Regex-based hierarchical classifier validated
#'       against LIPID MAPS nomenclature. No external dependencies. Covers
#'       all major lipid classes in untargeted lipidomics (GPL, SL, GL,
#'       FA, ST, acylcarnitines, oxylipins, bile acids). Default.}
#'     \item{\code{"lipidAnnotator"}}{Uses the \pkg{lipidAnnotator} package
#'       (available on GitHub, archived on Zenodo). Must be installed
#'       separately. Provides enhanced structural annotation.}
#'   }
#' @param verbose Logical(1). Print annotation summary (class distribution
#'   and count of unclassified lipids). Default: \code{TRUE}.
#'
#' @return The input \code{data.frame} with five columns appended:
#'   \describe{
#'     \item{\code{LipidClass}}{Abbreviated class (e.g. "PC", "TG", "Cer").}
#'     \item{\code{LipidClass_Full}}{Descriptive class name
#'       (e.g. "Ceramide", "Ether-PC").}
#'     \item{\code{LipidCategory_LMAPS}}{LIPID MAPS structural category
#'       (e.g. "Glycerophospholipids", "Sphingolipids").}
#'     \item{\code{LipidCategory_functional}}{Functional category, with
#'       Oxylipins and Bile Acids as standalone groups rather than nested
#'       under Fatty Acyls.}
#'     \item{\code{LipidCategory}}{Simplified category for plotting:
#'       same as \code{LipidCategory_functional} except Saccharolipids
#'       are shown as "Glycolipids".}
#'   }
#'   Lipids that cannot be classified receive \code{LipidClass = "Unknown"}.
#'
#' @seealso \code{\link{run_lsea}}, \code{\link{parse_lipid_chains}}
#'
#' @importFrom utils head
#'
#' @examples
#' df <- data.frame(
#'   LipidName = c("PC 36:2", "TG(54:3)", "SM d18:1/16:0",
#'                 "Cer(d18:1/24:0)", "LPC 18:0", "CE 18:1"),
#'   logFC     = c(1.2, -0.8, 0.5, -1.1, 0.3, 0.9),
#'   stringsAsFactors = FALSE
#' )
#'
#' annotated <- annotate_lipids(df)
#' annotated[, c("LipidName", "LipidClass", "LipidCategory")]
#'
#' @export
annotate_lipids <- function(
    data,
    lipid_col     = "LipidName",
    shorthand_col = "Shorthand",
    method        = c("internal", "lipidAnnotator"),
    verbose       = TRUE
) {

  method <- match.arg(method)

  # -- Input checks -----------------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.", call. = FALSE)
  if (!lipid_col %in% names(data))
    stop("Column ", sQuote(lipid_col), " not found in 'data'.", call. = FALSE)

  # Determine which column to annotate from
  # Prefer shorthand (more standardised) when available
  annot_col <- if (!is.null(shorthand_col) &&
                   shorthand_col %in% names(data)) {
    shorthand_col
  } else {
    lipid_col
  }

  # -- Dispatch ---------------------------------------------------------------
  if (method == "lipidAnnotator") {

    # annotate_lipid() is now part of easyLSEA (absorbed from lipidAnnotator).
    # Uses the full LIPID MAPS parser with canonical shorthand support.
    annot <- annotate_lipid(data[[annot_col]],
                            detail   = "full",
                            no_match = "ignore")

    # Map annotate_lipid() columns to easyLSEA standard
    data$LipidClass      <- annot$Class
    data$LipidClass_Full <- annot$lm_class_name

    # lm_category_name = LIPID MAPS structural category
    data$LipidCategory_LMAPS <- annot$lm_category_name

    # Functional category: Oxylipins and Bile Acids as standalone groups
    data$LipidCategory_functional <- mapply(
      .get_lipid_category_functional,
      annot$lm_category_name,
      annot$Class,
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )

    data$LipidCategory <- mapply(
      .get_lipid_category,
      data$LipidCategory_LMAPS,
      data$LipidCategory_functional,
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )

    # Replace NA LipidClass with "Unknown" for consistency with internal method
    data$LipidClass[is.na(data$LipidClass) |
                      data$LipidClass == ""] <- "Unknown"

  } else {

    # Internal regex-based classifier
    annot <- .annotate_internal(data[[annot_col]])
    data$LipidClass               <- annot$LipidClass
    data$LipidClass_Full          <- annot$LipidClass_Full
    data$LipidCategory_LMAPS      <- annot$LipidCategory_LMAPS
    data$LipidCategory_functional <- annot$LipidCategory_functional
    data$LipidCategory            <- annot$LipidCategory

  }

  # -- Summary ----------------------------------------------------------------
  if (verbose) {
    n_total   <- nrow(data)
    n_unknown <- sum(data$LipidClass == "Unknown", na.rm = TRUE)

    message(sprintf(
      "annotate_lipids(): %d lipids annotated | %d unclassified (%.1f%%)",
      n_total, n_unknown, 100 * n_unknown / n_total
    ))

    if (n_unknown > 0L) {
      unknown_names <- data[[lipid_col]][data$LipidClass == "Unknown"]
      message("  Unclassified: ",
              paste(head(unknown_names, 10L), collapse = ", "),
              if (n_unknown > 10L) paste0(" ... and ", n_unknown - 10L,
                                          " more") else "")
    }

    cls_tbl <- sort(table(data$LipidClass), decreasing = TRUE)
    message("  Class distribution: ",
            paste(names(cls_tbl), cls_tbl, sep = "=", collapse = " | "))
  }

  data
}
