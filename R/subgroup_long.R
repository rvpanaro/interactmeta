#' Reshape wide subgroup data to long format
#'
#' Shipped **interactmeta** datasets are already long (one row per trial ×
#' subgroup). If your source tables are wide (one row per trial, columns
#' suffixed by subgroup such as `yi_Younger` or `n_IV`), use this helper before
#' [cams()] or [escalc()]. Metadata on wide tables (`subgroups`,
#' `subgroup_id_cols`, `subgroup_value_cols`) drives the reshape; you can
#' override those arguments when needed.
#'
#' @param data A wide data frame, or an already long data frame (returned as-is).
#' @param subgroups Modifier labels in column-suffix order (default: dataset attribute).
#' @param id_cols Trial-level columns repeated on each row (default: attribute).
#' @param values Column stems to expand (default: attribute).
#' @return A long data frame with `subgroup` and the value columns.
#' @export
#' @examples
#' # Shipped `steroids` is already long; `subgroup_long()` is for wide source tables.
#' long_steroids <- subgroup_long(steroids)
#' head(long_steroids[, c("study", "subgroup", "n", "treat.events")])
subgroup_long <- function(data,
                          subgroups = NULL,
                          id_cols = NULL,
                          values = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if ("subgroup" %in% names(data)) {
    return(data)
  }
  if (is.null(subgroups)) subgroups <- attr(data, "subgroups")
  if (is.null(id_cols)) id_cols <- attr(data, "subgroup_id_cols")
  if (is.null(values)) values <- attr(data, "subgroup_value_cols")
  if (is.null(subgroups) || is.null(values)) {
    stop(
      "Wide subgroup data must define `subgroups` and value columns; ",
      "supply them as arguments or set attributes on `data`.",
      call. = FALSE
    )
  }
  if (is.null(id_cols)) {
    id_cols <- character(0)
  }
  id_cols <- intersect(id_cols, names(data))
  subgroups <- as.character(subgroups)
  values <- as.character(values)

  pieces <- lapply(seq_len(nrow(data)), function(i) {
    base <- if (length(id_cols)) data[i, id_cols, drop = FALSE] else data[i, , drop = FALSE][0, , drop = FALSE]
    lapply(subgroups, function(sg) {
      out <- base
      for (v in values) {
        col <- paste0(v, "_", sg)
        if (!col %in% names(data)) {
          stop("Missing wide column `", col, "`.", call. = FALSE)
        }
        out[[v]] <- data[[col]][i]
      }
      out$subgroup <- sg
      out
    })
  })
  long <- dplyr::bind_rows(unlist(pieces, recursive = FALSE))
  long$subgroup <- factor(long$subgroup, levels = subgroups)
  rownames(long) <- NULL
  long
}
