#' Format an interactmeta effect object according to `format`
#' @noRd
comparison_format_effect_object <- function(obj,
                                      format = c("list", "summary"),
                                      type,
                                      fit = NULL) {
  format <- match.arg(format)
  if (identical(format, "list")) {
    return(obj)
  }
  if (identical(format, "summary")) {
    if (is.null(obj$summary)) {
      stop(
        "No summary available; call with `summary = TRUE` or use another format.",
        call. = FALSE
      )
    }
    return(tibble::as_tibble(obj$summary))
  }
  stop("Unknown format.", call. = FALSE)
}
