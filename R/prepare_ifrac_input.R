is_ifrac_output <- function(x) {
  if (!is.data.frame(x)) {
    return(FALSE)
  }
  cols <- grep("^ifrac_", names(x), value = TRUE)
  cols <- cols[!grepl("^ifrac2_", cols)]
  length(cols) > 0L
}

prepare_ifrac_input <- function(x,
                                reference = NULL,
                                trial = "trial",
                                subgroup = "subgroup",
                                contrast = "all",
                                contrast_sep = "-",
                                big_var = 100^2,
                                include_overall = FALSE,
                                include_interaction = FALSE) {
  if (is_ifrac_output(x)) {
    out <- tibble::as_tibble(x)
    # Already designed: reference must agree with the stored design (ifrac_<nonref>).
    ref <- interactmeta_reference(out, reference = reference, subgroup_var = subgroup)
    out <- interactmeta_ensure_reference_column(out, ref)
    # Reject ifrac_<reference> columns (pi must be non-reference only).
    ifrac_cols <- grep("^ifrac_", names(out), value = TRUE)
    ifrac_cols <- ifrac_cols[!grepl("^ifrac2_", ifrac_cols)]
    bad <- ifrac_cols[sub("^ifrac_", "", ifrac_cols) == ref]
    if (length(bad)) {
      stop(
        "Design column(s) ",
        paste(bad, collapse = ", "),
        " name the reference subgroup `",
        ref,
        "`. CAMS `ifrac_*` columns must be non-reference only ",
        "(pi = IF_nonref, reference weight = 1 - pi).",
        call. = FALSE
      )
    }
    level_col <- if ("subgroup" %in% names(out)) "subgroup" else subgroup
    return(add_ifrac2_design(out, subgroup_col = level_col))
  }

  if (is.null(reference)) {
    stop(
      "`reference` is required when `x` is not already output from `ifrac()`.",
      call. = FALSE
    )
  }

  designed <- ifrac(
    x,
    trial = trial,
    subgroup = subgroup,
    reference = reference,
    contrast = contrast,
    contrast_sep = contrast_sep,
    big_var = big_var,
    include_overall = include_overall,
    include_interaction = include_interaction
  )
  level_col <- if ("subgroup" %in% names(designed)) "subgroup" else subgroup
  add_ifrac2_design(designed, subgroup_col = level_col)
}
