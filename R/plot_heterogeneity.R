#' Heterogeneity plots
#'
#' Posterior distributions of CAMS between-trial standard deviations.
#' Bayesian CAMS uses base-R tau-marginals ([bayesmeta::plot.bayesmeta]
#' `which = 4` style) from the stored `bmr` / `bayesmeta` objects (`τ_m`,
#' `τ_γ`). `"both"` is those two `which` plots side by side (`par(mfrow)`),
#' with a shared y-axis. Frequentist CAMS uses the same two panels as
#' profile (restricted) log-likelihoods of `τ` from the stored `rma`
#' objects ([metafor::profile.rma]), also side by side with a shared
#' y-axis.
#'
#' @param x A `interactmeta_heterogeneity` object from [heterogeneity()], or a
#'   [`interactmeta`] from [cam()] (in which case [heterogeneity()] is called).
#' @param which Which components to draw: `"both"` (default; overall /
#'   combined and interaction **side by side** with a shared y-axis),
#'   `"overall"` / `"combined"` (`τ_m`), or `"interaction"` (`τ_γ`).
#' @param probs Quantiles to annotate on ggplot densities (default
#'   `c(0.5, 0.75, 0.9)`). Ignored by the bayesmeta and metafor paths.
#' @param bins Number of equal-width histogram bins for MCMC draws (default
#'   `40`). Ignored by the bayesmeta and metafor paths and when tau is a
#'   discrete support.
#' @param digits Decimal places for ggplot quantile labels (default `2`).
#' @param colours Length-2 character vector `c(overall, interaction)` for the
#'   ggplot path. Defaults to density slate tokens (not subgroup red/blue).
#' @param base_size Base font size passed to [interactmeta_theme_report()].
#' @param title Optional ggplot title (default `NULL`). Ignored by the
#'   bayesmeta and metafor paths.
#' @param xlab Optional shared x-axis title. The default (`NULL`) labels
#'   each panel `"heterogeneity τ"` (combined / `τ_m`) or
#'   `"heterogeneity τ_γ"` (interaction). A supplied value is used on
#'   ggplot panels only.
#' @param printplot If `TRUE` (default), print the ggplot. Ignored by the
#'   bayesmeta and metafor paths (always drawn).
#' @param prior For the bayesmeta path, overlay the tau prior (dashed). Same
#'   meaning as [bayesmeta::plot.bayesmeta]'s `prior` (default `FALSE`).
#' @param ... Passed to [heterogeneity()] when `x` is a fit (`prob`,
#'   `component`, …). Ignored for `interactmeta_heterogeneity` objects except that
#'   unused names are dropped quietly.
#'
#' @return For Bayesian or frequentist CAMS with stored native objects,
#'   `x` invisibly after drawing base-R plots. Otherwise a **ggplot2**
#'   object (invisibly when `printplot = TRUE`).
#'
#' @examples
#' \dontrun{
#' fit <- cam(es, reference = "Younger")
#' plot_heterogeneity(fit)
#' plot_heterogeneity(fit, which = "combined")
#' }
#'
#' @seealso [heterogeneity()], [bayesmeta::plot.bayesmeta],
#'   [metafor::profile.rma], `vignette("further")`
#' @export
plot_heterogeneity <- function(x,
                               which = c("both", "overall", "interaction", "combined"),
                               probs = c(0.5, 0.75, 0.9),
                               bins = 40L,
                               digits = 2,
                               colours = NULL,
                               base_size = interactmeta_base_size(default = 11L),
                               title = NULL,
                               xlab = NULL,
                               printplot = TRUE,
                               prior = FALSE,
                               ...) {
  which <- plot_heterogeneity_which(which)
  het <- plot_heterogeneity_resolve(x, ...)
  objs <- heterogeneity_native_tau_objects(het, which = which)
  if (!is.null(objs)) {
    heterogeneity_plot_bayesmeta(objs, prior = isTRUE(prior))
    return(invisible(het))
  }
  rma_objs <- heterogeneity_native_rma_objects(het, which = which)
  if (!is.null(rma_objs)) {
    heterogeneity_plot_metafor_profile(rma_objs)
    return(invisible(het))
  }
  p <- plot_heterogeneity_ggplot(
    het,
    which = which,
    probs = probs,
    bins = bins,
    digits = digits,
    colours = colours,
    base_size = base_size,
    title = title,
    xlab = xlab
  )
  if (isTRUE(printplot)) {
    print(p)
    return(invisible(p))
  }
  p
}

#' @rdname plot_heterogeneity
#' @method plot interactmeta_heterogeneity
#' @export
plot.interactmeta_heterogeneity <- function(x,
                                            which = c("both", "overall", "interaction", "combined"),
                                            prior = FALSE,
                                            ...) {
  plot_heterogeneity(x, which = which, prior = prior, ...)
}

#' @noRd
plot_heterogeneity_which <- function(which) {
  which <- as.character(which)[[1L]]
  if (which %in% c("combined", "combined_population")) {
    which <- "overall"
  }
  match.arg(which, c("both", "overall", "interaction"))
}

#' Native bayesmeta / bmr objects for tau plots (no refit)
#' @noRd
heterogeneity_native_tau_objects <- function(het, which = "both") {
  fit <- attr(het, "interactmeta", exact = TRUE)
  if (is.null(fit) || !identical(as.character(fit$engine %||% "")[[1L]], "bayesmeta")) {
    return(NULL)
  }
  bm <- fit$bayesmeta
  if (is.null(bm)) {
    return(NULL)
  }
  objs <- list()
  if (which %in% c("both", "overall") && !is.null(bm$overall)) {
    objs$Overall <- bm$overall
  }
  if (which %in% c("both", "interaction") && !is.null(bm$interaction)) {
    objs$Interaction <- bm$interaction
  }
  if (!length(objs)) {
    NULL
  } else {
    objs
  }
}

#' Native metafor rma objects for tau profile plots (no refit)
#' @noRd
heterogeneity_native_rma_objects <- function(het, which = "both") {
  fit <- attr(het, "interactmeta", exact = TRUE)
  if (is.null(fit) || !identical(as.character(fit$engine %||% "")[[1L]], "metafor")) {
    return(NULL)
  }
  mf <- fit$metafor
  if (is.null(mf)) {
    return(NULL)
  }
  objs <- list()
  if (which %in% c("both", "overall") && !is.null(mf$overall)) {
    objs$Overall <- mf$overall
  }
  if (which %in% c("both", "interaction") && !is.null(mf$interaction)) {
    objs$Interaction <- mf$interaction
  }
  if (!length(objs)) {
    NULL
  } else {
    objs
  }
}

#' Side-by-side bayesmeta tau-marginal plots (shared ylim)
#' @noRd
heterogeneity_plot_bayesmeta <- function(objs, prior = FALSE, ...) {
  check_bayesmeta()
  force(list(...))
  panels <- lapply(names(objs), function(nm) {
    d <- heterogeneity_bayesmeta_tau_prep(objs[[nm]], prior = prior)
    d$xlab <- heterogeneity_tau_xlab(nm)
    d
  })
  ylim <- heterogeneity_plot_bayesmeta_ylim(panels)
  n <- length(panels)
  if (n > 1L) {
    old_par <- graphics::par(mfrow = n2mfrow_traceplot(n), oma = c(0, 0, 0, 0))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  for (d in panels) {
    heterogeneity_plot_bayesmeta_tau_draw(d, ylim = ylim)
  }
  invisible(objs)
}

#' @noRd
heterogeneity_plot_bayesmeta_ylim <- function(panels) {
  ymax <- max(vapply(panels, function(d) d$maxdens, numeric(1)))
  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- 1
  }
  c(0, ymax)
}

#' Panel x-axis label: combined `τ` vs interaction `τ_γ`
#' @noRd
heterogeneity_tau_xlab <- function(nm) {
  switch(
    as.character(nm)[[1L]],
    Overall = expression("heterogeneity " * tau),
    Interaction = expression("heterogeneity " * tau[gamma]),
    expression("heterogeneity " * tau)
  )
}

#' ggplot facet labeller for panel-specific tau x-axis labels
#' @noRd
heterogeneity_tau_facet_labeller <- function(labels) {
  list(panel = lapply(as.character(labels$panel), heterogeneity_tau_xlab))
}

#' Tau-marginal panel matching [bayesmeta::plot.bayesmeta] `which = 4`
#'
#' Works for both `bayesmeta` and `bmr` objects (same `dposterior(tau=)` /
#' `qposterior(tau=)` / `$summary[ , "tau"]` interface).
#' @noRd
heterogeneity_plot_bayesmeta_tau <- function(x,
                                            xlab = NULL,
                                            prior = FALSE,
                                            ylim = NULL) {
  d <- heterogeneity_bayesmeta_tau_prep(x, prior = prior)
  if (!is.null(xlab)) {
    d$xlab <- xlab
  }
  heterogeneity_plot_bayesmeta_tau_draw(d, ylim = ylim)
}

#' @noRd
heterogeneity_bayesmeta_tau_prep <- function(x, prior = FALSE) {
  tau_hi <- tryCatch(
    as.numeric(x$qposterior(tau = 0.995))[[1L]],
    error = function(e) NA_real_
  )
  if (!is.finite(tau_hi) || tau_hi <= 0) {
    tau_hi <- 1
  }
  taurange <- c(0, tau_hi * 1.1)
  tau <- seq(
    max(c(0, taurange[[1L]] - 0.1 * diff(taurange))),
    taurange[[2L]] + 0.1 * diff(taurange),
    length.out = 200L
  )
  dens <- as.numeric(x$dposterior(tau = tau))
  maxdens <- max(dens[is.finite(dens)], na.rm = TRUE)
  if (!is.finite(maxdens) || maxdens <= 0) {
    maxdens <- 1
  }
  dens[!is.finite(dens)] <- 10 * maxdens
  lo <- as.numeric(x$summary["95% lower", "tau"])[[1L]]
  med <- as.numeric(x$summary["median", "tau"])[[1L]]
  hi <- as.numeric(x$summary["95% upper", "tau"])[[1L]]
  d_lo <- NA_real_
  d_hi <- NA_real_
  indi <- NULL
  if (is.finite(lo) && is.finite(hi) && hi > lo) {
    indi <- (tau >= lo) & (tau <= hi)
    d_lo <- as.numeric(x$dposterior(tau = lo))[[1L]]
    d_hi <- as.numeric(x$dposterior(tau = hi))[[1L]]
    if (!is.finite(d_lo)) {
      d_lo <- 10 * maxdens
    }
    if (!is.finite(d_hi)) {
      d_hi <- 10 * maxdens
    }
  }
  d_med <- NA_real_
  if (is.finite(med)) {
    d_med <- as.numeric(x$dposterior(tau = med))[[1L]]
  }
  pdens <- NULL
  if (isTRUE(prior)) {
    prior_ok <- isTRUE(x$tau.prior.proper) || is.function(x$dprior)
    if (isTRUE(prior_ok)) {
      pdens <- tryCatch(
        as.numeric(x$dprior(tau = tau)),
        error = function(e) NULL
      )
      if (!is.null(pdens) && length(pdens) != length(tau)) {
        pdens <- NULL
      }
    }
  }
  list(
    x = x,
    tau = tau,
    dens = dens,
    taurange = taurange,
    maxdens = maxdens,
    lo = lo,
    med = med,
    hi = hi,
    d_lo = d_lo,
    d_hi = d_hi,
    d_med = d_med,
    indi = indi,
    pdens = pdens,
    prior = isTRUE(prior),
    xlab = heterogeneity_tau_xlab("Overall")
  )
}

#' @noRd
heterogeneity_plot_bayesmeta_tau_draw <- function(d, ylim = NULL) {
  if (is.null(ylim)) {
    ylim <- c(0, d$maxdens)
  }
  graphics::plot(
    c(d$taurange[[1L]], d$taurange[[2L]]),
    ylim,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  graphics::polygon(
    c(0, d$tau, max(d$tau)),
    c(0, d$dens, 0),
    border = NA,
    col = grDevices::grey(0.9)
  )
  if (!is.null(d$indi) && any(d$indi)) {
    graphics::polygon(
      c(rep(d$lo, 2L), d$tau[d$indi], rep(d$hi, 2L)),
      c(0, d$d_lo, d$dens[d$indi], d$d_hi, 0),
      border = NA,
      col = grDevices::grey(0.8)
    )
  }
  if (is.finite(d$med) && is.finite(d$d_med)) {
    graphics::lines(rep(d$med, 2L), c(0, d$d_med), col = grDevices::grey(0.6))
  }
  graphics::lines(d$tau, d$dens, col = "black")
  graphics::abline(h = 0, v = 0, col = grDevices::grey(0.4))
  if (!is.null(d$pdens)) {
    graphics::lines(d$tau, d$pdens, col = "black", lty = "dashed")
  }
  xlab <- d$xlab
  if (is.null(xlab)) {
    xlab <- heterogeneity_tau_xlab("Overall")
  }
  graphics::mtext(
    side = 1,
    line = graphics::par("mgp")[[1L]],
    xlab
  )
  graphics::mtext(
    side = 2,
    line = graphics::par("mgp")[[2L]],
    expression("marginal posterior density")
  )
  graphics::axis(1)
  graphics::box()
  invisible(d$x)
}

#' Side-by-side metafor tau profile-likelihood plots (shared ylim)
#' @noRd
heterogeneity_plot_metafor_profile <- function(objs, steps = 20L) {
  panels <- lapply(names(objs), function(nm) {
    d <- heterogeneity_metafor_tau_prep(objs[[nm]], steps = steps)
    if (is.null(d)) {
      return(NULL)
    }
    d$xlab <- heterogeneity_tau_xlab(nm)
    d
  })
  panels <- Filter(Negate(is.null), panels)
  if (!length(panels)) {
    stop(
      "Cannot profile tau for this frequentist CAMS fit.",
      call. = FALSE
    )
  }
  ylim <- heterogeneity_plot_metafor_ylim(panels)
  n <- length(panels)
  if (n > 1L) {
    old_par <- graphics::par(mfrow = n2mfrow_traceplot(n), oma = c(0, 0, 0, 0))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  for (d in panels) {
    heterogeneity_plot_metafor_tau_draw(d, ylim = ylim)
  }
  invisible(objs)
}

#' @noRd
heterogeneity_plot_metafor_ylim <- function(panels) {
  lo <- min(vapply(panels, function(d) d$ylim[[1L]], numeric(1)))
  hi <- max(vapply(panels, function(d) d$ylim[[2L]], numeric(1)))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    c(-1, 0)
  } else {
    c(lo, hi)
  }
}

#' Profile (restricted) log-likelihood of tau from a stored rma fit
#' @noRd
heterogeneity_metafor_tau_prep <- function(x, steps = 20L) {
  pr <- heterogeneity_metafor_profile_tau(x, steps = steps)
  if (is.null(pr)) {
    return(NULL)
  }
  tau <- if (!is.null(pr$tau)) {
    as.numeric(pr$tau)
  } else {
    sqrt(pmax(as.numeric(pr$tau2), 0))
  }
  ll <- as.numeric(pr$ll)
  ok <- is.finite(tau) & is.finite(ll)
  tau <- tau[ok]
  ll <- ll[ok]
  if (!length(tau)) {
    return(NULL)
  }
  vc <- sqrt(max(as.numeric(x$tau2)[[1L]], 0))
  if (!is.finite(vc)) {
    vc <- NA_real_
  }
  maxll <- as.numeric(pr$maxll)[[1L]]
  if (!is.finite(maxll)) {
    maxll <- max(ll)
  }
  if (is.finite(vc) && min(tau) <= vc && max(tau) >= vc) {
    pos <- which(tau >= vc)[[1L]]
    tau <- c(tau[seq_len(pos - 1L)], vc, tau[pos:length(tau)])
    ll <- c(ll[seq_len(pos - 1L)], maxll, ll[pos:length(ll)])
  }
  xlim <- c(0, max(c(as.numeric(pr$xlim), tau, vc), na.rm = TRUE))
  if (!all(is.finite(xlim)) || xlim[[2L]] <= xlim[[1L]]) {
    xlim <- c(0, max(tau, 1))
  }
  ylim <- as.numeric(pr$ylim)
  if (length(ylim) != 2L || !all(is.finite(ylim))) {
    ylim <- range(c(maxll, ll), na.rm = TRUE) + c(-0.1, 0.1)
  }
  list(
    x = x,
    tau = tau,
    ll = ll,
    vc = vc,
    maxll = maxll,
    xlim = xlim,
    ylim = ylim,
    ylab = heterogeneity_metafor_profile_ylab(pr$method %||% x$method),
    xlab = heterogeneity_tau_xlab("Overall")
  )
}

#' @noRd
heterogeneity_metafor_profile_ylab <- function(method) {
  if (identical(as.character(method %||% "")[[1L]], "REML")) {
    "restricted log-likelihood"
  } else {
    "log-likelihood"
  }
}

#' @noRd
heterogeneity_metafor_profile_tau <- function(x, steps = 20L) {
  if (!inherits(x, "rma")) {
    return(NULL)
  }
  steps <- max(5L, as.integer(steps)[[1L]])
  withCallingHandlers(
    tryCatch(
      stats::profile(
        x,
        plot = FALSE,
        progbar = FALSE,
        sqrt = TRUE,
        steps = steps
      ),
      error = function(e) NULL
    ),
    warning = function(w) {
      if (grepl(
        "larger than the log-likelihood",
        conditionMessage(w),
        ignore.case = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' @noRd
heterogeneity_plot_metafor_tau_draw <- function(d, ylim = NULL) {
  if (is.null(ylim)) {
    ylim <- d$ylim
  }
  graphics::plot(
    d$xlim,
    ylim,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  if (is.finite(d$vc)) {
    graphics::abline(v = d$vc, lty = "dotted", col = grDevices::grey(0.6))
  }
  if (is.finite(d$maxll)) {
    graphics::abline(h = d$maxll, lty = "dotted", col = grDevices::grey(0.6))
  }
  graphics::lines(d$tau, d$ll, col = "black")
  graphics::abline(v = 0, col = grDevices::grey(0.4))
  xlab <- d$xlab
  if (is.null(xlab)) {
    xlab <- heterogeneity_tau_xlab("Overall")
  }
  graphics::mtext(
    side = 1,
    line = graphics::par("mgp")[[1L]],
    xlab
  )
  graphics::mtext(
    side = 2,
    line = graphics::par("mgp")[[2L]],
    d$ylab
  )
  graphics::axis(1)
  graphics::box()
  invisible(d$x)
}

#' Side-by-side heterogeneity ggplot used by [plot_heterogeneity()]
#' @noRd
plot_heterogeneity_ggplot <- function(het,
                                      which = c("both", "overall", "interaction", "combined"),
                                      probs = c(0.5, 0.75, 0.9),
                                      bins = 40L,
                                      digits = 2,
                                      colours = NULL,
                                      base_size = interactmeta_base_size(default = 11L),
                                      title = NULL,
                                      xlab = NULL) {
  if (!inherits(het, "interactmeta_heterogeneity")) {
    stop("`het` must be a `interactmeta_heterogeneity` object from `heterogeneity()`.", call. = FALSE)
  }
  which <- plot_heterogeneity_which(which)
  probs <- sort(unique(as.numeric(probs)))
  probs <- probs[is.finite(probs) & probs > 0 & probs < 1]
  if (!length(probs)) {
    stop("`probs` must contain values in (0, 1).", call. = FALSE)
  }
  bins <- max(5L, as.integer(bins)[[1L]])
  digits <- max(0L, as.integer(digits)[[1L]])
  colours <- plot_heterogeneity_colours(colours)
  const <- interactmeta_plot_constants()

  pieces <- plot_heterogeneity_apply_colours(
    plot_heterogeneity_pieces(het, which = which),
    colours
  )
  hist_df <- dplyr::bind_rows(lapply(pieces, function(pc) {
    plot_heterogeneity_hist(
      pc$values,
      panel = pc$panel,
      bins = bins,
      fill = pc$fill,
      weight = pc$weight
    )
  }))
  dens_df <- dplyr::bind_rows(lapply(pieces, function(pc) {
    plot_heterogeneity_density(
      pc$values,
      panel = pc$panel,
      colour = pc$line,
      weight = pc$weight
    )
  }))
  ann_df <- dplyr::bind_rows(lapply(pieces, function(pc) {
    plot_heterogeneity_annotations(
      pc$values,
      panel = pc$panel,
      probs = probs,
      digits = digits,
      colour = pc$line,
      weight = pc$weight
    )
  }))
  panel_levels <- vapply(pieces, function(pc) pc$panel, character(1))
  if (nrow(hist_df)) {
    hist_df$panel <- factor(hist_df$panel, levels = panel_levels)
  } else {
    hist_df$panel <- factor(character(), levels = panel_levels)
  }
  if (nrow(dens_df)) {
    dens_df$panel <- factor(dens_df$panel, levels = panel_levels)
  } else {
    dens_df$panel <- factor(character(), levels = panel_levels)
  }
  if (nrow(ann_df)) {
    ann_df$panel <- factor(ann_df$panel, levels = panel_levels)
  } else {
    ann_df$panel <- factor(character(), levels = panel_levels)
  }

  xmax <- max(c(hist_df$xmax, dens_df$x, ann_df$x), na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) {
    xmax <- 1
  }
  ymax <- max(c(hist_df$ymax, dens_df$y, ann_df$y), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- 1
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = const$axis_line,
      linewidth = const$line_hair
    ) +
    ggplot2::geom_rect(
      data = hist_df,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$fill
      ),
      colour = NA,
      alpha = 0.45
    ) +
    ggplot2::geom_line(
      data = dens_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$panel,
        colour = .data$colour
      ),
      linewidth = const$line_mid
    ) +
    ggplot2::geom_segment(
      data = ann_df[ann_df$vline, , drop = FALSE],
      ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = .data$y),
      colour = const$null_line,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = ann_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y_lab,
        label = .data$label,
        colour = .data$colour
      ),
      vjust = -0.35,
      size = base_size / 3.6,
      fontface = "plain"
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.06))
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::coord_cartesian(
      xlim = c(0, xmax * 1.05),
      ylim = c(0, ymax * 1.2),
      clip = "off"
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    interactmeta_theme_report(base_size = base_size, variant = "density") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 10, 8, 8, unit = "pt")
    )
  if (!is.null(xlab)) {
    p <- p + ggplot2::labs(x = xlab)
    if (length(panel_levels) > 1L) {
      p <- p + ggplot2::facet_wrap(~panel, ncol = 2L)
    }
  } else if (length(panel_levels) > 1L) {
    p <- p +
      ggplot2::facet_wrap(
        ~panel,
        ncol = 2L,
        strip.position = "bottom",
        labeller = heterogeneity_tau_facet_labeller
      ) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.placement = "outside"
      )
  } else if (length(panel_levels) == 1L) {
    p <- p + ggplot2::labs(x = heterogeneity_tau_xlab(panel_levels[[1L]]))
  }

  p
}

# --- internals --------------------------------------------------------------

#' @noRd
plot_heterogeneity_resolve <- function(x, ...) {
  if (inherits(x, "interactmeta_construction_heterogeneity")) {
    stop(
      "plot_heterogeneity() requires the CAMS information-fraction construction. ",
      "DA, AD, and SWADA heterogeneity is not the CAMS intercept/slope pair, ",
      "and a Bayesian posterior does not add CAMS information-fraction functionality.",
      call. = FALSE
    )
  }
  interactmeta_abort_non_cams(x, "plot_heterogeneity()")
  if (inherits(x, "interactmeta_heterogeneity")) {
    return(x)
  }
  if (inherits(x, "interactmeta")) {
    dots <- list(...)
    # Drop plot-only args that may arrive via plot.interactmeta_heterogeneity → ...
    drop <- c(
      "which", "probs", "bins", "digits", "colours", "base_size",
      "title", "xlab", "printplot", "prior"
    )
    dots <- dots[setdiff(names(dots), drop)]
    return(do.call(heterogeneity, c(list(fit = x), dots)))
  }
  stop(
    "`x` must be a `interactmeta_heterogeneity` from `heterogeneity()` or an `interactmeta` from `cam()`.",
    call. = FALSE
  )
}

#' @noRd
plot_heterogeneity_colours <- function(colours = NULL) {
  const <- interactmeta_plot_constants()
  default <- c(
    overall = const$density_fill,
    interaction = const$density_point
  )
  if (is.null(colours)) {
    return(default)
  }
  colours <- as.character(colours)
  if (length(colours) == 1L) {
    colours <- rep(colours, 2L)
  }
  if (length(colours) < 2L) {
    stop("`colours` must be length 2 (overall, interaction).", call. = FALSE)
  }
  stats::setNames(colours[1:2], c("overall", "interaction"))
}

#' @noRd
plot_heterogeneity_prior_label <- function(type) {
  if (identical(type, "intercept")) {
    "Overall"
  } else {
    "Interaction"
  }
}

#' @noRd
plot_heterogeneity_pieces <- function(het, which = "both") {
  draws <- het$draws
  if (is.null(draws) || !nrow(draws)) {
    stop("`het` has no posterior draws.", call. = FALSE)
  }
  types <- unique(as.character(draws$type))
  want_int <- which %in% c("both", "overall")
  want_slp <- which %in% c("both", "interaction")
  has_int <- want_int && "intercept" %in% types
  has_slp <- want_slp && "slope" %in% types

  pieces <- list()
  if (has_int) {
    ii <- draws$type == "intercept"
    vals <- draws$value[ii]
    pieces[[length(pieces) + 1L]] <- list(
      panel = "Overall",
      type = "intercept",
      values = vals,
      weight = plot_heterogeneity_piece_weight(draws, ii, length(vals)),
      fill = NA_character_,
      line = NA_character_,
      title = plot_heterogeneity_prior_label("intercept")
    )
  }
  if (has_slp) {
    sl <- draws[draws$type == "slope", , drop = FALSE]
    # One interaction-block parameter (first if several).
    param <- unique(as.character(sl$parameter))[[1L]]
    keep <- sl$parameter == param
    vals <- sl$value[keep]
    pieces[[length(pieces) + 1L]] <- list(
      panel = "Interaction",
      type = "slope",
      values = vals,
      weight = plot_heterogeneity_piece_weight(sl, keep, length(vals)),
      fill = NA_character_,
      line = NA_character_,
      title = plot_heterogeneity_prior_label("slope")
    )
  }
  if (!length(pieces)) {
    stop(
      "No matching heterogeneity draws for `which = \"", which, "\"`.",
      call. = FALSE
    )
  }
  pieces
}

#' @noRd
plot_heterogeneity_apply_colours <- function(pieces, colours) {
  colours <- plot_heterogeneity_colours(colours)
  lapply(pieces, function(pc) {
    key <- if (identical(pc$type, "intercept")) "overall" else "interaction"
    pc$line <- colours[[key]]
    pc$fill <- colours[[key]]
    pc
  })
}

#' @noRd
plot_heterogeneity_piece_weight <- function(df, i, n) {
  if (!is.data.frame(df) || !".weight" %in% names(df)) {
    return(NULL)
  }
  w <- as.numeric(df$.weight[i])
  if (length(w) != n || !any(is.finite(w) & w > 0)) {
    return(NULL)
  }
  w
}

#' @noRd
plot_heterogeneity_normalize_weight <- function(values, weight) {
  n <- length(values)
  if (is.null(weight) || length(weight) != n) {
    return(rep(1 / max(n, 1L), n))
  }
  w <- pmax(as.numeric(weight), 0)
  w[!is.finite(w)] <- 0
  if (!sum(w)) {
    return(rep(1 / max(n, 1L), n))
  }
  w / sum(w)
}

#' Discrete tau support (quadrature nodes) vs many unique MCMC draws
#' @noRd
plot_heterogeneity_is_discrete <- function(values, bins = 40L) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  n <- length(values)
  if (n < 2L) {
    return(FALSE)
  }
  n_u <- length(unique(values))
  n_u <= max(2L * as.integer(bins)[[1L]], 80L) && n_u < 0.5 * n
}

#' @noRd
plot_heterogeneity_hist <- function(values,
                                    panel = "Overall",
                                    bins = 40L,
                                    fill = "grey",
                                    weight = NULL) {
  empty <- data.frame(
    xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric(),
    panel = character(), fill = character(), stringsAsFactors = FALSE
  )
  values <- as.numeric(values)
  ok <- is.finite(values) & values >= 0
  values <- values[ok]
  if (!length(values)) {
    return(empty)
  }
  weight <- plot_heterogeneity_normalize_weight(values, if (is.null(weight)) NULL else weight[ok])
  if (length(unique(values)) == 1L) {
    return(empty)
  }
  if (plot_heterogeneity_is_discrete(values, bins = bins)) {
    agg <- stats::aggregate(list(w = weight), by = list(x = values), FUN = sum)
    o <- order(agg$x)
    u <- agg$x[o]
    w <- agg$w[o] / sum(agg$w)
    edges <- plot_heterogeneity_node_edges(u)
    width <- pmax(edges$right - edges$left, .Machine$double.eps)
    y <- w / width
    return(data.frame(
      xmin = edges$left,
      xmax = edges$right,
      ymin = 0,
      ymax = y,
      panel = panel,
      fill = fill,
      stringsAsFactors = FALSE
    ))
  }
  br <- graphics::hist(values, breaks = bins, plot = FALSE)
  data.frame(
    xmin = br$breaks[-length(br$breaks)],
    xmax = br$breaks[-1L],
    ymin = 0,
    ymax = br$density,
    panel = panel,
    fill = fill,
    stringsAsFactors = FALSE
  )
}

#' Voronoi edges on [0, Inf) for unique support nodes
#' @noRd
plot_heterogeneity_node_edges <- function(u) {
  u <- sort(unique(as.numeric(u)))
  if (length(u) == 1L) {
    half <- max(u[[1L]] / 2, 0.02)
    return(list(left = max(0, u[[1L]] - half), right = u[[1L]] + half))
  }
  mid <- (u[-length(u)] + u[-1L]) / 2
  left0 <- max(0, u[[1L]] - (mid[[1L]] - u[[1L]]))
  rightk <- u[[length(u)]] + (u[[length(u)]] - mid[[length(mid)]])
  list(left = c(left0, mid), right = c(mid, rightk))
}

#' @noRd
plot_heterogeneity_safe_density <- function(values, n = 512L, weight = NULL) {
  values <- as.numeric(values)
  ok <- is.finite(values) & values >= 0
  values <- values[ok]
  if (length(values) < 2L || isTRUE(stats::sd(values) <= 0)) {
    return(NULL)
  }
  w <- NULL
  if (!is.null(weight)) {
    w <- plot_heterogeneity_normalize_weight(values, weight[ok])
    if (length(unique(w)) <= 1L) {
      w <- NULL
    }
  }
  tryCatch(
    if (is.null(w)) {
      stats::density(values, from = 0, n = as.integer(n)[[1L]])
    } else {
      stats::density(values, from = 0, n = as.integer(n)[[1L]], weights = w)
    },
    error = function(e) NULL
  )
}

#' @noRd
plot_heterogeneity_density <- function(values,
                                       panel = "Overall",
                                       colour = "black",
                                       n = 512L,
                                       weight = NULL) {
  d <- plot_heterogeneity_safe_density(values, n = n, weight = weight)
  if (is.null(d)) {
    return(data.frame(
      x = numeric(), y = numeric(), panel = character(), colour = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    x = d$x,
    y = d$y,
    panel = panel,
    colour = colour,
    stringsAsFactors = FALSE
  )
}

#' @noRd
plot_heterogeneity_weighted_quantile <- function(values, weight, probs) {
  o <- order(values)
  x <- values[o]
  w <- plot_heterogeneity_normalize_weight(x, weight[o])
  cw <- cumsum(w)
  vapply(probs, function(p) {
    hit <- which(cw >= p)[1L]
    if (!length(hit) || is.na(hit)) {
      x[[length(x)]]
    } else {
      x[[hit]]
    }
  }, numeric(1))
}

#' @noRd
plot_heterogeneity_annotations <- function(values,
                                           panel = "Overall",
                                           probs = c(0.5, 0.75, 0.9),
                                           digits = 2,
                                           colour = "black",
                                           weight = NULL) {
  values <- as.numeric(values)
  ok <- is.finite(values) & values >= 0
  values <- values[ok]
  empty <- data.frame(
    x = numeric(), y = numeric(), y_lab = numeric(), label = character(),
    panel = character(), colour = character(), vline = logical(),
    stringsAsFactors = FALSE
  )
  if (!length(values)) {
    return(empty)
  }
  if (!is.null(weight)) {
    weight <- weight[ok]
  }
  qs <- if (!is.null(weight) && length(weight) == length(values)) {
    plot_heterogeneity_weighted_quantile(values, weight, probs)
  } else {
    stats::quantile(values, probs = probs, names = FALSE, type = 7)
  }
  keep <- !duplicated(round(qs, digits + 4L))
  qs <- qs[keep]
  if (!length(qs)) {
    return(empty)
  }
  d <- plot_heterogeneity_safe_density(values, n = 512L, weight = weight)
  y_at <- if (is.null(d)) {
    rep(1, length(qs))
  } else {
    vapply(qs, function(q) {
      stats::approx(d$x, d$y, xout = q, rule = 2)$y
    }, numeric(1))
  }
  data.frame(
    x = as.numeric(qs),
    y = y_at,
    y_lab = y_at,
    label = format(round(qs, digits), nsmall = digits, trim = TRUE),
    panel = panel,
    colour = colour,
    vline = c(TRUE, rep(FALSE, max(0L, length(qs) - 1L))),
    stringsAsFactors = FALSE
  )
}
