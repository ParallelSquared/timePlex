################################################################################
# Extended Data Fig. 6a-c: timePlex vs LF no-plex ratios against TIC and RT
################################################################################

quant_cols_default <- c("plex_Area", "coeff", "plex_Area_rtnorm", "coeff_rtnorm")

read_jmod_ids <- function(results_dir, time_channel = NULL, q_cut = 0.01,
                          quant_cols=quant_cols_default) {
  pq <- file.path(results_dir, "filtered_IDs.parquet")                                        # JMod-timePlex_beta
  if (!file.exists(pq)) pq <- file.path(results_dir, "filtered_IDs_parquet_columns.parquet")  # older JMod
  ids <- arrow::read_parquet(
    pq,
    col_select = dplyr::any_of(c("untag_prec", "time_channel", "is_decoy",
                                 "Qvalue", "rt", "plex_Area", "coeff", "plex_Area_rtnorm", "coeff_rtnorm"))
  )
  if (!is.null(time_channel)) {
    ids <- filter(ids, .data$time_channel == .env$time_channel)
  }
  ids %>%
    dplyr::filter(!is_decoy, Qvalue <= q_cut) %>%
    dplyr::arrange(Qvalue) %>%
    dplyr::distinct(untag_prec, .keep_all = TRUE) %>%
    dplyr::select(untag_prec, rt, all_of(quant_cols))
}

## MS1 TIC trace from an mzML: rt (min), tic.
read_ms1_tic <- function(mzml) {
  ms <- mzR::openMSfile(mzml)
  on.exit(mzR::close(ms))
  hdr <- mzR::header(ms)
  hdr <- hdr[hdr$msLevel == 1L, ]
  tibble::tibble(rt = hdr$retentionTime / 60, tic = hdr$totIonCurrent) %>%
    arrange(rt)
}

## TIC interpolated at the given RTs (min); clamps to the edge outside the trace.
tic_at_rt <- function(tic, rt) {
  approx(tic$rt, tic$tic, xout = rt, rule = 2)$y
}

## Precursors seen in both runs, with each run's TIC at its own precursor RT and
## log2(no-plex / timePlex) for every quant column (NA unless both sides > min_quant).
pair_tp_noplex <- function(tp_ids, noplex_ids, tic_tp, tic_noplex,
                           quant_cols = quant_cols_default, min_quant = 1) {
  paired <- dplyr::inner_join(tp_ids, noplex_ids, by = "untag_prec",
                       suffix = c("_tp", "_noplex")) %>%
    dplyr::mutate(
      tic_tp         = tic_at_rt(tic_tp, rt_tp),
      tic_noplex     = tic_at_rt(tic_noplex, rt_noplex),
      log2_tic_ratio = log2(tic_noplex / tic_tp)
    )
  for (q in quant_cols) {
    num <- paired[[paste0(q, "_noplex")]]
    den <- paired[[paste0(q, "_tp")]]
    paired[[paste0("log2_", q)]] <- ifelse(num > min_quant & den > min_quant, log2(num / den), NA_real_)
  }
  paired
}


## ---- plots ------------------------------------------------------------------

ratio_ylab <- "Log2, ratio (no-plex / timePlex)"

sd_label <- function(x) {
  x <- x[is.finite(x)]
  sprintf("precursor ratio sd = %.3f  (n = %s)", sd(x), format(length(x), big.mark = ","))
}

ratio_points <- function(colour) paper_raster(geom_point(colour = colour, alpha = 0.2, size = 0.5, stroke = 0))
ratio_zero   <- function() geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.25)

## a/b. tic = "noplex": TIC of the no-plex run at its precursor RT; "tp": TIC of the
## timePlex run at its precursor RT. The x axis zooms to the central `pct` of the
## points (all points are kept; Pearson r of log10 TIC vs log2 ratio uses all).
plot_ratio_vs_tic <- function(paired, quant = "plex_Area", tic = "noplex", pct = 0.99,
                              ylim = c(-4, 4)) {
  tic_col <- paste0("tic_", tic)
  run_lab <- c(noplex = "no-plex", tp = "timePlex")[[tic]]
  y_col   <- paste0("log2_", quant)
  d <- dplyr::filter(paired, is.finite(.data[[y_col]]), .data[[tic_col]] > 0)
  probs <- c((1 - pct) / 2, 1 - (1 - pct) / 2)
  r <- cor(log10(d[[tic_col]]), d[[y_col]])
  ggplot(d, aes(x = .data[[tic_col]], y = .data[[y_col]])) +
    ratio_points("steelblue") +
    ratio_zero() +
    scale_x_log10(breaks = scales::breaks_log(6)) +
    coord_cartesian(xlim = quantile(d[[tic_col]], probs), ylim = ylim) +
    labs(title = sprintf("Precursor ratio vs %s TIC", run_lab),
         subtitle = sprintf("Pearson r = %.3f  (n = %s)", r, format(nrow(d), big.mark = ",")),
         x = sprintf("Log10, MS1 TIC at %s precursor RT", run_lab), y = ratio_ylab) +
    paper_theme_box(grid = TRUE) +
    theme(plot.title = element_text(hjust = 0))
}

## c
plot_ratio_tic_vs_rt <- function(paired, quant = "plex_Area", ylim = c(-4, 4)) {
  d <- dplyr::filter(paired, is.finite(.data[[paste0("log2_", quant)]]))
  long <- bind_rows(
    tibble::tibble(rt = d$rt_tp, value = d[[paste0("log2_", quant)]], series = "Precursor ratio"),
    tibble::tibble(rt = d$rt_tp, value = d$log2_tic_ratio, series = "TIC ratio")
  )
  ggplot(long, aes(x = rt, y = value, colour = series)) +
    paper_raster(geom_point(alpha = 0.2, size = 0.5, stroke = 0)) +
    ratio_zero() +
    scale_colour_manual(values = c("Precursor ratio" = "steelblue", "TIC ratio" = "darkorange")) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2, shape = 15))) +
    coord_cartesian(ylim = ylim) +
    labs(title = "Covariation of TIC ratio and precursor ratio",
         subtitle = sd_label(d[[paste0("log2_", quant)]]),
         x = "Retention time (min)", y = ratio_ylab, colour = NULL) +
    paper_theme_box(grid = TRUE) +
    theme(plot.title = element_text(hjust = 0),
          legend.position = "inside", legend.position.inside = c(0.02, 0.98),
          legend.justification = c(0, 1), legend.background = element_blank())
}

## d-g
plot_ratio_vs_rt <- function(paired, quant, title = quant, colour = "steelblue",
                             ylim = c(-4, 4)) {
  d <- dplyr::filter(paired, is.finite(.data[[paste0("log2_", quant)]]))
  ggplot(d, aes(x = rt_tp, y = .data[[paste0("log2_", quant)]])) +
    ratio_points(colour) +
    ratio_zero() +
    coord_cartesian(ylim = ylim) +
    labs(title = title, subtitle = sd_label(d[[paste0("log2_", quant)]]),
         x = "Retention time (min)", y = ratio_ylab) +
    paper_theme_box(grid = TRUE) +
    theme(plot.title = element_text(hjust = 0))
}

## Extended Data Fig. 6 a-g: a / b = ratio vs no-plex / timePlex TIC, c spans the full
## width. same_precursors = TRUE restricts every panel to precursors with a
## finite ratio for all four quant columns, so the sds are compared on one set.
plot_tp_noplex_figure <- function(paired, same_precursors = TRUE, ylim = c(-4, 4)) {
  paper_figure("ED6")
  if (same_precursors) {
    paired <- dplyr::filter(paired, if_all(paste0("log2_", quant_cols_default), is.finite))
  }
  patchwork::wrap_plots(
    plot_ratio_vs_tic(paired, "plex_Area", tic = "noplex", ylim = ylim),
    plot_ratio_vs_tic(paired, "plex_Area", tic = "tp", ylim = ylim),
    plot_ratio_tic_vs_rt(paired, "plex_Area", ylim = ylim),
    plot_ratio_vs_rt(paired, "plex_Area", "MS1, before normalization", "steelblue", ylim),
    plot_ratio_vs_rt(paired, "plex_Area_rtnorm", "MS1, after normalization", "firebrick", ylim),
    plot_ratio_vs_rt(paired, "coeff", "MS2, before normalization", "steelblue", ylim),
    plot_ratio_vs_rt(paired, "coeff_rtnorm", "MS2, after normalization", "firebrick", ylim),
    design = "AB\nCC\nDE\nFG"
  ) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = paper_size("tag"), face = "bold"),
          plot.title = element_text(size = 8, hjust = 0),
          plot.subtitle = element_text(size = 7))
}

## The assembled figure plus every panel on its own (a-g, no letters)
save_tp_noplex_figure <- function(fig, path, width = 7.09, height = 9, out_dir = dirname(path)) {
  paper_save(fig, path, width = width, height = height)
  panel_files <- file.path(out_dir, paste0("ExtDataFig6", letters[1:7], ".pdf"))
  panel_w <- c(3.5, 3.5, 7, 3.5, 3.5, 3.5, 3.5)
  for (i in seq_along(panel_files)) paper_save(fig[[i]], panel_files[i], width = panel_w[i], height = 2.2)
  c(path, panel_files)
}
