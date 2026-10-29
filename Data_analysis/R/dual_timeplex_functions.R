################################################################################
# Fig. 6: dual-timePlex at >1,200 samples/day
################################################################################


## ---- Axis-gap compression shared by the QC panels ---------------------------
## One squash() mapping, computed once from the run time intervals, is applied
## to everything drawn (bars, bins, axis breaks, slash marks), so positions and
## labels cannot disagree. Every between-run gap wider than gap_min_hr is
## compressed.

#' Build the squash mapping from per-run time intervals [s, e] (hours).
#' Returns identity when nothing needs truncating.
DualTP_gap_squash <- function(s, e, truncate_gap = TRUE, gap_min_hr = 1,
                              gap_frac = 0.05, n_breaks = 5, gap_pad_hr = 0) {
  o <- order(s); s <- s[o]; e <- pmax(cummax(e[o]), s)
  gaps <- if (length(s) > 1) s[-1] - head(e, -1) else numeric(0)
  big  <- which(gaps > gap_min_hr)
  if (!truncate_gap || !length(big))
    return(list(truncate = FALSE, sq = identity,
                block_of = function(x) rep(1L, length(x))))
  
  bstart <- s[c(1L, big + 1L)]                 # block boundaries in TRUE hours
  bend   <- e[c(big, length(e))]
  kept_span <- sum(bend - bstart)
  GAP <- gap_frac * max(kept_span, 1e-9) + gap_pad_hr
  off <- cumsum(c(0, bstart[-1] - head(bend, -1) - GAP))
  
  block_of <- function(x) pmax(1L, findInterval(x, bstart))
  sq <- function(x) x - off[block_of(x)]
  
  blk_breaks <- function(lo, hi, n) {
    b <- scales::pretty_breaks(n)(c(lo, hi)); b[b >= lo & b <= hi]
  }
  tb <- unlist(mapply(blk_breaks, bstart, bend,
                      MoreArgs = list(n = n_breaks), SIMPLIFY = FALSE))
  step <- min(diff(sort(unique(tb))))
  acc  <- 10^floor(log10(step))
  
  list(truncate = TRUE, sq = sq, block_of = block_of, GAP = GAP,
       kept_span = kept_span,
       brk = sq(tb), lbl = scales::number(tb, accuracy = acc),
       x_break = sq(head(bend, -1)) + GAP / 2,      # junction centres, squashed
       sub = sprintf("x axis truncated: %.1f h of idle time removed across %d gap%s",
                     sum(gaps[big]) - length(big) * GAP, length(big),
                     ifelse(length(big) > 1, "s", "")))
}

## slash-mark data.frame for every compressed junction
.DualTP_slashes <- function(g, ymax) {
  do.call(rbind, lapply(g$x_break, function(xb)
    data.frame(xc = xb + c(-1, 1) * g$GAP * 0.16, dx = g$GAP * 0.30,
               y0 = 0.005 * ymax, y1 = 0.07 * ymax)))
}

## shared coloring: hue = plexDIA channel (CARTO "Bold"), lightness = timePlex
## column (earliest darkest). Returns list(vals, key_of).
.DualTP_fill_scale <- function(pd_chans, tp_chans) {
  bold <- c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74",
            "#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")
  if (length(pd_chans) > length(bold))
    stop("More plexDIA channels (", length(pd_chans), ") than Bold palette colors")
  base_cols <- setNames(bold[seq_along(pd_chans)], pd_chans)
  shade <- function(col, f) {
    to <- if (f < 0) "black" else "white"
    colorRampPalette(c(col, to))(101)[round(abs(f) * 100) + 1]
  }
  fs <- if (length(tp_chans) == 1) 0 else seq(-0.3, 0.5, length.out = length(tp_chans))
  combos <- expand.grid(time_channel = tp_chans, plex_channel = pd_chans)
  combos <- combos[order(combos$plex_channel, combos$time_channel), ]
  combos$key <- paste0("d", combos$plex_channel, " tp", combos$time_channel)
  combos$col <- mapply(function(p, tc) shade(base_cols[as.character(p)],
                                             fs[match(tc, tp_chans)]),
                       combos$plex_channel, combos$time_channel)
  list(vals = setNames(combos$col, combos$key), levels = combos$key)
}


## ---- QC multipanel (Fig. 6c-f) ----------------------------------------------
## c RT histogram, d RT offset vs the middle column, e proteins per run, f
## observed vs expected fold change per timePlex channel. c-e share one x axis
## per block of consecutive runs (real acquisition hours); "..." marks the
## idle time between blocks.

DualTP_QC_multipanel <- function(dat, runs = NULL, binwidth = 0.5,
                                 quant_col = "coeff_ChannelUnique_Area_rtnorm",
                                 ## panel d: per-sample yeast log2 offsets (any
                                 ## common scale); a pair's expected FC is
                                 ## offset[num] - offset[den]
                                 sample_exp = c(A = -2, B = 0, C = -3),
                                 ## cross-sample pairs, oriented POSITIVE
                                 ## (A/C = +1, B/A = +2, B/C = +3); like-pairs
                                 ## (A/A, B/B, C/C, expected 0) added automatically
                                 dens_pairs = data.frame(
                                   numerator   = c("A", "B", "B"),
                                   denominator = c("C", "A", "C")),
                                 dens_grouper = "seqcharge",
                                 middle_channel = 1,
                                 truncate_gap = TRUE, gap_min_hr = 1,
                                 gap_frac = 0.05, gap_extra_hr = 0.25,
                                 n_breaks = 5,
                                 set_colors = c("timePlex A" = "#FF4040",
                                                "timePlex B" = "#7143FF"),
                                 tp_greens = c("#1d3319", "#437a3b", "#7cc969"),
                                 tp_labels = c("\u03940 min", "\u03942.5 min", "\u03945 min"),
                                 extra = "DualTP", out_dir = "Figs/Fig6") {
  paper_figure(out_dir)
  if (!is.null(runs)) dat <- dat %>% dplyr::filter(Run %in% runs)

  ## ---- shared x geometry ----------------------------------------------------
  dat$t_hr <- dat$t_start_hr + dat$rt / 60
  run_iv <- dat %>% dplyr::group_by(Run, plex, injection_index) %>%
    dplyr::summarise(s = min(t_hr, na.rm = TRUE), e = max(t_hr, na.rm = TRUE),
                     .groups = "drop") %>% dplyr::arrange(injection_index)
  g <- DualTP_gap_squash(run_iv$s, run_iv$e, truncate_gap, gap_min_hr, gap_frac,
                         n_breaks, gap_pad_hr = gap_extra_hr)
  dat$blk <- g$block_of(dat$t_hr)
  run_iv$blk     <- g$block_of(run_iv$s)
  run_iv$center  <- (run_iv$s + run_iv$e) / 2
  run_iv$run_lab <- paste("Run", run_iv$injection_index)
  run_iv$set_lab <- sub("^timePlex ", "Set ", run_iv$plex)
  n_blk <- length(unique(run_iv$blk))
  bar_w <- if (nrow(run_iv) > 1) 0.8 * min(diff(sort(run_iv$center))) else 0.02
  # every time panel spans the same hours per block, so c-e line up
  blk_lim <- run_iv %>% dplyr::group_by(blk) %>%
    dplyr::summarise(lo = min(s) - bar_w / 2, hi = max(e) + bar_w / 2, .groups = "drop") %>%
    tidyr::pivot_longer(c(lo, hi), values_to = "x")

  fills <- .DualTP_fill_scale(sort(unique(dat$plex_channel)),
                              tp_chans <- sort(unique(dat$time_channel)))
  dat$fill_key <- factor(paste0("d", dat$plex_channel, " tp", dat$time_channel),
                         levels = fills$levels)

  ## green scale for timePlex columns: darkest = first, lightest = last
  greens <- if (length(tp_chans) == length(tp_greens)) tp_greens
  else colorRampPalette(c(tp_greens[1], tp_greens[length(tp_greens)]))(length(tp_chans))
  names(greens) <- tp_chans
  tp_labs <- setNames(tp_labels[seq_along(tp_chans)], tp_chans)

  ## one facet per block of runs; strips hidden, "..." drawn in the gap
  gap_in <- 0.35
  time_facets <- list(
    facet_grid(. ~ blk, scales = "free_x", space = "free_x"),
    geom_blank(data = blk_lim, aes(x = x), inherit.aes = FALSE),
    coord_cartesian(clip = "off")
  )
  time_theme <- theme(strip.text = element_blank(), strip.background = element_blank(),
                      panel.spacing.x = unit(gap_in, "in"))
  gap_dots <- function(y) {
    if (n_blk < 2) return(NULL)
    geom_text(data = data.frame(blk = 1, y = y), aes(x = Inf, y = y), label = "...",
              inherit.aes = FALSE, hjust = -1.1, fontface = "bold", size = paper_geom_text(11))
  }
  hours <- scale_x_continuous(breaks = scales::breaks_width(0.5), expand = expansion(mult = 0.01))

  ## ---- c: RT histogram, continuous acquisition time -------------------------
  ymax_a <- max(table(dat$blk, floor(dat$t_hr / (binwidth / 60))))
  pa <- ggplot(dat, aes(x = t_hr, fill = fill_key)) +
    geom_histogram(binwidth = binwidth / 60, position = "stack", color = NA) +
    geom_text(data = run_iv, aes(x = center, y = Inf, label = run_lab), inherit.aes = FALSE,
              vjust = -1.45, size = paper_geom_text(8)) +
    geom_text(data = run_iv, aes(x = center, y = Inf, label = set_lab, colour = plex), inherit.aes = FALSE,
              vjust = -0.25, size = paper_geom_text(8)) +
    gap_dots(ymax_a / 2) +
    scale_colour_manual(values = set_colors, guide = "none") +
    scale_fill_manual(values = fills$vals, drop = FALSE, guide = "none") +
    hours +
    scale_y_continuous(labels = function(x) comma(x / binwidth),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Acquisition time (hrs)", y = "Precursors / min") +
    time_facets + paper_theme() + time_theme +
    theme(plot.margin = margin(22, 3, 2, 3))

  ## ---- d: precursor-level RT offset vs the middle column ---------------------
  ## one dot per (Run, timePlex column, seqcharge): x = median empirical RT on
  ## the clock; y = difference in empirical rt between this column and the
  ## MIDDLE column for the same seqcharge (middle column's own dots sit at 0).
  mid_b <- dat %>% dplyr::filter(time_channel == middle_channel) %>%
    dplyr::group_by(Run, seqcharge) %>%
    dplyr::summarise(rt_mid = median(rt, na.rm = TRUE), .groups = "drop")
  rt_pts <- dat %>% dplyr::inner_join(mid_b, by = c("Run", "seqcharge")) %>%
    dplyr::group_by(Run, plex, time_channel, seqcharge) %>%
    dplyr::summarise(x = median(t_hr, na.rm = TRUE),
                     d_rt = median(rt - rt_mid, na.rm = TRUE),
                     .groups = "drop")
  rt_pts$blk <- g$block_of(rt_pts$x)

  rng_b <- range(rt_pts$d_rt, na.rm = TRUE)
  pb <- ggplot(rt_pts, aes(x = x, y = d_rt, colour = factor(time_channel))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.25) +
    paper_raster(geom_point(size = 0.8, alpha = 0.015, stroke = 0)) +   # all ~500k precursor points
    gap_dots(mean(rng_b)) +
    scale_colour_manual(values = greens, guide = "none") +
    hours +
    scale_y_continuous(breaks = paper_int_breaks(3)) +
    labs(x = "Acquisition time (hrs)", y = "\u0394RT from middle\ntimePlex channel") +
    time_facets + paper_theme() + time_theme +
    theme(plot.margin = margin(2, 3, 2, 3))

  ## ---- e: stacked protein bars per run --------------------------------------
  counts <- dat %>% dplyr::filter(Protein_Qvalue < 0.01) %>%
    dplyr::distinct(Run, plex, time_channel, plex_channel, Sample, prot, fill_key) %>%
    dplyr::group_by(Run, plex, time_channel, plex_channel, Sample, fill_key) %>%
    dplyr::count() %>% dplyr::ungroup() %>%
    dplyr::left_join(run_iv[, c("Run", "center", "blk")], by = "Run")
  totals_c <- counts %>% dplyr::group_by(Run, plex, blk, center) %>%
    dplyr::summarise(n = sum(n), .groups = "drop")

  pc <- ggplot(counts, aes(x = center, y = n, fill = fill_key)) +
    geom_col(width = bar_w, position = "stack", color = NA) +
    gap_dots(max(totals_c$n) / 2) +
    scale_fill_manual(values = fills$vals, drop = FALSE, guide = "none") +
    scale_x_continuous(breaks = run_iv$center, labels = run_iv$run_lab,
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Protein data points") +
    time_facets + paper_theme() + time_theme +
    theme(axis.text.x = element_text(size = 8), axis.ticks.x = element_blank(),
          plot.margin = margin(2, 3, 2, 3))

  ## ---- f: per-timechannel slope plots from CHANNEL-PAIRWISE ratios -----------
  ## WITHIN each (run x timePlex column), every ratio is between two individual
  ## channels: cross-sample pairs (dens_pairs) use ALL numerator-channel x
  ## denominator-channel combinations; like-pairs (A/A, B/B, C/C, expected 0)
  ## use all combinations of a sample's DIFFERENT channels -- so the IQR at
  ## every tick measures the same channel-level spread. Human precursors are
  ## used ONLY to center each (run x column x pair) ratio distribution and are
  ## then dropped; the curve is yeast medians per expected tick.
  dd <- DualTP_assign_species(dat)
  vals <- dd[[quant_col]]
  vals[vals <= 1] <- NA                  # quant <= 1 treated as missing
  dd$quant_val <- log2(vals)
  dd <- dd[is.finite(dd$quant_val), ]
  dd <- dd %>%
    dplyr::group_by(Run, time_channel, plex_channel) %>%
    dplyr::mutate(quant_val = quant_val - median(quant_val, na.rm = TRUE)) %>%
    dplyr::ungroup()
  dd <- dd %>% dplyr::distinct(Run, time_channel, plex_channel,
                               .data[[dens_grouper]], .keep_all = TRUE)

  ## all channel x channel combinations within (run x column x precursor)
  sdt <- as.data.table(dd)[, c("Run", "plex", "time_channel", "Sample", "species",
                               "plex_channel", "quant_val", dens_grouper), with = FALSE]
  sp <- merge(sdt, sdt,
              by = c("Run", "plex", "time_channel", "species", dens_grouper),
              allow.cartesian = TRUE, suffixes = c(".x", ".y"))

  rat <- list()
  for (i in seq_len(nrow(dens_pairs))) {
    a <- dens_pairs$numerator[i]; b <- dens_pairs$denominator[i]
    r <- sp[Sample.x == a & Sample.y == b]
    if (!nrow(r)) next
    rat[[paste0(a, b)]] <- data.frame(
      Run = r$Run, plex = r$plex, time_channel = r$time_channel,
      species = r$species, ratio = r$quant_val.x - r$quant_val.y,
      pair = paste0(a, "/", b),
      exp_fc = unname(sample_exp[a] - sample_exp[b]))
  }
  ## like-pair nulls: a sample's DIFFERENT channels against each other
  lk <- sp[Sample.x == Sample.y & plex_channel.x < plex_channel.y]
  like <- data.frame(Run = lk$Run, plex = lk$plex, time_channel = lk$time_channel,
                     species = lk$species, ratio = lk$quant_val.x - lk$quant_val.y,
                     pair = paste0(lk$Sample.x, "/", lk$Sample.y), exp_fc = 0)
  rat <- dplyr::bind_rows(c(rat, list(like)))

  ## human alignment per (run x column x pair), then humans are DROPPED
  rat <- rat %>%
    dplyr::group_by(Run, time_channel, pair) %>%
    dplyr::mutate(ratio = ratio - median(ratio[species == "Human"], na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(species == "Yeast", !is.na(exp_fc))

  slope_tab <- rat %>%
    dplyr::group_by(Run, plex, time_channel, exp_fc) %>%
    dplyr::summarise(med = median(ratio, na.rm = TRUE),
                     q25 = quantile(ratio, 0.25, na.rm = TRUE),
                     q75 = quantile(ratio, 0.75, na.rm = TRUE),
                     n = dplyr::n(), .groups = "drop")
  slope_tab$Run <- factor(slope_tab$Run, levels = run_iv$Run)
  ## one facet column per run, on the acquisition-time axis of c-e: each run's
  ## slope plot spans its bar in e (bar_w hours around the run centre), and
  ## neighbouring columns meet halfway between runs. With no spacing inside a
  ## block, the c-e gap between blocks and the same block ends, f maps hours to
  ## inches exactly like c-e, so every plot sits under its run.
  labs_runs <- run_iv$run_lab
  slope_tab$run_lab <- factor(labs_runs[as.integer(slope_tab$Run)], levels = labs_runs)
  new_blk <- diff(run_iv$blk) != 0
  tcs <- sort(unique(slope_tab$time_channel))
  slope_tab$tc_lab <- factor(tp_labs[as.character(slope_tab$time_channel)], levels = tp_labs[as.character(tcs)])

  ticks_d <- sort(unique(slope_tab$exp_fc))
  rng_d <- range(c(slope_tab$q25, slope_tab$q75), na.rm = TRUE)
  to_hr <- function(fc, centre) centre + (fc - mean(rng_d)) * bar_w / diff(rng_d)   # fold change -> hours
  slope_tab$x <- to_hr(slope_tab$exp_fc, run_iv$center[as.integer(slope_tab$Run)])

  ## column edges: halfway between runs; block ends as in c-e (1% expansion)
  edges <- do.call(rbind, lapply(split(run_iv, run_iv$blk), function(b){
    lim <- range(blk_lim$x[blk_lim$blk == b$blk[1]])
    cut <- c(lim[1] - 0.01 * diff(lim), (head(b$center, -1) + tail(b$center, -1)) / 2, lim[2] + 0.01 * diff(lim))
    data.frame(run_lab = rep(b$run_lab, 2), x = c(head(cut, -1), tail(cut, -1)))
  }))
  edges$run_lab <- factor(edges$run_lab, levels = labs_runs)
  per_run <- data.frame(run_lab = factor(labs_runs, levels = labs_runs), centre = run_iv$center)
  ab_df <- merge(per_run, data.frame(tc_lab = factor(levels(slope_tab$tc_lab), levels = levels(slope_tab$tc_lab))))
  axis_df <- data.frame(per_run, tc_lab = factor(levels(slope_tab$tc_lab)[length(tcs)], levels = levels(slope_tab$tc_lab)))
  ## "..." in the gap after block 1, in the middle row, placed as in c-e
  dots_d <- data.frame(run_lab = factor(labs_runs[which(new_blk)], levels = labs_runs),
                       tc_lab = factor(levels(slope_tab$tc_lab)[ceiling(length(tcs) / 2)],
                                       levels = levels(slope_tab$tc_lab)))
  pd <- ggplot(slope_tab, aes(x = x)) +
    geom_blank(data = edges, aes(x = x), inherit.aes = FALSE) +
    geom_segment(data = ab_df, inherit.aes = FALSE,
                 aes(x = to_hr(rng_d[1], centre), xend = to_hr(rng_d[2], centre),
                     y = rng_d[1], yend = rng_d[2]),
                 linetype = "dashed", color = "grey55", linewidth = 0.25) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = factor(time_channel)),
                alpha = 0.45, colour = "grey30", linewidth = 0.2) +
    geom_line(aes(y = med), colour = "grey40", linewidth = 0.4) +
    geom_point(aes(y = med), shape = 21, colour = "black", fill = "white",
               size = 1.5, stroke = 0.35) +
    # x axis line under each run's plot only
    geom_segment(data = axis_df, inherit.aes = FALSE,
                 aes(x = to_hr(rng_d[1], centre), xend = to_hr(rng_d[2], centre), y = -Inf, yend = -Inf),
                 linewidth = 0.25) +
    geom_text(data = dots_d, aes(x = Inf, y = mean(rng_d)), label = "...",
              inherit.aes = FALSE, hjust = -1.1, fontface = "bold", size = paper_geom_text(11)) +
    facet_grid(tc_lab ~ run_lab, switch = "y", scales = "free_x", space = "free_x") +
    scale_fill_manual(values = greens, guide = "none") +
    # ticks 0-3 under each run, labelled in fold change
    scale_x_continuous(breaks = function(lim) {
                         centres <- run_iv$center[run_iv$center > lim[1] & run_iv$center < lim[2]]
                         as.vector(outer(ticks_d, centres, to_hr))
                       },
                       labels = function(b) {
                         centres <- run_iv$center[sapply(b, function(v) which.min(abs(run_iv$center - v)))]
                         format(round(mean(rng_d) + (b - centres) * diff(rng_d) / bar_w))
                       }) +
    scale_y_continuous(breaks = ticks_d) +
    coord_cartesian(ylim = rng_d + c(-1, 1) * 0.05 * diff(rng_d), expand = FALSE, clip = "off") +
    labs(x = expression(paste("Expected ", Log[2], " fold change")),     # spaces in strings, not ~:
         y = expression(paste(Log[2], " precursor fold change"))) +       # fewer pieces in Illustrator
    paper_theme() +
    theme(strip.text.x = element_blank(), strip.background = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = paper_size("text")),
          axis.text = element_text(size = 6), axis.line.x = element_blank(),
          panel.spacing.x = unit(ifelse(new_blk, gap_in, 0), "in"), panel.spacing.y = unit(3, "pt"),
          plot.margin = margin(2, 3, 3, 3))

  ## ---- legends above c: plexDIA x timePlex color grid, timePlex greens -------
  ## placed and sized as in the manuscript: the grid from ~15% to ~61% of the
  ## width (small tiles, 8-pt labels, 9-pt titles), the greens from ~67%
  ## (offsets are written above the tiles; no x guide)
  plex_chans <- sort(unique(dat$plex_channel))
  key <- expand.grid(time_channel = tp_chans, plex_channel = plex_chans)
  key$fill_key <- paste0("d", key$plex_channel, " tp", key$time_channel)
  leg_fill <- ggplot(key, aes(x = factor(plex_channel), y = factor(time_channel, levels = rev(tp_chans)),
                              fill = fill_key)) +
    geom_tile(width = 0.4, height = 0.6) +
    annotate("text", x = seq_along(plex_chans), y = Inf, label = paste0("\u0394", plex_chans),
             vjust = -0.35, size = paper_geom_text(8)) +
    scale_fill_manual(values = fills$vals, guide = "none") +
    scale_y_discrete(labels = rev(tp_labs)) +
    labs(x = NULL, y = "timePlex\nchannel", title = "plexDIA mass offsets (Da)") +
    guides(x = "none") +
    coord_cartesian(clip = "off") +
    paper_theme() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, lineheight = 0.9),
          plot.title = element_text(size = paper_size("title"), hjust = 0.5, margin = margin(0, 0, 12, 0)),
          plot.margin = margin(0, 0, 0, 0))
  leg_h <- 0.7
  ## green key in inch coordinates: title right-aligned left of the dots, rows 0.13 in apart
  gk <- data.frame(tc = factor(tp_chans, levels = tp_chans),
                   y = 0.33 + 0.13 * (rev(seq_along(tp_chans)) - (length(tp_chans) + 1) / 2))
  leg_green <- ggplot(gk, aes(x = 0.62, y = y)) +
    geom_point(aes(colour = tc), size = 2.2) +
    geom_text(aes(x = 0.71, label = tp_labs[as.character(tc)]), hjust = 0, size = paper_geom_text(8)) +
    annotate("text", x = 0.53, y = 0.33, label = "timePlex\nchannel", hjust = 1,
             lineheight = 0.9, size = paper_geom_text(paper_size("title"))) +
    scale_colour_manual(values = greens, guide = "none") +
    coord_cartesian(xlim = c(0, 1.3), ylim = c(0, leg_h), expand = FALSE, clip = "off") +
    theme_void()
  ## both keys as free elements: aligned in one patchwork row, the tile panel collapses
  legends <- patchwork::wrap_plots(patchwork::plot_spacer(), patchwork::wrap_elements(full = leg_fill),
                                   patchwork::plot_spacer(), patchwork::wrap_elements(full = leg_green),
                                   patchwork::plot_spacer(), nrow = 1,
                                   widths = grid::unit(c(1.05, 3.3, 0.3, 1.3, 1.14), "in"))

  ## ---- panels and the aligned c-f block -------------------------------------
  fig_w <- 7.09
  heights <- c(c = 1.75, d = 1.3, e = 1.25, f = 1.9)
  ## the legend row is a fixed 0.7 in; the panels share the rest
  pc_c <- patchwork::wrap_plots(patchwork::wrap_elements(full = legends), pa, ncol = 1,
                                heights = grid::unit(c(leg_h, 1.75), c("in", "null")))
  out <- list(
    fig6c = paper_panel(pc_c, file.path(out_dir, "Fig6c.pdf"), width = fig_w, height = leg_h + heights[["c"]]),
    fig6d = paper_panel(pb, file.path(out_dir, "Fig6d.pdf"), width = fig_w, height = heights[["d"]]),
    fig6e = paper_panel(pc, file.path(out_dir, "Fig6e.pdf"), width = fig_w, height = heights[["e"]]),
    fig6f = paper_panel(pd, file.path(out_dir, "Fig6f.pdf"), width = fig_w, height = heights[["f"]]))
  ## f's row labels sit left of its y axis; kept out of the column alignment,
  ## they no longer push the y axes of c-e away from their panels
  cf <- patchwork::wrap_plots(patchwork::wrap_elements(full = legends), pa, pb, pc,
                              patchwork::free(pd, type = "space", side = "l"), ncol = 1,
                              heights = grid::unit(c(leg_h, heights), c("in", rep("null", length(heights)))))
  out$fig6cf <- paper_panel(cf, file.path(out_dir, "Fig6c-f.pdf"), width = fig_w,
                            height = leg_h + sum(heights))

  c(out, list(counts = counts, rt_pts = rt_pts, slope_tab = slope_tab))
}
