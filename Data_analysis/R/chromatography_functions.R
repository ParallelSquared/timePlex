################################################################################
# Chromatography: RT accuracy, carryover between timePlex channels, peak widths
################################################################################


## ---- Library vs fine-tuned RTs (Fig. 2f) ------------------------------------

pD_plot_RTs_TP <- function(dat=benchmarking_dat8, run="JD0311", filter_level=0.01, use_best_channel=T, lib_fpath,
                           tp_labels = c("\u03940 min", "\u03945 min", "\u039411 min"),
                           out_file = "Figs/Fig2/Fig2f.pdf"){
  paper_figure(out_file)

  dat_specific <- dat[dat$Run==run,]

  # original library RT (Tr_recalibrated), matched on stripped sequence + charge
  lib <- fread(lib_fpath, select = c("PeptideSequence", "PrecursorCharge", "transition_group_id", "ProteinName", "Tr_recalibrated"))
  lib$seqcharge <- paste0(lib$PeptideSequence, "_", lib$PrecursorCharge)
  lib <- lib %>% distinct(transition_group_id, ProteinName, .keep_all=T) %>% dplyr::select("seqcharge", "Tr_recalibrated")
  dat_specific <- dat_specific %>% left_join(lib, by = "seqcharge")

  if(use_best_channel==T){
    dat_specific <- dat_specific %>% dplyr::filter(Best_pDChannel_Qvalue<0.01 & is_decoy==F) %>% dplyr::filter(Qvalue < filter_level)
  } else{
    dat_specific <- dat_specific %>% dplyr::filter(Qvalue<0.01 & is_decoy==F)
  }

  chan0 <- dat_specific[dat_specific$time_channel==1,] %>% dplyr::select("seqcharge","lib_rt") %>%
    dplyr::rename("lib_rt_col1" = "lib_rt")
  dat_specific_unified <- dat_specific %>% inner_join(chan0, by =c("seqcharge"="seqcharge"))

  dat_specific_unified <- dat_specific_unified[!is.na(dat_specific_unified$Tr_recalibrated),]

  dat_specific_unified <- dat_specific_unified %>%
    dplyr::group_by(time_channel) %>%
    dplyr::mutate(lowess_fit = lowess(Tr_recalibrated, rt, f = 0.1)$y)

  dat_specific_unified$d_rt_lowess <- dat_specific_unified$rt-dat_specific_unified$lowess_fit

  ## ---- Fig. 2f: empirical RT vs the original library iRT (left) and vs the
  ## fine-tuned library RT of the middle channel (right), per timePlex channel
  tp_cols <- c("#16330E", "#4D9934", "#E1F7DA")
  rt_panel <- function(x, xlab, legend){
    ggplot(dat_specific_unified, aes(y = rt, x = .data[[x]])) +
      paper_raster(geom_point(alpha = 0.1, shape = 21, size = 0.84, stroke = 0.14, aes(fill = as.factor(time_channel)))) +
      scale_fill_manual(values = tp_cols, labels = tp_labels, name = "timePlex") +
      labs(x = xlab, y = "Empirical RT, min") +
      guides(fill = guide_legend(override.aes = list(alpha = 0.6, size = 1.3))) +
      paper_theme() +
      theme(legend.position = if (legend) "inside" else "none",
            legend.position.inside = c(0.02, 1), legend.justification = c(0, 1),
            legend.key.spacing.y = unit(0, "pt"))
  }
  arrow_panel <- ggplot() +
    annotate("segment", x = 0.08, xend = 0.92, y = 0, yend = 0, linewidth = 0.4,
             arrow = grid::arrow(length = unit(3, "pt"), type = "closed")) +
    annotate("text", x = 0.5, y = 0.06, label = "Transfer\nlearning", vjust = 0, lineheight = 0.9,
             size = paper_geom_text(7)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1), expand = FALSE) +
    theme_void()

  p <- patchwork::wrap_plots(rt_panel("Tr_recalibrated", "Original library iRT (min)", legend = TRUE) +
                               coord_cartesian(xlim = c(-50, 150)),
                             arrow_panel,
                             rt_panel("lib_rt_col1", "Finetuned library RT (min)", legend = FALSE),
                             nrow = 1, widths = c(1, 0.32, 1))
  fig2f <- paper_panel(p, out_file, width = 4.54, height = 2.05)

  dat_specific_unified$d_RT <- dat_specific_unified$rt-dat_specific_unified$lib_rt_col1
  chan1 <- dat_specific_unified[dat_specific_unified$time_channel==1,] %>% dplyr::select("seqcharge","rt") %>% dplyr::rename("rt_first"="rt")

  dat_specific_unified2 <- dat_specific_unified %>% inner_join(chan1, by =c("seqcharge"="seqcharge"))
  dat_specific_unified2$delta_RT_first <- dat_specific_unified2$rt - dat_specific_unified2$rt_first

  invisible(list(fig2f = fig2f, data = dat_specific_unified2))
}

## ---- Carryover into unintended timePlex channels (Extended Data Fig. 1) -----

compare_carryover <- function(dat_all_carryover, meta, use_best_channel = F, filter_level = 0.01,
                              quant_col = "plex_Area", min_intensity = 1, min_reps = 2,
                              out_dir = "Figs/ExtDataFig1") {
  paper_figure(out_dir)

  ## Metadata: one row per run; the condition comes from Notes (smart quotes stripped),
  ## ordered by DDM % so "No DDM" is panel a
  if (is.character(meta)) meta <- data.frame(data.table::fread(meta))
  meta <- meta %>%
    dplyr::filter(!is.na(Run), trimws(Run) != "") %>%   # spreadsheet exports can end in blank rows
    dplyr::mutate(Run       = trimws(Run),
                  Condition = trimws(gsub("[\"“”]", "", V9)),
                  ddm_pct   = ifelse(grepl("^no", Condition, ignore.case = T), 0,
                                     as.numeric(stringr::str_extract(Condition, "[0-9.]+")))) %>%
    dplyr::arrange(ddm_pct, Rep)
  cond_levels <- unique(meta$Condition)
  meta$Condition <- factor(meta$Condition, levels = cond_levels)

  if (!quant_col %in% names(dat_all_carryover)) {
    stop(sprintf("column '%s' not in the data; add it to the fread(select = ...) in Jmod_read_folder", quant_col))
  }

  dat <- dat_all_carryover
  dat <- dat %>% dplyr::mutate("seqcharge" = untag_prec)
  dat <- dat %>% dplyr::mutate("prot" = protein)
  dat <- dat %>% dplyr::group_by(Run, time_channel, seqcharge) %>%
    dplyr::mutate("Best_pDChannel_Qvalue" = min(Qvalue, na.rm = T)) %>% ungroup()
  if (use_best_channel == T) {
    dat <- dat %>% dplyr::filter(Best_pDChannel_Qvalue < 0.01) %>% dplyr::filter(Qvalue < filter_level)
  } else {
    dat <- dat %>% dplyr::filter(Qvalue < 0.01)
  }

  dat <- dat %>% dplyr::filter(Run %in% meta$Run)
  missing_runs <- setdiff(meta$Run, unique(dat$Run))
  if (length(missing_runs) > 0) {
    warning(sprintf("runs in meta but not in data: %s", paste(missing_runs, collapse = ", ")))
  }
  dat$quant <- dat[[quant_col]]
  dat$time_channel <- as.integer(round(dat$time_channel))

  ## One row per run x precursor: intensity and RT in each time channel
  ## (if a precursor has several rows in a channel, keep its most intense one)
  dat_d <- dat %>%
    dplyr::arrange(dplyr::desc(quant)) %>%
    dplyr::distinct(Run, untag_prec, time_channel, .keep_all = T) %>%
    dplyr::select(Run, untag_prec, time_channel, quant, rt) %>%
    tidyr::pivot_wider(names_from = time_channel, values_from = c(quant, rt),
                       names_glue = "{.value}_ch{time_channel}")
  for (col in c("quant_ch0", "quant_ch1", "rt_ch0", "rt_ch1")) {
    if (!col %in% names(dat_d)) dat_d[[col]] <- NA_real_
  }

  ## time_channel 0 = timePlex ch 1, time_channel 1 = timePlex ch 2.
  ## A precursor counts as quantifiable in a channel only if its intensity is > min_intensity.
  dat_d <- dat_d %>%
    dplyr::mutate(
      ch0 = ifelse(!is.na(quant_ch0) & quant_ch0 > min_intensity, quant_ch0, NA_real_),
      ch1 = ifelse(!is.na(quant_ch1) & quant_ch1 > min_intensity, quant_ch1, NA_real_),
      rat = log2(ch1 / ch0)) %>%
    dplyr::inner_join(meta %>% dplyr::select(Run, Condition, Rep), by = "Run")

  empty_conds <- setdiff(cond_levels, as.character(unique(dat_d$Condition)))
  if (length(empty_conds) > 0) {
    warning(sprintf("no data for condition(s): %s", paste(empty_conds, collapse = ", ")))
  }
  cond_levels <- setdiff(cond_levels, empty_conds)

  ## Quantifiable coverage per run, then mean +/- SD over replicates
  set_levels <- c("timePlex ch 1", "timePlex ch 2", "Both")
  set_cols   <- c("timePlex ch 1" = "#00BFC4", "timePlex ch 2" = "#F87664", "Both" = "grey60")

  counts <- dat_d %>%
    dplyr::group_by(Run, Condition, Rep) %>%
    dplyr::summarise(`timePlex ch 1` = sum(!is.na(ch0)),
                     `timePlex ch 2` = sum(!is.na(ch1)),
                     Both            = sum(!is.na(ch0) & !is.na(ch1)), .groups = "drop") %>%
    tidyr::pivot_longer(dplyr::all_of(set_levels), names_to = "Set", values_to = "n") %>%
    dplyr::mutate(Set = factor(Set, levels = set_levels))

  count_sum <- counts %>%
    dplyr::group_by(Condition, Set) %>%
    dplyr::summarise(mean_n = mean(n), sd_n = sd(n), n_rep = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(sd_n = dplyr::coalesce(sd_n, 0))

  ## Ratios: only precursors quantified in both channels in >= min_reps reps of a condition,
  ## one point per precursor = its median log2 ratio (and median ch 1 RT) over those reps
  dat_both <- dat_d %>% dplyr::filter(!is.na(rat))
  dat_prec <- dat_both %>%
    dplyr::group_by(Condition, untag_prec) %>%
    dplyr::summarise(rat = median(rat), rt_ch0 = median(rt_ch0), n_reps = dplyr::n(),
                     .groups = "drop") %>%
    dplyr::filter(n_reps >= min_reps)

  rat_sum <- dat_prec %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(n = dplyr::n(), median_log2 = median(rat), .groups = "drop")

  ## Shared axes across conditions so the panels are directly comparable
  bar_max     <- max(count_sum$mean_n + count_sum$sd_n, counts$n) * 1.12
  rat_lim     <- range(dat_prec$rat, na.rm = T)
  rt_lim      <- range(dat_prec$rt_ch0, na.rm = T)
  hist_breaks <- seq(rat_lim[1], rat_lim[2], length.out = 61)

  ## One row per condition (a = No DDM, b = 0.015% DDM, ...): coverage | ratio vs RT | histogram
  rows <- lapply(seq_along(cond_levels), function(i) {
    cond <- cond_levels[i]
    cs   <- count_sum %>% dplyr::filter(Condition == cond)
    cr   <- counts %>% dplyr::filter(Condition == cond)
    db   <- dat_prec %>% dplyr::filter(Condition == cond)
    med  <- median(db$rat)
    med_lab <- sprintf('paste("median ", Log[2], " ratio = %s  (ch 2 = %.2f%% of ch 1)")',
                       sprintf("%.2f", med), 100 * 2^med)

    p_bar <- ggplot(cs, aes(x = Set, y = mean_n, fill = Set)) +
      geom_col(width = 0.7, color = "black", linewidth = 0.2) +
      geom_errorbar(aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n), width = 0.25, linewidth = 0.25) +
      geom_point(data = cr, aes(x = Set, y = n), inherit.aes = F, size = 0.6,
                 position = position_jitter(width = 0.12, height = 0, seed = 1)) +
      geom_text(aes(y = mean_n + sd_n, label = scales::comma(round(mean_n))),
                vjust = -0.5, size = paper_geom_text(7)) +
      scale_fill_manual(values = set_cols, guide = "none") +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0))) +
      coord_cartesian(ylim = c(0, bar_max)) +
      labs(title = cond, x = NULL, y = "Quantified precursors") +
      paper_theme() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            plot.title = element_text(hjust = 0),
            plot.margin = margin(2, 2, 2, 10))     # room for the panel letter left of the y title

    p_dots <- ggplot(db, aes(x = rt_ch0, y = rat)) +
      paper_raster(geom_point(alpha = 0.25, size = 0.6, stroke = 0)) +
      geom_hline(yintercept = med, linetype = "dashed", color = "red", linewidth = 0.4) +
      annotate("text", x = rt_lim[1], y = rat_lim[2], hjust = 0, vjust = 1, size = paper_geom_text(),
               label = med_lab, parse = TRUE) +
      coord_cartesian(xlim = rt_lim, ylim = rat_lim) +
      labs(x = "RT in timePlex ch 1 (min)", y = bquote(Log[2]*", ch 2 / ch 1")) +
      paper_theme()

    p_hist <- ggplot(db, aes(y = rat)) +
      geom_histogram(breaks = hist_breaks, alpha = 0.6, color = "black", linewidth = 0.2, fill = "grey") +
      geom_hline(yintercept = med, linetype = "dashed", color = "red", linewidth = 0.4) +
      coord_cartesian(ylim = rat_lim) +
      scale_x_continuous(n.breaks = 3, labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
      paper_theme() +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      labs(x = "Count")

    # x titles right under their own axis, not below the bar chart's angled labels
    p <- patchwork::wrap_plots(p_bar, patchwork::free(p_dots, type = "label"),
                               patchwork::free(p_hist, type = "label"), ncol = 3, widths = c(1.4, 3.5, 1)) +
      # small margins above and below, so the rows sit a little closer together
      patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 5.5, 2, 5.5)))
    paper_panel(p, file.path(out_dir, paste0("ExtDataFig1", letters[i], ".pdf")), width = 7.09, height = 1.78)
  })
  names(rows) <- paste0("ed1", letters[seq_along(rows)])

  invisible(c(rows, list(counts = counts, count_summary = count_sum, ratio_summary = rat_sum,
                         dat_prec = dat_prec, dat_d = dat_d)))
}


## ---- Peak widths per column and implementation (Extended Data Fig. 5a) ------

Jmod_FWHMs <- function(uPAC=benchmarking_dat_uPAC$Amount_8, IO=benchmarking_dat_IO$Amount_8,
                       system_labels = c(IO = "IonOpticks Aurora", uPAC = "Thermo \u00b5PAC Neo"),
                       out_file = "Figs/ExtDataFig5/ExtDataFig5a.pdf"){
  paper_figure(out_file)

  uPAC$system <- "uPAC"
  IO$system <- "IO"
  joined <- rbind.fill(uPAC, IO)
  joined <- joined %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    )) %>% dplyr::filter(Type =="LF\n3-timePlex")
  joined <- joined %>% dplyr::filter(Qvalue<0.01)
  joined$FWHM <- joined$MS1_FWHM*60
  joined$Column <- joined$Column

  med_joined <- joined %>% dplyr::group_by(seqcharge,system,Column) %>%
    dplyr::summarize(med_FWHM = median(FWHM,na.rm=T),
                     med_quant = median(MS1_Area,na.rm=T)) %>% dplyr::ungroup() %>% na.omit() %>%
    dplyr::add_count(seqcharge) %>% dplyr::filter(n==6)

  y_d0d4 <- med_joined$med_FWHM %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(0.0,.995)))

  med_joined$system <- factor(system_labels[med_joined$system], levels = system_labels)

  p <- ggplot(med_joined, aes(x = med_FWHM, y = as.factor(Column))) +
    facet_grid(rows = vars(system)) +
    geom_density_ridges2(
      alpha = 0.75,
      scale = 1.0,
      color = "black",
      fill="grey70",
      linewidth = 0.25
    ) +
    stat_summary(
      aes(y = as.factor(Column)),
      fun = median,
      geom = "crossbar",
      width = 0.33,
      color = "black",
      linewidth = 0.25,
      middle.linewidth = 0.5
    ) + labs(x="FWHM (seconds)", y = "LC column")+
    coord_cartesian(xlim = y_d0d4_val) +
    paper_theme_box() +
    theme(legend.position = "none")

  invisible(list(ed5a = paper_panel(p, out_file, width = 2.364, height = 2.85), fwhm = med_joined))
}
