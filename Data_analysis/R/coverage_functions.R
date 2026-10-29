################################################################################
# Identification coverage: precursors and proteins per run, sample and channel
################################################################################


## ---- Missing-species check (Fig. 3a) ----------------------------------------

Missing_species_hist_quant <- function(dat, extra="", use_best_channel=F, filter_level=0.01,
                                       channel_labels = c("\u03940 min, Human only",
                                                          "\u03944 min, Human & Yeast",
                                                          "\u03948 min, Human & Yeast"),
                                       out_file = "Figs/Fig3/Fig3a.pdf"){
  paper_figure(out_file)
  dat$seqcharge <- dat$untag_prec
  dat$prot <- dat$protein

  dat <- dat %>% dplyr::group_by(Run, time_channel, seqcharge) %>%
    dplyr::mutate("Best_pDChannel_Qvalue" = min(Qvalue,na.rm=T)) %>% ungroup()

  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(Best_pDChannel_Qvalue<0.01 & Qvalue<filter_level)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<0.01)
  }


  dat_lim <- dat %>% dplyr::add_count(untag_prec, channel)
  dat_lim <- dat_lim %>% dplyr::filter(n==1 & Proteotypic==T) #removes homologous peptides
  dat_lim <- dat_lim %>% dplyr::mutate(species=ifelse(grepl("HUMAN",ProteinName),"Human",
                                                      ifelse(grepl("YEAST",ProteinName),"Yeast","Arabidopsis")))

  # one facet per timePlex channel (sorted), labelled by offset and sample composition
  dat_lim <- dat_lim %>%
    dplyr::filter(MS1_Area_rtnorm > 1) %>%
    dplyr::mutate(species     = factor(species, levels = c("Human", "Yeast", "Arabidopsis")),
                  channel_lab = factor(channel_labels[match(channel, sort(unique(channel)))],
                                       levels = channel_labels))

  dat_lim_counts <- dat_lim %>%
    dplyr::group_by(species, channel_lab) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    dplyr::mutate(label = sprintf('paste("%s: ", italic(n), " = %s")', species, scales::comma(count)),
                  line  = as.integer(species))

  p <- ggplot(dat_lim, aes(x = log10(MS1_Area_rtnorm), fill = species)) +
    geom_histogram(alpha = 0.75, position = "identity", colour = "black", linewidth = 0.08) +
    facet_grid(~channel_lab) +
    scale_fill_manual(values = paper_cols$species, name = NULL) +
    scale_x_continuous(limits = c(0, NA), breaks = seq(0, 8, 2)) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = bquote(Log[10]*", MS1 precursor intensity"), y = "#Precursors") +
    guides(fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.12))) +
    paper_theme_box(base_size = 7.5, frame = 0.12) +
    theme(axis.text = element_text(size = 7), legend.text = element_text(size = 7),
          legend.key.size = unit(9, "pt"), legend.key.spacing.y = unit(2, "pt"))

  # species counts as evenly spaced lines in the top-left of each facet
  yr <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  dat_lim_counts$y <- yr[2] - (0.06 + 0.14 * (dat_lim_counts$line - 1)) * diff(yr)
  p <- p + geom_text(data = dat_lim_counts, aes(x = 0.3, y = y, label = label), hjust = 0, vjust = 1,
                     parse = TRUE, size = paper_geom_text(6.5), inherit.aes = FALSE)

  fig3a <- paper_panel(p, out_file, width = 7.09, height = 1.45)
  invisible(list(fig3a = fig3a, counts = dat_lim_counts))
}


## ---- 3-plex benchmarks (Fig. 3d-e, Extended Data Fig. 3) --------------------

tp_coverage <- function(dat, extra = "", use_best_channel = FALSE, filter_level = 0.01,
                        sample_id = c("digit", "full"),
                        out_files = c(prec = "Figs/Fig3/Fig3d.pdf", prot = "Figs/Fig3/Fig3e.pdf"),
                        width = 2.92, height = 2.4) {
  paper_figure(out_files)

  dat <- dat %>% dplyr::mutate("Sample" = ifelse(Sample %in% c("A","B","C"), Sample,
                                                                   ifelse(grepl("C1|D1",Sample),"A",
                                                                          ifelse(grepl("C2|D2",Sample),"B","C"))))

  sample_id <- match.arg(sample_id, c("digit", "full"))

  narrow_types <- c("3-plexDIA\n3-timePlex", "LF\n3-timePlex", "3-plexDIA")
  x_limits     <- c("LF\nNo-plex", "3-plexDIA", "LF\n3-timePlex", "3-plexDIA\n3-timePlex")
  base_colors  <- unname(paper_cols$samples3)

  ## ------------------------------------------------------------------
  ## helpers
  ## ------------------------------------------------------------------

  # Sample IDs look like C1/C2/C3/D1/D2/D3: the digit is the multiplexed
  # sample, the letter is the column/prep replicate. Collapsing to the digit
  # gives the 3 fill levels the figure is designed around and turns C/D into
  # replicates, which is what makes the SE meaningful.
  set_sample_levels <- function(df) {
    v    <- as.character(df$Sample)
    digs <- sub("^[^0-9]*", "", v)
    use_digit <- sample_id == "digit" &&
      !any(digs == "" | is.na(digs)) &&
      length(unique(digs)) <= 3
    if (sample_id == "digit" && !use_digit) {
      warning("tp_coverage(): could not reduce Sample to <=3 numeric IDs (",
              paste(sort(unique(v)), collapse = ", "),
              "); using full sample IDs instead.", call. = FALSE)
    }
    if (use_digit) v <- digs

    lv <- sort(unique(v))
    lv <- c(intersect(c("A", "B", "C", "1", "2", "3"), lv),
            setdiff(lv, c("A", "B", "C", "1", "2", "3")))
    df$Sample <- factor(v, levels = lv)
    df
  }

  make_palette <- function(f) {
    lv   <- levels(factor(f))
    lv
    cols <- if (length(lv) <= length(base_colors)) base_colors[seq_along(lv)]
    else colorRampPalette(base_colors)(length(lv))
    stats::setNames(cols, lv)
  }

  summarise_counts <- function(dat_counts) {
    dat_counts %>%
      dplyr::group_by(Sample, Amount, Type, Label) %>%
      dplyr::summarise(
        mean_count = mean(n, na.rm = TRUE),
        n_rep      = dplyr::n(),
        se         = stats::sd(n, na.rm = TRUE) / sqrt(dplyr::n()),
        .groups    = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(Sample), Type) %>%
      dplyr::group_by(Type) %>%
      dplyr::mutate(stacked_mean = cumsum(mean_count)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        stacked_mean = ifelse(Type == "LF\nNo-plex", mean_count, stacked_mean),
        er_max       = stacked_mean + se,
        er_min       = stacked_mean - se
      )
  }

  # One label per bar: the stack total for stacked bars, each bar for LF no-plex
  make_text <- function(dat_summary) {
    dat_summary %>%
      dplyr::group_by(Type, Sample, Label) %>%
      dplyr::summarise(
        text_value = round(max(stacked_mean)),
        text_pos   = max(c(er_max, stacked_mean), na.rm = TRUE),
        .groups    = "drop"
      ) %>%
      dplyr::group_by(Type) %>%
      dplyr::filter(Type == "LF\nNo-plex" | text_pos == max(text_pos)) %>%
      dplyr::ungroup()
  }

  build_plot <- function(dat_summary, df_text, ylab, my_colors) {
    err <- dat_summary %>% dplyr::filter(!is.na(er_min), !is.na(er_max))

    ggplot(dat_summary, aes(x = Type, y = mean_count, fill = Sample)) +
      scale_x_discrete(limits = intersect(x_limits, dat_summary$Type[dat_summary$mean_count > 0])) +

      geom_bar(data = dat_summary %>% dplyr::filter(Type %in% narrow_types),
               stat = "identity", color = "black", width = 0.35, linewidth = 0.25) +

      geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
               stat = "identity", position = position_dodge(width = 1.17),
               color = "black", width = 1, linewidth = 0.25) +

      geom_errorbar(data = err %>% dplyr::filter(Type %in% narrow_types),
                    aes(ymin = er_min, ymax = er_max),
                    width = 0.18, linewidth = 0.25) +

      geom_errorbar(data = err %>% dplyr::filter(Type == "LF\nNo-plex"),
                    aes(ymin = er_min, ymax = er_max),
                    width = 0.5, linewidth = 0.25,
                    position = position_dodge(1.17)) +

      geom_text(data = df_text %>% dplyr::filter(Type %in% narrow_types),
                aes(x = Type, y = text_pos, label = comma(text_value)),
                inherit.aes = FALSE, vjust = -0.3, size = paper_geom_text(paper_size("text") - 1)) +

      geom_text(data = df_text %>% dplyr::filter(Type == "LF\nNo-plex"),
                aes(x = Type, y = text_pos, label = comma(text_value), group = Sample),
                vjust = -0.3, size = paper_geom_text(paper_size("text") - 1), position = position_dodge(1.17)) +

      scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.06))) +
      scale_fill_manual(values = my_colors) +
      labs(y = ylab, x = NULL, fill = "Samples") +
      paper_theme() +
      # sizes as in the manuscript (Fig. 3d-e and Extended Data Fig. 3): larger
      # bar labels on x, sample letters as large as the legend title
      theme(axis.text.y            = element_text(size = 7),
            axis.text.x            = element_text(size = paper_size("title") + 1),
            axis.title.y           = element_text(size = 9),
            legend.text            = element_text(size = paper_size("title")),
            legend.position        = "inside",
            legend.position.inside = c(0.04, 0.97),
            legend.justification   = c(0, 1),
            legend.direction       = "horizontal",
            legend.key.spacing.x   = unit(2.5, "pt")) +
      guides(fill = guide_legend(title = "Samples",
                                 label.position = "bottom",
                                 title.position = "top", title.hjust = 0.5))
  }

  ## ------------------------------------------------------------------
  ## filtering / annotation
  ## ------------------------------------------------------------------

  if (isTRUE(use_best_channel)) {
    # NOTE: BestChannel cut kept at the original hard-coded 0.01, not filter_level
    dat <- dat %>%
      dplyr::filter(Best_pDChannel_Qvalue < 0.01) %>%
      dplyr::filter(Qvalue < filter_level)
  } else {
    dat <- dat %>% dplyr::filter(Qvalue < 0.01)
  }

  dat <- dat %>%
    dplyr::mutate(Type = dplyr::case_when(
      timePlex == TRUE  & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE  & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE                              ~ "LF\nNo-plex"
    )) %>%
    dplyr::filter(!Repeatability)          # drops NA, unlike dat[dat$Repeatability == F, ]

  if ("Implementation" %in% colnames(dat)) {
    dat <- dat %>%
      dplyr::mutate(Type = paste0(Type, "_", Implementation)) %>%
      dplyr::filter(Type == "3-plexDIA" | Type == "LF\nNo-plex")
  }


  ## ------------------------------------------------------------------
  ## precursor panel
  ## ------------------------------------------------------------------

  dat_counts_prec <- dat %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    dplyr::ungroup()

  dat_counts_prec <- dat_counts_prec %>% dplyr::mutate(Label = ifelse(Type == "3-plexDIA", "mTRAQ", Label))


  if (!"3-plexDIA" %in% dat$Type) {
    missing_3plexDIA <- dat %>%
      dplyr::distinct(Sample, Column, Amount, Label) %>%
      dplyr::mutate(Type = "3-plexDIA", n = 0, Run = "Placeholder_Run")
    dat_counts_prec <- dplyr::bind_rows(dat_counts_prec, missing_3plexDIA)
  }

  my_colors       <- make_palette(dat_counts_prec$Sample)

  summary_prec <- summarise_counts(dat_counts_prec)
  text_prec    <- make_text(summary_prec)

  p_prec <- build_plot(summary_prec, text_prec, "Precursor data points / Run", my_colors)

  ## ------------------------------------------------------------------
  ## protein panel
  ## ------------------------------------------------------------------

  dat_counts_prot <- dat %>%
    dplyr::filter(Protein_Qvalue < 0.01) %>%
    dplyr::distinct(Sample, Run, Column, Amount, Label, Type, prot, .keep_all = TRUE) %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    dplyr::ungroup()

  dat_counts_prot <- dat_counts_prot %>% dplyr::mutate(Label = ifelse(Type == "3-plexDIA", "mTRAQ", Label))


  dat_counts_prot <- set_sample_levels(dat_counts_prot)
  my_colors_prot  <- make_palette(dat_counts_prot$Sample)

  summary_prot <- summarise_counts(dat_counts_prot)
  text_prot    <- make_text(summary_prot)

  p_prot <- build_plot(summary_prot, text_prot, "Protein data points / Run", my_colors_prot)

  invisible(list(prec   = paper_panel(p_prec, out_files[["prec"]], width, height),
                 prot   = paper_panel(p_prot, out_files[["prot"]], width, height),
                 counts = dat_counts_prot))
}


## ---- 3-column vs single-column spray (Extended Data Fig. 2) -----------------

tp_coverage_flow <- function(dat, extra = "", filter_level = 0.01,
                             comparisons = list("~200 nL/min per column" = c("600_3", "200_1"),
                                                "~66 nL/min per column"  = c("200_3", "66_1")),
                             out_dir = "Figs/ExtDataFig2") {
  paper_figure(out_dir)

  # Setup comes from the meta-sheet Notes column: "600 nL/min, 3 columns spraying" -> "600_3".
  # Each comparison is a pair; intensity ratios are first / second setup (3 cols / 1 col).

  dat <- dat %>% dplyr::mutate("Sample" = ifelse(Sample %in% c("A","B","C"), Sample,
                                                 ifelse(grepl("C1|D1",Sample),"A",
                                                        ifelse(grepl("C2|D2",Sample),"B","C"))))

  my_colors <- paper_cols$samples3
  bar_width <- 0.7
  dodge     <- position_dodge(width = 0.8, preserve = "single")

  ## ------------------------------------------------------------------
  ## helpers
  ## ------------------------------------------------------------------

  # "600_3" -> "600 nL/min\n3 columns"
  setup_label <- function(key) {
    flow  <- sub("_.*", "", key)
    ncols <- as.integer(sub(".*_", "", key))
    paste0(flow, " nL/min\n", ncols, ifelse(ncols == 1, " column", " columns"))
  }

  summarise_counts <- function(dat_counts) {
    dat_counts %>%
      dplyr::mutate(Sample = factor(Sample, levels = names(my_colors))) %>%
      dplyr::group_by(Comparison, Setup, Sample) %>%
      dplyr::summarise(
        mean_count = mean(n, na.rm = TRUE),
        n_rep      = dplyr::n(),
        se         = stats::sd(n, na.rm = TRUE) / sqrt(dplyr::n()),
        .groups    = "drop"
      ) %>%
      dplyr::mutate(er_max   = mean_count + se,
                    er_min   = mean_count - se,
                    text_pos = pmax(mean_count, er_max, na.rm = TRUE))
  }

  # Bars per setup, one facet per flow-rate comparison; the sample legend is
  # only shown when there is more than one sample
  build_plot <- function(dat_summary, ylab) {
    ggplot(dat_summary, aes(x = Setup, y = mean_count, fill = Sample)) +
      geom_col(position = dodge, color = "black", width = bar_width, linewidth = 0.25) +
      geom_errorbar(data = dat_summary %>% dplyr::filter(!is.na(se)),
                    aes(ymin = er_min, ymax = er_max),
                    width = bar_width / 2, linewidth = 0.25, position = dodge) +
      geom_text(aes(y = text_pos, label = comma(mean_count, accuracy = 1), group = Sample),
                vjust = -0.4, size = paper_geom_text(), position = dodge) +
      facet_grid(. ~ Comparison, scales = "free_x", space = "free_x") +
      scale_x_discrete(labels = setup_label) +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
      scale_fill_manual(values = my_colors) +
      labs(y = ylab, x = NULL, fill = "Samples") +
      paper_theme() +
      theme(strip.background = element_blank(),
            strip.text       = element_text(face = "bold"),
            legend.position  = if (dplyr::n_distinct(dat_summary$Sample) > 1) "top" else "none")
  }

  ## ------------------------------------------------------------------
  ## filtering / annotation
  ## ------------------------------------------------------------------

  if ("is_decoy" %in% colnames(dat)) dat <- dat %>% dplyr::filter(!as.logical(is_decoy))

  cmp <- data.frame(Comparison = rep(names(comparisons), lengths(comparisons)),
                    Setup      = unlist(comparisons, use.names = FALSE),
                    stringsAsFactors = FALSE)
  cmp$pos <- stats::ave(seq_along(cmp$Setup), cmp$Comparison, FUN = seq_along)

  dat <- dat %>%
    dplyr::filter(!Repeatability, timePlex == FALSE, Label == "LF") %>%   # LF no-plex only
    dplyr::mutate(Setup = paste0(stringr::str_match(V9, "([0-9.]+)\\s*nL/min")[, 2], "_",
                                 stringr::str_match(V9, "([0-9]+)\\s*columns?\\s+spraying")[, 2]))

  unknown <- setdiff(cmp$Setup, unique(dat$Setup))
  if (length(unknown) > 0) {
    stop("No runs for setup(s) ", paste(unknown, collapse = ", "),
         ". Setups found in Notes: ", paste(sort(unique(dat$Setup)), collapse = ", "))
  }

  dat <- dat %>%
    dplyr::inner_join(cmp, by = "Setup", relationship = "many-to-many") %>%
    dplyr::mutate(Comparison = factor(Comparison, levels = names(comparisons)),
                  Setup      = factor(Setup, levels = unique(cmp$Setup)))

  ## ------------------------------------------------------------------
  ## precursor panel (a)
  ## ------------------------------------------------------------------

  dat_counts_prec <- dat %>% dplyr::count(Comparison, Setup, Sample, Run, name = "n")
  summary_prec    <- summarise_counts(dat_counts_prec)

  p_prec <- build_plot(summary_prec, "# Precursors / Run")

  ## ------------------------------------------------------------------
  ## protein panel (b)
  ## ------------------------------------------------------------------

  dat_counts_prot <- dat %>%
    dplyr::filter(Protein_Qvalue < 0.01) %>%
    dplyr::distinct(Comparison, Setup, Sample, Run, prot) %>%
    dplyr::count(Comparison, Setup, Sample, Run, name = "n")
  summary_prot    <- summarise_counts(dat_counts_prot)

  p_prot <- build_plot(summary_prot, "# Proteins / Run")

  ## ------------------------------------------------------------------
  ## MS1 intensity change, first / second setup of each pair (c)
  ## ------------------------------------------------------------------

  quant <- dat %>%
    dplyr::filter(MS1_Area > 1) %>%
    dplyr::group_by(Comparison, pos, Sample, seqcharge) %>%
    dplyr::summarise(quant = mean(MS1_Area), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pos, values_from = quant, names_prefix = "quant_") %>%
    dplyr::filter(!is.na(quant_1), !is.na(quant_2)) %>%
    dplyr::mutate(ratio      = quant_1 / quant_2,
                  log2_ratio = log2(ratio))

  quant_summary <- quant %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(Numerator   = comparisons[[as.character(Comparison[1])]][1],
                     Denominator = comparisons[[as.character(Comparison[1])]][2],
                     n_prec      = dplyr::n(),
                     median_log2 = stats::median(log2_ratio),
                     fold        = 2^median_log2,
                     .groups     = "drop") %>%
    dplyr::mutate(label = sprintf("%.2fx\nn = %s", fold, comma(n_prec)))

  ratio_labels <- sapply(comparisons, function(s) paste0(gsub("\n", " ", setup_label(s[1])), " /\n",
                                                         gsub("\n", " ", setup_label(s[2]))))
  y_lim <- as.numeric(stats::quantile(quant$log2_ratio, probs = c(0.005, 0.995)))

  p_int <- ggplot(quant, aes(y = log2_ratio)) +
    geom_histogram(fill = "grey65", color = "black", linewidth = 0.2, binwidth = 0.1) +
    geom_hline(data = quant_summary, aes(yintercept = median_log2),
               color = "red", linetype = "dashed", linewidth = 0.5) +
    geom_text(data = quant_summary, aes(x = Inf, y = Inf, label = label),
              hjust = 1.05, vjust = 1.2, size = paper_geom_text()) +
    facet_grid(. ~ Comparison, labeller = labeller(Comparison = ratio_labels)) +
    coord_cartesian(ylim = y_lim) +
    scale_x_continuous(labels = comma, n.breaks = 3) +
    labs(x = "# Precursors", y = bquote(Log[2]*", MS1 intensity ratio")) +
    paper_theme() +
    theme(strip.background = element_blank(),
          strip.text       = element_text(size = paper_size("text"), face = "bold"),
          panel.spacing.x  = unit(1, "lines"))

  ## ------------------------------------------------------------------
  ## console output
  ## ------------------------------------------------------------------

  print(summary_prec %>% dplyr::select(Comparison, Setup, Sample, n_rep, mean_count, se) %>%
          dplyr::rename(precursors = mean_count), n = Inf)
  print(summary_prot %>% dplyr::select(Comparison, Setup, Sample, n_rep, mean_count, se) %>%
          dplyr::rename(proteins = mean_count), n = Inf)
  print(quant_summary %>% dplyr::select(-label))

  list(ed2a = paper_panel(p_prec, file.path(out_dir, "ExtDataFig2a.pdf"), width = 3.64, height = 2.12),
       ed2b = paper_panel(p_prot, file.path(out_dir, "ExtDataFig2b.pdf"), width = 3.64, height = 2.12),
       ed2c = paper_panel(p_int,  file.path(out_dir, "ExtDataFig2c.pdf"), width = 3.45, height = 2.9),
       prec = dat_counts_prec, prot = dat_counts_prot,
       intensity = quant, intensity_summary = quant_summary)
}


## ---- Optimized vs unoptimized spectral library (Extended Data Fig. 4) -------

compare_coverage <- function(mTRAQ, LF, out_dir = "Figs/ExtDataFig4"){
  paper_figure(out_dir)

  mTRAQ <- mTRAQ %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))

  LF <- LF %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))

  mTRAQcounts <- mTRAQ %>% dplyr::filter(is_decoy==F & Qvalue <0.5 & Best_pDChannel_Qvalue<0.01) %>%
    dplyr::add_count(Run) %>% dplyr::select("Run","n","Type") %>% dplyr::distinct(Run,.keep_all=T)

  LFcounts <- LF %>% dplyr::filter(is_decoy==F & Qvalue <0.5 & Best_pDChannel_Qvalue<0.01) %>%
    dplyr::add_count(Run) %>% dplyr::select("Run","n") %>% dplyr::distinct(Run,.keep_all=T) %>%
    dplyr::rename("LF_n"="n")

  joined <- mTRAQcounts %>% inner_join(LFcounts, by = "Run")

  joined$change <- (joined$n/joined$LF_n - 1)*100   # optimized vs unoptimized library, % of the unoptimized count

  summary_stats <- joined %>%
    dplyr::group_by(Type) %>%
    dplyr::summarise(mean_change = mean(change),
              se = sd(change) / sqrt(n()))

  summary_stats <- summary_stats %>%
    dplyr::mutate(label   = sprintf("%+.1f%%", mean_change),
                  label_y = ifelse(mean_change >= 0, mean_change + se, mean_change - se),
                  vjust   = ifelse(mean_change >= 0, -0.6, 1.6))

  plex_labs <- c("3-plexDIA" = "3-plexDIA\n\n(3-plex)", "3-plexDIA\n3-timePlex" = "3-plexDIA\n3-timePlex\n(9-plex)")

  p_change <- ggplot(joined, aes(x = Type)) +
    geom_col(data = summary_stats, aes(y = mean_change),
             width = 0.6, fill = "gray70", color = "black", linewidth = 0.25) +
    geom_errorbar(data = summary_stats, aes(ymin = mean_change - se, ymax = mean_change + se),
                  width = 0.2, linewidth = 0.25) +
    geom_text(data = summary_stats, aes(y = label_y, label = label, vjust = vjust),
              size = paper_geom_text()) +
    scale_x_discrete(labels = plex_labs) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = NULL, y = "Percent change in # precursors\nidentified with optimized library") +
    paper_theme() +
    theme(axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12))


  LFcounts <- LF %>% dplyr::filter(is_decoy==F & Qvalue <0.5 & Best_pDChannel_Qvalue<0.01) %>%
    dplyr::add_count(Run) %>% dplyr::select("Run","n","Type") %>% dplyr::distinct(Run,.keep_all=T)

  LFcounts$lib <- "Unoptimized library"
  mTRAQcounts$lib <- "Optimized library"

  both <- rbind(LFcounts, mTRAQcounts)
  both$lib <- factor(both$lib, levels = c("Unoptimized library", "Optimized library"))


  summary_stats_both <- both %>% dplyr::group_by(Type,lib) %>% dplyr::summarise(mean_n = mean(n),
                                                                  se = sd(n) / sqrt(n()))

  pd <- position_dodge(width = 0.6)

  p_counts <- ggplot(both, aes(x = Type, fill = lib)) +
    geom_col(data = summary_stats_both, aes(y = mean_n),
             position = pd, width = 0.6, color = "black", alpha = 0.7, linewidth = 0.25) +
    geom_errorbar(data = summary_stats_both, aes(ymin = mean_n - se, ymax = mean_n + se),
                  position = pd, width = 0.2, linewidth = 0.25) +
    geom_beeswarm(aes(y = n, fill = lib), shape = 21, color = "black", stroke = 0.3,
                  dodge.width = 0.6, size = 1.4, alpha = 0.5, cex = 2, show.legend = FALSE) +
    scale_x_discrete(labels = plex_labs) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("Unoptimized library" = "#F8766D", "Optimized library" = "#00BFC4"),
                      name = "Data searched with:") +
    labs(x = NULL, y = "Precursor data points / run") +
    paper_theme() +
    # legend as in the manuscript: above a, left-aligned, large keys
    theme(axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12),
          legend.position = "top", legend.justification = "left", legend.direction = "vertical",
          legend.title.position = "top", legend.key.size = unit(14, "pt"),
          legend.margin = margin(0, 0, 4, 20))

  invisible(list(ed4a = paper_panel(p_counts, file.path(out_dir, "ExtDataFig4a.pdf"), width = 3.545, height = 3.9),
                 ed4b = paper_panel(p_change, file.path(out_dir, "ExtDataFig4b.pdf"), width = 3.545, height = 3.3),
                 change = summary_stats, counts = summary_stats_both))
}


## ---- Coverage per column (Extended Data Fig. 5b-c) --------------------------

columns_specific_coverage <- function(dat, use_best_channel = F, filter_level=0.01, extra="IO",
                                      title = NULL, out_file = paste0("points_tp_perColumnperSample", extra, ".pdf")){
  paper_figure(out_file)


  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(Best_pDChannel_Qvalue<0.01) %>% dplyr::filter(Qvalue<filter_level)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<0.01)
  }

  dat <- dat %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))
  dat <- dat[dat$Repeatability ==F,]


  dat <- dat %>% dplyr::mutate(Label = ifelse(Type == "3-plexDIA","mTRAQ",Label))
  dat <- dat[dat$Type=="LF\n3-timePlex",]



  dat_counts <- dat %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    dplyr::ungroup()

  dat_counts$Column <- as.factor(as.numeric(dat_counts$Column))
  # LF timePlex samples C1/C2/C3 are samples A/B/C
  dat_counts$Sample <- ifelse(grepl("C1|D1", dat_counts$Sample), "A",
                              ifelse(grepl("C2|D2", dat_counts$Sample), "B",
                                     ifelse(grepl("C3|D3", dat_counts$Sample), "C", dat_counts$Sample)))
  colors = c("#edae49", "#d1495b","#00798c")
  p <- ggplot(dat_counts, aes(x=Column, y=n,fill=Sample)) +
    geom_point(alpha=0.5, size=2.6, shape=21, stroke = 0.3) +
    facet_wrap(~Sample) +
    scale_fill_manual(values=colors) +
    scale_y_continuous(labels = scales::comma, limits = c(0, max(dat_counts$n)+100)) +
    labs(x = "LC column", y = "# precursors", title = title) +
    paper_theme_box() +
    theme(legend.position = "none", plot.title = element_text(size = 10))

  invisible(list(panel = paper_panel(p, out_file, width = 2.363, height = 2.85), counts = dat_counts))
}


## ---- 9-plexDIA x 3-timePlex (Fig. 5c) ---------------------------------------

tp_coverage_3x9 <- function(dat,extra="", use_best_channel = F, filter_level=0.01, plex3x9=F,
                            out_files = c(prec = "Figs/Fig5/Fig5c_precursors.pdf",
                                          prot = "Figs/Fig5/Fig5c_proteins.pdf")){
  paper_figure(out_files)

  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(Best_pDChannel_Qvalue<0.01) %>% dplyr::filter(Qvalue < filter_level)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<0.01)
  }

  if(plex3x9==T){
    dat <- dat %>%
      dplyr::mutate(Type = case_when(
        timePlex == TRUE & Label != "LF" ~ "9-plexDIA\n3-timePlex",
        timePlex == FALSE & Label != "LF" ~ "9-plexDIA",
        TRUE ~ "LF\nNo-plex"
      ))
  } else{
    dat <- dat %>%
      dplyr::mutate(Type = case_when(
        timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
        timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
        timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
        TRUE ~ "LF\nNo-plex"
      ))
  }

  dat <- dat[dat$Repeatability ==F,]
  dat_keep <- dat



  dat$Label <- paste0(dat$Label,"_",dat$time_channel)

  # Samples J/K/L/M of the 9-plex are shown as A/B/C/D
  sample_levels <- if (plex3x9) c("J", "K", "L", "M") else c("A", "B", "C")
  sample_labels <- if (plex3x9) c("A", "B", "C", "D") else c("A", "B", "C")
  my_colors     <- unname(paper_cols$samples4)

  # Stacked bars, one segment per sample x mass-tag x timePlex channel, with
  # the stack total on top
  summarise_stack <- function(dat_counts, text_offset){
    dat_summary <- dat_counts %>%
      dplyr::group_by(Sample, Label, Amount, Type) %>%
      dplyr::summarise(
        mean_count = mean(n, na.rm = TRUE),
        se = sd(n, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>% dplyr::ungroup()

    dat_summary <- dat_summary %>%
      dplyr::arrange(desc(Sample), Type) %>%
      dplyr::group_by(Type) %>%
      dplyr::mutate(
        stacked_mean = cumsum(mean_count)) %>% dplyr::ungroup() %>%
      dplyr::group_by(Label, Sample, Type) %>%
      dplyr::mutate(stacked_mean = ifelse(Type=="LF\nNo-plex", mean_count,stacked_mean))%>%
      dplyr::mutate(
        er_max = stacked_mean + se,
        er_min = stacked_mean - se
      ) %>%
      dplyr::ungroup()

    df_text <- dat_summary %>%
      dplyr::group_by(Type, Sample,Label) %>%
      dplyr::summarise(
        text_value = round(max(er_max) - max(se)),
        text_pos = max(er_max) + text_offset
      )

    df_text_max <- df_text %>%
      dplyr::group_by(Type) %>%
      dplyr::mutate(label = ifelse(text_value == max(text_value, na.rm = TRUE),
                                   comma(text_value), "")) %>%
      dplyr::ungroup()
    list(summary = dat_summary, text = df_text_max)
  }

  stack_plot <- function(s, ylab, legend){
    stacked_types <- c("9-plexDIA\n3-timePlex", "9-plexDIA")
    ggplot(s$summary %>% dplyr::filter(!grepl("No-plex",Type)), aes(x = Type, y = mean_count, fill = Sample)) +
      scale_x_discrete(limits = c("9-plexDIA", "9-plexDIA\n3-timePlex")) +
      geom_bar(data = s$summary %>% dplyr::filter(Type %in% stacked_types),
               stat = "identity", color = "black", width = 0.35, linewidth = 0.25) +
      geom_errorbar(data = s$summary %>% dplyr::filter(Type %in% stacked_types),
                    aes(ymin = er_min, ymax = er_max), width = 0.18, linewidth = 0.25) +
      geom_text(data = s$text %>% dplyr::filter(Type %in% stacked_types),
                aes(x = Type, y = text_pos, label = label), inherit.aes = FALSE,
                vjust = 0, size = paper_geom_text(5)) +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.06))) +
      scale_fill_manual(values = my_colors, labels = sample_labels, drop = FALSE) +
      labs(y = ylab, x = NULL, fill = "Samples") +
      paper_theme() +
      theme(legend.position        = if (legend) "inside" else "none",
            legend.position.inside = c(0.03, 1),
            legend.justification   = c(0, 1),
            legend.direction       = "horizontal",
            legend.title           = element_text(size = 6),
            legend.text            = element_text(size = 6),
            legend.key.size        = unit(7, "pt"),
            legend.key.spacing.x   = unit(2, "pt")) +
      guides(fill = guide_legend(title = "Samples", label.position = "bottom",
                                 title.position = "top", title.hjust = 0.5))
  }

  dat_counts <- dat %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    dplyr::ungroup()
  dat_counts$Sample <- factor(dat_counts$Sample, levels = sample_levels)

  p_prec <- stack_plot(summarise_stack(dat_counts, 0), "Precursors data points / Run", legend = TRUE)

  dat_counts <- dat %>% dplyr::filter(Protein_Qvalue<0.01) %>% dplyr::distinct(Sample, Run, Column, Amount, Label, Type, prot,.keep_all=T) %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>% dplyr::ungroup()
  dat_counts$Sample <- factor(dat_counts$Sample, levels = sample_levels)

  p_prot <- stack_plot(summarise_stack(dat_counts, 0), "Protein data points / Run", legend = FALSE)

  invisible(list(prec   = paper_panel(p_prec, out_files[["prec"]], width = 1.6, height = 2.2),
                 prot   = paper_panel(p_prot, out_files[["prot"]], width = 1.52, height = 2.2),
                 counts = dat_counts))
}


## ---- Dual-timePlex IDs per channel and set, FWHM per column (Extended Data Fig. 10) ----

prec_prot_counts <- function(dat_all_DualTP, runs=NULL,
                             tp_labels = c(T0 = "\u03940 min", T1 = "\u03942.5 min", T2 = "\u03945 min"),
                             out_dir = "Figs/ExtDataFig10"){
  paper_figure(out_dir)
  ## ---- settings (edit to taste) ----
  prot_q  <- 0.01   # Protein_Qvalue
  chan_q  <- 0.01   # Best_tP_specific_channelQvalue
  prec_q  <- 0.5      # Qvalue (set e.g. 0.5 to also cap run-level precursor q)... JMod outputs at Qvalue < 0.01 for LF no plex 
  n_plex  <- 9      # plexDIA channels per timePlex channel
  n_cols  <- 6      # LC columns (2 Sets x 3 timePlex channels); FWHM precursors must be in all

  if (!is.null(runs)) dat_all_DualTP <- dat_all_DualTP %>% dplyr::filter(Run %in% runs)

  ## ---- confident IDs, annotated with Set / tag / LC column ----
  # Set = A if Run number is odd, B if even
  dat_all_DualTP$untag_prec <- dat_all_DualTP$seqcharge
  dat_all_DualTP$protein <- dat_all_DualTP$prot
  dat_filt <- dat_all_DualTP %>%
    dplyr::filter(!is_decoy,
           Protein_Qvalue < prot_q,
           Best_pDChannel_Qvalue < chan_q,
           Qvalue < prec_q) %>%
    dplyr::mutate(run_num = as.integer(gsub("\\D", "", Run)),
           Set     = ifelse(run_num %% 2 == 1, "Set A", "Set B"),
           tag     = sub("^.*_", "", channel),
           tP      = paste0("T", time_channel))

  ## ---- IDs per sample for each Run x timePlex channel ----
  # total channel-level IDs across all mass tags / n_plex = average IDs per sample
  ids_per_run <- dat_filt %>%
    dplyr::group_by(Run, Set, tP) %>%
    dplyr::summarise(Precursors = n_distinct(untag_prec, tag) / n_plex,
              Proteins   = n_distinct(protein, tag) / n_plex, .groups = "drop") %>%
    pivot_longer(c(Precursors, Proteins), names_to = "level", values_to = "n")

  ## ---- average across runs for each LC column (Set x timePlex channel) ----
  ids_avg <- ids_per_run %>%
    dplyr::group_by(level, Set, tP) %>%
    dplyr::summarise(mean = mean(n), sd = sd(n), n_runs = n(), .groups = "drop")

  ## --- mean per sample for Set A vs Set B (across runs and timePlex channels) ---
  print(ids_per_run %>% dplyr::group_by(level, Set) %>%
          dplyr::summarise(mean = mean(n), sd = sd(n), .groups = "drop"))

  ## ---- bar plots: b = precursors, c = proteins ----
  set_cols <- c("Set A" = "#F8766D", "Set B" = "#00BFC4")
  plot_bar <- function(lvl, legend) {
    dodge <- position_dodge(width = 0.8)
    ggplot(dplyr::filter(ids_avg, level == lvl), aes(tP, mean, fill = Set)) +
      geom_col(position = dodge, width = 0.75, alpha = 0.8) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                    position = dodge, width = 0.25, linewidth = 0.3) +
      scale_x_discrete(labels = tp_labels) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
      scale_fill_manual(values = set_cols, name = NULL) +
      labs(x = "timePlex channel", y = paste("#", lvl, "/ sample")) +
      # the legend is drawn invisibly in c so that b and c keep the same panel height
      guides(fill = guide_legend(override.aes = list(alpha = if (legend) 0.8 else 0))) +
      paper_theme_box(grid = TRUE) +
      theme(legend.position = "top", legend.key.size = unit(10, "pt"),
            legend.text = element_text(colour = if (legend) "black" else NA),
            plot.margin = margin(3, 3, 3, 3))
  }

  ## ---- FWHM: median per precursor per LC column, intersected across all columns ----
  fwhm_col <- dat_filt %>%
    dplyr::filter(!is.na(MS1_FWHM)) %>%
    dplyr::group_by(untag_prec, Set, tP) %>%
    dplyr::summarise(med_FWHM = median(MS1_FWHM) * 60, .groups = "drop")   # min -> s

  n_found <- nrow(distinct(fwhm_col, Set, tP))
  if (n_found != n_cols)
    warning(sprintf("Found %d LC columns with FWHM values, expected %d", n_found, n_cols))

  med_fwhm <- fwhm_col %>%
    dplyr::add_count(untag_prec) %>%
    dplyr::filter(n == n_cols)
  message(sprintf("FWHM: %d precursors present in all %d LC columns",
                  n_distinct(med_fwhm$untag_prec), n_cols))
  x_lim <- quantile(med_fwhm$med_FWHM, probs = c(0, 0.995), names = FALSE)
  fwhm_lab <- med_fwhm %>%
    dplyr::group_by(Set, tP) %>%
    dplyr::summarise(med = signif(median(med_FWHM), 2), .groups = "drop") %>%
    dplyr::mutate(lab = sprintf(paste0("%.", pmax(0, 1 - floor(log10(abs(med)))), "f s"), med))

  p_fwhm <- ggplot(med_fwhm, aes(x = med_FWHM, y = tP)) +
    facet_grid(rows = vars(Set)) +
    ggridges::geom_density_ridges2(
      alpha = 0.75,
      scale = 1.0,
      color = "black",
      fill  = "grey70",
      linewidth = 0.25
    ) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.33,
      color = "black",
      linewidth = 0.25,
      middle.linewidth = 0.5
    ) +
    geom_text(data = fwhm_lab, inherit.aes = FALSE,
              aes(x = x_lim[2], y = tP, label = lab),
              hjust = 1, vjust = -0.4, size = paper_geom_text(paper_size("text") - 1)) +
    scale_y_discrete(labels = tp_labels) +
    labs(x = "FWHM (seconds)", y = "timePlex channel") +
    coord_cartesian(xlim = x_lim) +
    paper_theme_box() +
    theme(legend.position = "none")

  invisible(list(ed10a = paper_panel(p_fwhm, file.path(out_dir, "ExtDataFig10a.pdf"), width = 2.0, height = 2.5),
                 ed10b = paper_panel(plot_bar("Precursors", legend = TRUE),
                                     file.path(out_dir, "ExtDataFig10b.pdf"), width = 2.6, height = 2.5),
                 ed10c = paper_panel(plot_bar("Proteins", legend = FALSE),
                                     file.path(out_dir, "ExtDataFig10c.pdf"), width = 2.49, height = 2.5),
                 ids = ids_avg, fwhm = med_fwhm))
}
