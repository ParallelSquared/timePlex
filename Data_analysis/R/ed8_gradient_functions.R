################################################################################
# Extended Data Fig. 8: quantitative accuracy across the gradient for one
# 3-plexDIA x 3-timePlex run, per column (timePlex channel)
################################################################################

ed8_columns <- c("0" = "1 (Δ0 min)", "1" = "2 (Δ4 min)", "2" = "3 (Δ8 min)")   # time_channel -> column (IO offsets)
ed8_column_fills <- c("0" = "#1E4D20", "1" = "#4E9C48", "2" = "#C5E2BD")                     # dark to light green, as the preprint


## ---- Data ---------------------------------------------------------------------

# Precursor ratios numerator/denominator within each plexDIA set (one per column) of
# `run`: log2 difference of the normalized quant (norm), centered on the human median of
# that column as in the Fig. 4 ratio plots, with the precursor's RT (median over channels)
ed8_ratios <- function(normed_quant, run = "JD0413", numerator = "B", denominator = "C"){
  d <- normed_quant[normed_quant$Run == run & !normed_quant$Repeatability, ]
  d$Sample  <- ifelse(grepl("C1|D1", d$Sample), "A", ifelse(grepl("C2|D2", d$Sample), "B", "C"))
  d$species <- ifelse(grepl("HUMAN", d$ProteinName), "Human", ifelse(grepl("YEAST", d$ProteinName), "Yeast", NA))
  d[!is.na(d$species), ] %>%
    dplyr::distinct(seqcharge, run_chan, .keep_all = TRUE) %>%
    dplyr::group_by(time_channel, seqcharge, species) %>%
    dplyr::summarise(num = norm[Sample == numerator][1], den = norm[Sample == denominator][1],
                     rt = median(rt), .groups = "drop") %>%
    dplyr::filter(is.finite(num), is.finite(den)) %>%
    dplyr::mutate(ratio = num - den) %>%
    dplyr::group_by(time_channel) %>%
    dplyr::mutate(ratio = ratio - median(ratio[species == "Human"])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(run = run, comparison = paste0(numerator, "/", denominator))
}

# Interquartile range and median of the ratios along RT, per column and species:
# quantiles in `bin`-min windows holding at least `min_n` ratios, smoothed with loess
ed8_iqr_bands <- function(r, bin = 2, min_n = 20, span = 0.4){
  r %>%
    dplyr::mutate(rt_bin = (floor(rt / bin) + 0.5) * bin) %>%
    dplyr::group_by(time_channel, species, rt_bin) %>%
    dplyr::summarise(n = dplyr::n(), q25 = quantile(ratio, 0.25), q50 = median(ratio),
                     q75 = quantile(ratio, 0.75), .groups = "drop") %>%
    dplyr::filter(n >= min_n) %>%
    dplyr::group_by(time_channel, species) %>%
    dplyr::filter(dplyr::n() >= 5) %>%
    dplyr::mutate(dplyr::across(c(q25, q50, q75),
                                ~ stats::predict(stats::loess(.x ~ rt_bin, weights = n, span = span)))) %>%
    dplyr::ungroup()
}


## ---- Plots --------------------------------------------------------------------

# a: number of precursor ratios along RT, stacked by column. As in the
# manuscript, the column legend sits right of the title, above the histogram
ed8_counts_plot <- function(r, xlim, binwidth = 1){
  r$column <- factor(as.character(r$time_channel), levels = names(ed8_columns))
  ggplot(r, aes(x = rt, fill = column)) +
    geom_histogram(binwidth = binwidth, boundary = 0, position = "stack", colour = NA) +
    scale_fill_manual(values = ed8_column_fills, labels = ed8_columns, name = "Column (timePlex channel)") +
    scale_x_continuous(breaks = seq(0, 200, 20)) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(xlim = xlim) +
    labs(x = NULL, y = "# Precursor ratios",
         title = "3-plexDIA x 3-timePlex",
         subtitle = bquote(italic(n)*" = "*.(scales::comma(nrow(r)))*" precursor ratios")) +
    paper_theme() + paper_drop_x() +
    theme(plot.title = element_text(size = 16, hjust = 0, margin = margin(30, 0, 2, 0)),   # legend title above it
          plot.subtitle = element_text(size = 12),
          legend.position = "inside", legend.position.inside = c(0.6, 1.03), legend.justification = c(0, 0),
          legend.title = element_text(size = 13), legend.key.size = unit(13, "pt"),
          legend.key.spacing.y = unit(3, "pt"))
}

# b: ratios along RT, one facet per column, with the per-species IQR (ribbon), median
# (line) and expected ratio (dashed)
ed8_ratio_rt_plot <- function(r, xlim, exp_FC = 3, ylim = c(-6, 6), alphaScaler = 1.4){
  bands <- ed8_iqr_bands(r)
  strip <- function(x) paste0("Column ", ed8_columns[as.character(x)])
  r$Alpha <- ifelse(r$species == "Human", 0.04, 0.16) * alphaScaler
  sp <- paper_cols$species[c("Human", "Yeast")]
  ggplot(r, aes(x = rt)) +
    paper_raster(geom_point(aes(y = ratio, colour = species, alpha = Alpha), size = 1.5, stroke = 0)) +
    expected_lines(exp_FC) +
    geom_ribbon(data = bands, aes(x = rt_bin, ymin = q25, ymax = q75, fill = species), alpha = 0.75) +
    geom_line(data = bands, aes(x = rt_bin, y = q50, colour = species), linewidth = 0.3) +
    facet_wrap(~ time_channel, ncol = 1, labeller = as_labeller(strip)) +
    scale_colour_manual(values = sp, guide = "none") +
    scale_fill_manual(values = sp, name = NULL, guide = species_guide()) +
    scale_alpha_identity() +
    scale_x_continuous(breaks = seq(0, 200, 20)) +
    scale_y_continuous(breaks = paper_int_breaks(4)) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = "Retention time (min)", y = bquote(Log[2]*", "*.(r$comparison[1]))) +
    paper_theme_box() +
    # species legend inside the top facet (column 1), bottom right, as in the manuscript
    theme(strip.text = element_text(size = 12),
          legend.position = "inside", legend.position.inside = c(0.985, 0.725), legend.justification = c(1, 0),
          legend.key.size = unit(13, "pt"), legend.key.spacing.y = unit(2, "pt"))
}

# Extended Data Fig. 8 (a above b, RT axes aligned) plus each panel on its own;
# portrait, like the manuscript's version, with a shorter a
paper_ed8 <- function(r, file = "Figs/ExtDataFig8/ExtDataFig8_full.pdf", width = 7.09,
                      heights = c(2.1, 5.4)){
  paper_figure(file)
  xlim <- range(r$rt)
  a <- ed8_counts_plot(r, xlim)
  b <- ed8_ratio_rt_plot(r, xlim)
  fig <- patchwork::wrap_plots(a, b, ncol = 1, heights = heights) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = paper_size("tag"), face = "bold"),
          plot.tag.location = "plot")                # letters in the top-left corner
  paper_save(fig, file, width = width, height = sum(heights) + 0.3)
  panel_files <- file.path(dirname(file), c("ExtDataFig8a.pdf", "ExtDataFig8b.pdf"))
  paper_save(a, panel_files[1], width = width, height = heights[1])
  paper_save(b, panel_files[2], width = width, height = heights[2])
  c(file, panel_files)
}
