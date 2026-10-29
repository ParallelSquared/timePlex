################################################################################
# Fig. 4 v2: c-e with the samples per run of each method in the facet headers
################################################################################

# Samples measured per run for each data type in Fig. 4c-e
fig4_samples_per_run <- c(np = 1, pd = 3, tp = 3, tp_pd = 9)

# Intersected precursor ratios per comparison, as in quant_all_comparisons_v2()
quant_comparisons_data <- function(normed_quant, grouper = "seqcharge", min_num_reps = 2, alphascaler = 1){

  normed_quant <- normed_quant %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tp_pd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))

  normed_quant <- normed_quant[normed_quant$Repeatability==F,] #remove Runs that are just for computing reproducibility
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",ProteinName),"Human",
                                                                    ifelse(grepl("YEAST",ProteinName),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))

  normed_quant <- normed_quant %>% dplyr::mutate("Sample" = ifelse(grepl("C1|D1",Sample),"A",
                                                                   ifelse(grepl("C2|D2",Sample),"B","C")))

  np <- normed_quant[normed_quant$Type=="np",]

  #take middle replicate for all so that number of samples matches number of samples for LF no-plex
  pD <- normed_quant[normed_quant$Type=="pd",] %>% dplyr::filter(grepl("JD0425|JD0426|JD0427",run_chan))
  tp_pd <- normed_quant[normed_quant$Type=="tp_pd",] %>% dplyr::filter(Rep == 2)
  tp <- normed_quant[normed_quant$Type=="tp",] %>% dplyr::filter(Rep == 2)

  joined <- rbind(np, pD, tp_pd, tp)
  normed_quant_avg_ints <- joined %>% dplyr::group_by(Type, .data[[grouper]], Sample) %>%
    dplyr::summarize("quant_avg" = mean(quant_val,na.rm=T)) %>% dplyr::rename("denom_sample"="Sample")

  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()

  joined <- joined %>%
    left_join(human_medians, by = "run_chan") %>%
    dplyr::mutate(norm = norm-median_Human) %>%
    dplyr::select(-median_Human)

  joined_quant <- joined %>% dplyr::group_by(Type, Sample,.data[[grouper]],species) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_num_reps-1)) %>% #require atleast 1/3
    dplyr::summarize(med_quant = median(norm,na.rm=T))

  formula_str <- paste(grouper, "species", sep = "+")
  joined_quant_d <- dcast(
    joined_quant,
    formula = as.formula(paste(formula_str, "~ Type+Sample")),
    value.var = "med_quant"
  )
  for (ty in c("np", "tp", "pd", "tp_pd")) {
    joined_quant_d[[paste0(ty, "_AB")]] <- joined_quant_d[[paste0(ty, "_A")]] - joined_quant_d[[paste0(ty, "_B")]]
    joined_quant_d[[paste0(ty, "_AC")]] <- joined_quant_d[[paste0(ty, "_A")]] - joined_quant_d[[paste0(ty, "_C")]]
    joined_quant_d[[paste0(ty, "_BC")]] <- joined_quant_d[[paste0(ty, "_B")]] - joined_quant_d[[paste0(ty, "_C")]]
  }

  joined_quant_m <- joined_quant_d %>% dplyr::select(.data[[grouper]],species,
                                                     np_AB, np_AC, np_BC,
                                                     tp_AB, tp_AC, tp_BC,
                                                     pd_AB, pd_AC, pd_BC,
                                                     tp_pd_AB, tp_pd_AC, tp_pd_BC) %>% melt()

  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()

  joined_quant_m <- joined_quant_m %>%
    left_join(human_medians, by = "variable") %>%
    dplyr::mutate(value = value-median_Human) %>%
    dplyr::select(-median_Human) %>% na.omit()

  joined_quant_m <- joined_quant_m %>%
    dplyr::mutate(denom_sample = ifelse(grepl("AB",variable),"B","C"))
  joined_quant_m$Type <- sub("_[^_]*$", "", joined_quant_m$variable)

  joined_quant_m <- joined_quant_m %>% left_join(normed_quant_avg_ints, by =c("Type","denom_sample",paste0(grouper)))
  allcounts <- joined_quant_m %>% na.omit() %>% dplyr::count(variable)

  comparisons <- list(AB = list(num = "A", den = "B", exp_FC = -2, panel = "c"),
                      AC = list(num = "A", den = "C", exp_FC = 1,  panel = "d"),
                      BC = list(num = "B", den = "C", exp_FC = 3,  panel = "e"))
  rows <- lapply(names(comparisons), function(cmp){
    joined_cmp <- joined_quant_m[grepl(cmp, joined_quant_m$variable),]
    int_cmp <- joined_cmp %>% na.omit() %>% dplyr::count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
    joined_cmp_int <- joined_cmp %>% dplyr::filter(.data[[grouper]]%in%int_cmp[[grouper]]) %>% na.omit()
    joined_cmp_int$type <- factor(sub(paste0("_", cmp, "$"), "", joined_cmp_int$variable), levels = names(fig4_samples_per_run))
    joined_cmp_int$Alpha <- ifelse(grepl("Human", joined_cmp_int$species), 0.04, 0.12) * alphascaler
    list(int = joined_cmp_int, n = nrow(int_cmp),
         x_lim = as.numeric(quantile(sort(log2(joined_cmp_int$quant_avg)), probs = c(.005, .995))),
         y_lim = as.numeric(quantile(sort(joined_cmp_int$value), probs = c(.005, .995))))
  })
  names(rows) <- names(comparisons)

  # precursor ratios quantified per data type, and intersected across all four (the bars of
  # the original Fig. 4c-e)
  counts <- do.call(rbind, lapply(names(comparisons), function(cmp){
    cc <- allcounts[grepl(paste0("_", cmp, "$"), allcounts$variable), ]
    data.frame(comparison = cmp, type = c(sub(paste0("_", cmp, "$"), "", cc$variable), "Intersected"),
               n = c(cc$n, rows[[cmp]]$n))
  }))
  list(rows = rows, comparisons = comparisons, counts = counts)
}

# Fig. 4c-e v2: per comparison the scatter plots of the four data types (samples
# per run in the facet headers) and their box plots
quant_all_comparisons_v3 <- function(normed_quant, grouper = "seqcharge", alphascaler = 1, min_num_reps = 2,
                                     out_dir = "Figs/Fig4"){
  d <- quant_comparisons_data(normed_quant, grouper, min_num_reps, alphascaler)
  rows <- d$rows; comparisons <- d$comparisons
  n_run <- fig4_samples_per_run
  type_labs <- stats::setNames(paste0(c("LF no-plex", "3-plexDIA", "3-timePlex", "3x3 plexDIA & timePlex"), "\n",
                                      n_run, ifelse(n_run == 1, " sample/run", " samples/run")), names(n_run))
  short_labs <- c(np = "LF np", pd = "3-pD", tp = "3-tP", tp_pd = "3x3")
  sp <- paper_cols$species[c("Human", "Yeast")]
  x_lim <- range(sapply(rows, `[[`, "x_lim"))
  n_pt <- round(0.85 * paper_size("title"))   # n line 15% below the y title (pdf() prints whole points)

  ce_grid <- function(letters = FALSE){
    tag <- function(lab, vjust = 1) ggplot() +
      annotate("text", x = 0, y = 1, label = lab, hjust = 0, vjust = vjust,
               fontface = "bold", size = paper_geom_text(paper_size("tag"))) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
      theme_void()

    plots <- list()
    for (k in seq_along(comparisons)){
      cmp <- names(comparisons)[k]; cc <- comparisons[[cmp]]; r <- rows[[cmp]]
      first <- k == 1; last <- k == length(comparisons)
      ratio_lab <- paste0(cc$num, "/", cc$den)

      scatter <- ggplot(r$int, aes(x = log2(quant_avg), y = value, colour = species, alpha = Alpha)) +
        expected_lines(cc$exp_FC) +
        paper_raster(geom_point(size = 0.8, stroke = 0)) +
        facet_grid(cols = vars(type), labeller = as_labeller(type_labs)) +
        scale_colour_manual(values = sp) +
        scale_alpha_identity() +
        scale_x_continuous(breaks = paper_int_breaks(2), labels = paper_int_labels) +
        scale_y_continuous(breaks = paper_int_breaks(1)) +
        coord_cartesian(xlim = x_lim, ylim = r$y_lim) +
        # two tight lines as markdown (plotmath atop() leaves a wide gap that squeezes the
        # panels)
        labs(x = bquote(Log[2]*", denominator"),
             y = paste0("Log<sub>2</sub>, ", ratio_lab, "<br><span style='font-size:", n_pt, "pt'>*n* = ",
                        scales::comma(r$n), "</span>")) +
        paper_theme_box() +
        theme(legend.position = "none", strip.text = element_text(size = paper_size("text")),
              axis.title.y = ggtext::element_markdown(lineheight = 1, margin = margin(r = 1)),
              panel.spacing.x = unit(2, "pt"), plot.margin = margin(1, 0, 1, 0))

      box <- ggplot(r$int, aes(x = type, y = value, fill = species)) +
        expected_lines(cc$exp_FC) +
        geom_boxplot(alpha = 0.85, outlier.shape = NA, linewidth = 0.25, show.legend = FALSE) +
        paper_box_outliers(alpha = alphascaler * 0.01) +
        scale_fill_manual(values = sp, name = NULL) +
        scale_x_discrete(labels = short_labs) +
        scale_y_continuous(breaks = paper_int_breaks(1)) +
        coord_cartesian(ylim = r$y_lim) +
        paper_theme_box() +
        theme(axis.title = element_blank(), axis.text.y = element_blank(), legend.position = "none",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              plot.margin = margin(1, 1, 1, 2))
      # species key above the first box plots, beside the facet headers
      if (first) box <- box +
        geom_point(data = data.frame(type = factor("np", levels = names(n_run)), value = r$y_lim[1], species = names(sp)),
                   shape = 22, size = 2.6, stroke = 0.3, alpha = 0) +
        guides(fill = species_guide()) +
        theme(legend.position = "inside", legend.position.inside = c(0.5, 1), legend.justification.inside = c(0.5, 0),
              legend.text = element_text(size = paper_size("text")), legend.key.size = unit(1.1 * paper_size("text"), "pt"))

      if (!first) scatter <- scatter + theme(strip.text = element_blank(), strip.background = element_blank())
      if (!last) {
        scatter <- scatter + paper_drop_x()
        box <- box + paper_drop_x()
      } else {
        box <- patchwork::free(box, type = "space", side = "b")
      }
      if (letters) plots <- c(plots, list(tag(if (first) " " else cc$panel)))   # c is placed level with f by paper_fig4_v2()
      plots <- c(plots, list(scatter, box))
    }
    widths <- if (letters) grid::unit(c(0.16, 4.05, 0.62), c("in", "null", "null")) else c(4.05, 0.62)
    # e ends level with f because the block has no bottom margin
    patchwork::wrap_plots(plots, ncol = length(widths), widths = widths) +
      patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 1, 0, 0)))
  }

  W <- 5.6; H <- 2.78   # ratio panels 0.74 in tall each
  ce_file <- file.path(out_dir, "Fig4c-e_v2.pdf")
  paper_save(ce_grid(), ce_file, width = W - 0.16, height = H)
  fig4ce <- list(grob = paper_grob(ce_grid(letters = TRUE)), width = W, height = H, file = ce_file)
  invisible(list(fig4ce = fig4ce, counts = d$counts))
}

# Extended Data Fig. 7: precursor ratios quantified per data type for Fig. 4c-e
# (A/B, A/C, B/C), laid out as the preprint's ratio-count figure
paper_ed7 <- function(counts, file = "Figs/ExtDataFig7/ExtDataFig7_full.pdf"){
  paper_figure(file)
  type_labs <- c(np = "LF no-plex", pd = "3-plexDIA", tp = "3-timePlex", tp_pd = "3x3 plexDIA & timePlex",
                 Intersected = "Intersected")
  type_cols <- c(np = "grey10", pd = "grey30", tp = "grey50", tp_pd = "grey70", Intersected = "grey90")   # as the original Fig. 4c-e bars
  panels <- lapply(c(AB = "A/B", AC = "A/C", BC = "B/C"), function(lab){
    d <- counts[counts$comparison == gsub("/", "", lab), ]
    d$type <- factor(d$type, levels = names(type_labs))
    p <- ggplot(d, aes(x = type, y = n, fill = type)) +
      geom_col(colour = "black", alpha = 0.7, linewidth = 0.25, width = 0.85) +
      geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = paper_geom_text(paper_size("text") - 1)) +
      scale_fill_manual(values = type_cols, guide = "none") +
      scale_x_discrete(labels = type_labs) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
      labs(x = NULL, y = paste0("Precursor ratios (", lab, ")")) +
      paper_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.margin = margin(16, 4, 2, 2))      # room for the panel letter above the y title
    p
  })
  W <- 7.09; w <- W / 3; h <- 2.7
  panels <- mapply(function(p, l) paper_panel(p, sub("_full\\.pdf$", paste0(l, ".pdf"), file), width = w, height = h),
                   panels, c("a", "b", "c"), SIMPLIFY = FALSE)
  xs <- c(0, w, 2 * w)
  paper_assemble(file, width = W, height = h, items = mapply(function(p, x, l) paper_place(p, x, 0, l, x, 0),
                                                             panels, xs, c("a", "b", "c"), SIMPLIFY = FALSE))
}

# Fig. 4 v2 assembly: as paper_fig4(), with f as tall as c-e
paper_fig4_v2 <- function(a, b, ce, f, file = "Figs/Fig4/Fig4_full_v2.pdf"){
  W <- 7.09
  f$height <- ce$height
  a <- paper_align_left(a, paper_plot_left(ce))   # a's plots start where c-e's do
  paper_assemble(file, width = W, height = a$height + ce$height, items = list(
    paper_place(a,  0,             0,        "a", 0, 0),
    paper_place(b,  W - b$width,   0,        "b", W - b$width, 0),
    paper_place(ce, 0,             a$height, "c", 0, a$height),
    paper_place(f,  W - f$width,   a$height, "f", W - f$width, a$height)))
}
