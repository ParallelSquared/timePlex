################################################################################
# Quantification: RT-dependent normalization and quantitative accuracy
################################################################################


## ---- Normalization ----------------------------------------------------------

prec_normalizeQuant <- function(dat, quant_col="MS1_Int", use_best_channel = F,
                                filter_level=0.01,
                                compute_unadjusted_quant = T){
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(Best_pDChannel_Qvalue<0.01) %>% dplyr::filter(Qvalue < filter_level)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<0.01)
  }
  dat <- dat[dat[[quant_col]]>1,]
  #normalization using only HUMAN proteins
  dat
  datH <- dat %>%
    dplyr::filter(!grepl("YEAST", prot) & !grepl("YEAST", ProteinName)) %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(med_H = median(.data[[quant_col]], na.rm = TRUE), .groups = "drop")
  
  dat <- dat %>%
    left_join(datH, by = "run_chan") %>%
    dplyr::mutate(med_norm = log2(.data[[quant_col]]) - log2(med_H)) %>%
    dplyr::group_by(seqcharge) %>%
    dplyr::mutate(seqcharge_mean = mean(.data[[quant_col]], na.rm = TRUE)) %>%
    dplyr::mutate(
                  seqcharge_offset = mean(med_norm, na.rm = TRUE)) %>%
    dplyr::mutate(norm = med_norm - mean(med_norm, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(norm) & norm != 0 & is.finite(norm))
  
  
  dat$run_timechannel <- paste0(dat$Run,"_",dat$time_channel)
  dat$quant_val <- (2^dat$norm)*dat$seqcharge_mean
  
  if(compute_unadjusted_quant ==T){
    
    raw_col <- sub("_rtnorm$", "", quant_col)  # e.g. "quant_rtnorm" -> "quant"
    
    dat <- dat[dat[[raw_col]] > 1, ]
    # normalization using only HUMAN proteins
    datH_raw <- dat %>%
      dplyr::filter(!grepl("YEAST", prot) & !grepl("YEAST", ProteinName)) %>%
      dplyr::group_by(run_chan) %>%
      dplyr::summarise(med_H_raw = median(.data[[raw_col]], na.rm = TRUE), .groups = "drop")
    
    dat <- dat %>%
      left_join(datH_raw, by = "run_chan") %>%
      dplyr::mutate(med_norm_raw = log2(.data[[raw_col]]) - log2(med_H_raw)) %>%
      dplyr::group_by(seqcharge) %>%
      dplyr::mutate(seqcharge_mean_raw   = mean(.data[[raw_col]], na.rm = TRUE),
                    seqcharge_offset_raw = mean(med_norm_raw, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(before_adj = med_norm_raw - seqcharge_offset_raw) %>%
      dplyr::filter(!is.na(before_adj) & before_adj != 0 & is.finite(before_adj))
    
    dat <- dat %>% dplyr::group_by(run_chan) %>%
      dplyr::mutate(rtnorm_factor_MS1_medscaled = log2(median(.data[[raw_col]],na.rm=T)*rtnorm_factor_MS1)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(seqcharge) %>%
      dplyr::mutate(
        rtnorm_factor_MS1_medscaled_offset_raw = mean(rtnorm_factor_MS1_medscaled, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rtnorm_factor_MS1_adj = rtnorm_factor_MS1_medscaled - rtnorm_factor_MS1_medscaled_offset_raw) %>%
      dplyr::filter(!is.na(rtnorm_factor_MS1_adj) & rtnorm_factor_MS1_adj != 0 & is.finite(rtnorm_factor_MS1_adj))
    
  }
  
  return(dat)
}


## ---- Ratio panels (scatter, box, density) shared by the plots below ---------

# Log2 ratio numerator - denominator, centered on the human median, with the
# log2 denominator intensity as `denominator` (scatter x axis)
ratio_data <- function(dat, numerator, denominator){
  dat$FCs <- dat[[paste0(numerator, denominator)]]
  dat$denominator <- dat[[denominator]]
  dat <- dat %>% dplyr::filter(!is.na(FCs) & is.finite(FCs)) %>% dplyr::filter(FCs != 0)
  dat$FCs_corrected <- dat$FCs - median(dat$FCs[grepl("Human", dat$species)])
  dat
}

# Expected log2 ratios as fine, short-dashed lines: 0 for human, exp_FC for
# yeast. pdf() scales dash units with the line width (0.75 pt at linewidth
# 0.35): "31" = 2.2 pt dash, 0.75 pt gap
expected_lines <- function(exp_FC, linewidth = 0.35){
  list(geom_hline(yintercept = exp_FC, linetype = "31", colour = paper_cols$expected[["Yeast"]],
                  linewidth = linewidth),
       geom_hline(yintercept = 0, linetype = "31", colour = paper_cols$expected[["Human"]],
                  linewidth = linewidth))
}

# Square, black-outlined legend keys for the species colors
species_guide <- function(){
  guide_legend(override.aes = list(shape = 22, size = 2.6, colour = "black", stroke = 0.3, alpha = 1,
                                   fill = unname(paper_cols$species[c("Human", "Yeast")])))
}

# Height of the default geom_density() curve at `at`, and its peak
density_at <- function(x, at){
  d <- stats::density(x[is.finite(x)])
  stats::approx(d$x, d$y, xout = at)$y
}
density_peak <- function(x) max(stats::density(x[is.finite(x)])$y)

# The three ratio panels of one comparison: scatter vs log2 denominator,
# box plots and densities per species
quant_panels <- function(dat, numerator = "A", denominator = "B", exp_FC = 2, alphaScaler = 1,
                         probs = c(0, 1), x_lim = NULL, y_lim = NULL, dens_lim = NULL,
                         x_axis = TRUE, legend = FALSE, n_label = TRUE, point_size = 1.2){
  dat <- ratio_data(dat, numerator, denominator)
  if (is.null(y_lim)) y_lim <- as.numeric(quantile(sort(dat$FCs_corrected), probs = probs))
  if (is.null(x_lim)) x_lim <- as.numeric(quantile(sort(dat$denominator[is.finite(dat$denominator)]), probs = probs))
  if (is.null(dens_lim)) dens_lim <- c(0, 1.05 * max(tapply(dat$FCs_corrected, dat$species, density_peak)))
  dat$Alpha <- ifelse(grepl("Human", dat$species), 0.04, 0.16) * alphaScaler
  sp <- paper_cols$species[c("Human", "Yeast")]
  drop_x <- if (x_axis) NULL else paper_drop_x()

  scatter <- ggplot(dat, aes(x = denominator, y = FCs_corrected, colour = species)) +
    expected_lines(exp_FC) +
    paper_raster(geom_point(aes(alpha = Alpha), size = point_size, stroke = 0)) +
    scale_colour_manual(values = sp, name = NULL) +
    scale_alpha_identity() +
    scale_x_continuous(breaks = paper_int_breaks(2), labels = paper_int_labels) +
    scale_y_continuous(breaks = paper_int_breaks(2)) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(x = bquote(Log[2]*", "*.(denominator)),
         y = bquote(Log[2]*", "*.(paste0(numerator, "/", denominator))),
         subtitle = if (n_label) bquote(italic(n)*" = "*.(scales::comma(nrow(dat)))*" precursor ratios")) +
    guides(colour = species_guide()) +
    paper_theme() + drop_x +
    theme(legend.position = if (legend) "inside" else "none",
          legend.position.inside = c(1, 1), legend.justification = c(1, 1),
          plot.margin = margin(1, 0, 1, 0))

  box <- ggplot(dat, aes(x = species, y = FCs_corrected, fill = species)) +
    expected_lines(exp_FC) +
    geom_boxplot(alpha = 0.85, outlier.shape = NA, linewidth = 0.25) +
    paper_box_outliers(alpha = 0.01 * alphaScaler) +
    scale_fill_manual(values = sp) +
    scale_y_continuous(breaks = paper_int_breaks(2)) +
    coord_cartesian(ylim = y_lim) +
    paper_theme() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks.x = element_blank(), plot.margin = margin(1, 1, 1, 1))

  dens_brk <- c(0, signif(0.8 * dens_lim[2], 1))
  dense <- ggplot(dat, aes(y = FCs_corrected, fill = species)) +
    geom_density(alpha = 0.75, linewidth = 0.25) +
    expected_lines(exp_FC) +
    scale_fill_manual(values = sp) +
    scale_x_continuous(breaks = dens_brk, labels = scales::label_number(accuracy = 0.1)) +
    scale_y_continuous(breaks = paper_int_breaks(2)) +
    coord_cartesian(xlim = dens_lim, ylim = y_lim, expand = FALSE) +
    labs(x = "density") +
    paper_theme() + drop_x +
    theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(),
          plot.margin = margin(1, 1, 1, 1))

  list(scatter = scatter, box = box, dense = dense, n = nrow(dat))
}

# Scatter | box | density for one comparison, as one plot (per-run checks)
quant <- function(dat, numerator = "A", denominator = "B", exp_FC = 2,
                  title_actual = "3-plexDIA & 3-timePlex: 8 ng HY (uPAC)", alphaScaler = 1,
                  probs = c(0, 1), x_lim=NULL, y_lim=NULL){
  p <- quant_panels(dat, numerator, denominator, exp_FC, alphaScaler, probs, x_lim, y_lim, legend = TRUE)
  patchwork::wrap_plots(p$scatter, p$box, p$dense, nrow = 1, widths = c(0.8, 0.2, 0.33)) +
    patchwork::plot_annotation(title = title_actual,
                               theme = theme(plot.title = element_text(size = paper_size("title"),
                                                                       face = "bold", hjust = 0.5)))
}

# Same filtering + human-median correction as quant(), returns the y quantile range
fc_limits <- function(dat, numerator, denominator, probs = c(.005, .995)){
  quant_col <- paste0(numerator, denominator)
  FCs <- dat[[quant_col]]
  keep <- !is.na(FCs) & is.finite(FCs) & FCs != 0
  FCs <- FCs[keep]
  H_med <- median(FCs[grepl("Human", dat$species[keep])])
  as.numeric(quantile(sort(FCs - H_med), probs = probs))
}

# Denominator (x-axis) quantile range, using the same row filter as quant()
den_limits <- function(dat, numerator, denominator, probs = c(.005, .995)){
  FCs <- dat[[paste0(numerator, denominator)]]
  keep <- !is.na(FCs) & is.finite(FCs) & FCs != 0
  den <- dat[[denominator]][keep]
  as.numeric(quantile(sort(den[is.finite(den)]), probs = probs))
}

# Precursor intensity relative to its channel median vs RT; with trend = TRUE
# the RT normalization factor is drawn in red
rt_intensity_panel <- function(d, y, trend = FALSE, x_axis = TRUE, ylab = "Log<sub>10</sub>, relative<br>intensity"){
  p <- ggplot(d, aes(x = rt, y = .data[[y]])) +
    paper_raster(ggpointdensity::geom_pointdensity(alpha = 0.08, size = 1, stroke = 0)) +
    scale_colour_viridis_c(option = "mako") +
    scale_y_continuous(limits = c(-3, 3), breaks = c(-2, 0, 2)) +
    labs(x = "RT (min)", y = ylab) +
    paper_theme() +
    theme(legend.position = "none", axis.title.y = ggtext::element_markdown(size = paper_size("text"), lineheight = 1.1),
          plot.margin = margin(1, 0, 1, 0))
  if (trend) p <- p + paper_raster(geom_point(aes(y = log10(1 / rtnorm_factor_MS1)), colour = "red", size = 0.4, stroke = 0))
  if (!x_axis) p <- p + paper_drop_x()
  p
}


## ---- Single run before vs after RT correction (Fig. 4a) ---------------------

# Fig. 4a: one channel's intensities vs RT and one comparison's ratios, before
# (top row) and after (bottom row) the RT correction. The rows share x axes,
# so only the bottom row carries x tick labels.
fig4a_before_after <- function(chan_dat, before_d, after_d, numerator, denominator, exp_FC,
                               x_lim, y_lim, alphascaler = 1){
  dens_lim <- c(0, 1.05 * max(sapply(list(before_d, after_d), function(d){
    r <- ratio_data(d, numerator, denominator)
    max(tapply(r$FCs_corrected, r$species, density_peak))
  })))
  qb <- quant_panels(before_d, numerator, denominator, exp_FC, alphascaler, x_lim = x_lim, y_lim = y_lim,
                     dens_lim = dens_lim, x_axis = FALSE, legend = TRUE, n_label = TRUE)
  qa <- quant_panels(after_d, numerator, denominator, exp_FC, alphascaler, x_lim = x_lim, y_lim = y_lim,
                     dens_lim = dens_lim, x_axis = TRUE, legend = FALSE, n_label = FALSE)
  smp <- unique(chan_dat$Sample); stopifnot(length(smp) == 1)
  smp <- ifelse(smp %in% c("A", "B", "C"), smp, ifelse(grepl("C1|D1", smp), "A", ifelse(grepl("C2|D2", smp), "B", "C")))
  ylab <- paste0("Log<sub>10</sub>, rel. intensity<br>(sample ", smp, ")")   # two lines: one is taller than the panel
  row_label <- function(txt){
    ggplot() + annotate("text", x = 0, y = 0, label = txt, angle = 90, fontface = "bold",
                        size = paper_geom_text(12)) + theme_void()
  }
  patchwork::wrap_plots(row_label("Before"), rt_intensity_panel(chan_dat, "plex_Area_norm", trend = TRUE, x_axis = FALSE, ylab = ylab),
                        qb$scatter, qb$box, qb$dense,
                        row_label("After"), rt_intensity_panel(chan_dat, "plex_Area_rtnorm_norm", ylab = ylab),
                        qa$scatter, qa$box, qa$dense,
                        ncol = 5, widths = grid::unit(c(0.22, 1.47, 2.13, 0.51, 0.49), c("in", rep("null", 4)))) +
    patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
}

# Per run: intensity vs RT for every channel and the within-run ratio plots
# (folders named by run). The run/channel/comparison chosen for Fig. 4a is
# also saved to out_dir as a manuscript panel.
single_run_quant <- function(normed_quant, use_best_channel = F,quant_col="coeff",alphascaler=1,
                             fig_run = "JD0422", fig_chan = "JD0422_2_LF", fig_comparison = "BC",
                             out_dir = "Figs/Fig4"){

  normed_quant <- normed_quant[normed_quant$Repeatability==F,] #remove Runs that are just for computing CVs
  normed_quant$quant_val_beforeAdj <- normed_quant[[quant_col]]
  normed_quant <- normed_quant %>% dplyr::group_by(run_chan) %>%
    dplyr::mutate(plex_Area_rtnorm_norm = log10(plex_Area_rtnorm) - log10(median(plex_Area_rtnorm,na.rm=T)),
                  plex_Area_norm= log10(plex_Area) - log10(median(plex_Area,na.rm=T))) %>% ungroup()

  normed_quant <- normed_quant %>% dplyr::mutate("Sample" = ifelse(grepl("C1|D1",Sample),"A",
                                                                   ifelse(grepl("C2|D2",Sample),"B","C")))

  normed_quant1 <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",ProteinName),"Human",
                                                                     ifelse(grepl("YEAST",ProteinName),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))
  normed_quant1_noplex <- normed_quant1[(normed_quant1$timePlex==F) & (normed_quant1$Label=="LF"),]
  normed_quant1_plex <- normed_quant1[(normed_quant1$timePlex==T)|(normed_quant1$Label!="LF"),]

  comparisons <- list(AB = list(num = "A", den = "B", exp_FC = -2),
                      AC = list(num = "A", den = "C", exp_FC = 1),
                      BC = list(num = "B", den = "C", exp_FC = 3))
  save_quant <- function(p, file) paper_save(p, file, width = 4.6, height = 1.8)
  fig4a <- NULL

  uni_runs <- unique(normed_quant1_plex$Run)
  for(i in 1:length(uni_runs)){
    subDir <- uni_runs[i]
    dir.create(file.path(subDir))
    run_of_interest <- normed_quant1[normed_quant1$Run==uni_runs[i],]

    uni_runchans <- unique(run_of_interest$run_chan)

    for(j in 1:length(unique(uni_runchans))){

      dat_singleSample <- run_of_interest[run_of_interest$run_chan==uni_runchans[j],]
      paper_save(rt_intensity_panel(dat_singleSample, "plex_Area_rtnorm_norm") +
                   labs(title = paste0("After adjustment (", uni_runchans[j], ")")),
                 paste0(subDir,"/AFTER_RT_adj_",uni_runchans[j],".pdf"), width = 2.6, height = 2)
      paper_save(rt_intensity_panel(dat_singleSample, "plex_Area_norm", trend = TRUE) +
                   labs(title = paste0("Before adjustment (", uni_runchans[j], ")")),
                 paste0(subDir,"/BEFORE_RT_adj_",uni_runchans[j],".pdf"), width = 2.6, height = 2)
    }
    uni_current <- unique(run_of_interest[1,])
    if(uni_current$timePlex==T){

      if(!grepl("LF",uni_current$Label)){ #plexDIA + timePlex

        for(k in 0:(length(unique(run_of_interest$time_channel))-1)){

          within_pD <- run_of_interest[run_of_interest$time_channel==k,]
          within_pD <- within_pD %>% dplyr::distinct(seqcharge, run_chan,.keep_all=T)
          within_pD_d <- dcast(within_pD, seqcharge+ProteinName+species~Sample, value.var = "norm")
          within_pD_d$AB <- within_pD_d$A-within_pD_d$B
          within_pD_d$BC <- within_pD_d$B-within_pD_d$C
          within_pD_d$AC <- within_pD_d$A-within_pD_d$C

          for (cmp in names(comparisons)){
            cc <- comparisons[[cmp]]
            save_quant(quant(within_pD_d, numerator = cc$num, denominator = cc$den, exp_FC = cc$exp_FC,
                             title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA, timePlex-",k),
                             alphaScaler = alphascaler),
                       paste0(subDir,"/TP_PD_Quant_",uni_runs[i],"_",cmp,"_",k,"_timeplex_bestchannel_",use_best_channel,"_",quant_col,".pdf"))
          }

        }

      } else if((uni_current$timePlex==F) & (!grepl("LF",uni_current$Label))){

        within_pD <- run_of_interest
        within_pD <- within_pD %>% dplyr::distinct(seqcharge, run_chan,.keep_all=T)
        within_pD_d <- dcast(within_pD, seqcharge+ProteinName+species~Sample, value.var = "norm")
        within_pD_d$AB <- within_pD_d$A-within_pD_d$B
        within_pD_d$BC <- within_pD_d$B-within_pD_d$C
        within_pD_d$AC <- within_pD_d$A-within_pD_d$C

        for (cmp in names(comparisons)){
          cc <- comparisons[[cmp]]
          save_quant(quant(within_pD_d, numerator = cc$num, denominator = cc$den, exp_FC = cc$exp_FC,
                           title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA"),
                           alphaScaler = alphascaler),
                     paste0(subDir,"/PD_Quant_",uni_runs[i],"_",cmp,"_bestchannel_",use_best_channel,"_",quant_col,".pdf"))
        }

      } else{


        within_run <- run_of_interest %>% dplyr::distinct(seqcharge, run_chan, .keep_all = T)

        #### after RT adjustment
        within_run$log2quant_val <- log2(within_run$quant_val)
        within_run_d <- dcast(within_run, seqcharge + ProteinName + species ~ Sample, value.var = "log2quant_val")
        within_run_d$AB <- within_run_d$A - within_run_d$B
        within_run_d$BC <- within_run_d$B - within_run_d$C
        within_run_d$AC <- within_run_d$A - within_run_d$C

        #### before RT adjustment
        within_run$log2quant_val_beforeAdj <- log2(within_run$quant_val_beforeAdj)
        within_run_d_before <- dcast(within_run, seqcharge + ProteinName + species ~ Sample, value.var = "log2quant_val_beforeAdj")
        within_run_d_before$AB <- within_run_d_before$A - within_run_d_before$B
        within_run_d_before$BC <- within_run_d_before$B - within_run_d_before$C
        within_run_d_before$AC <- within_run_d_before$A - within_run_d_before$C

        #### shared limits per comparison: min of both lower bounds, max of both upper bounds
        lim_probs <- c(.005, .995)
        shared_lim <- function(f, num, den){
          r <- rbind(f(within_run_d,        num, den, probs = lim_probs),
                     f(within_run_d_before, num, den, probs = lim_probs))
          c(min(r[, 1]), max(r[, 2]))
        }

        title_wr <- "3-timePlex: 8 ng HY... within run"
        fbase <- paste0(subDir, "/TP_Quant_", uni_runs[i], "_")

        for (cmp in names(comparisons)){
          cc <- comparisons[[cmp]]
          ylim_c <- shared_lim(fc_limits, cc$num, cc$den); xlim_c <- shared_lim(den_limits, cc$num, cc$den)

          #### after and before RT adjustment
          save_quant(quant(within_run_d, numerator = cc$num, denominator = cc$den, exp_FC = cc$exp_FC,
                           title_actual = title_wr, alphaScaler = alphascaler, x_lim = xlim_c, y_lim = ylim_c),
                     paste0(fbase, cmp, "_bestchannel_", use_best_channel, "_", quant_col, ".pdf"))
          save_quant(quant(within_run_d_before, numerator = cc$num, denominator = cc$den, exp_FC = cc$exp_FC,
                           title_actual = title_wr, alphaScaler = alphascaler, x_lim = xlim_c, y_lim = ylim_c),
                     paste0(fbase, cmp, "_bestchannel_", use_best_channel, "_beforeAdj", quant_col, ".pdf"))

          #### Fig. 4a
          if (uni_runs[i] == fig_run && cmp == fig_comparison){
            p <- fig4a_before_after(run_of_interest[run_of_interest$run_chan == fig_chan,],
                                    within_run_d_before, within_run_d, cc$num, cc$den, cc$exp_FC,
                                    x_lim = xlim_c, y_lim = ylim_c, alphascaler = alphascaler)
            fig4a <- paper_panel(p, file.path(out_dir, "Fig4a.pdf"), width = 5.45, height = 2.45)
          }
        }
      }

    }}



  for(p in 1:(length(unique(normed_quant1_noplex$Rep)))){
    subDir <- paste0("NoPlex_LF_Rep_",p)
    dir.create(file.path(subDir))

    reps_noplex <- normed_quant1_noplex %>% dplyr::distinct(Run,.keep_all=T) %>%
      dplyr::group_by(Rep,Sample) %>%
      dplyr::mutate(Rep = n()) %>%
      dplyr::ungroup()
    within_Rep <- normed_quant1_noplex[normed_quant1_noplex$Rep==p,]
    if(length(unique(within_Rep$Sample))==3){

      within_Rep <- dcast(within_Rep, seqcharge+ProteinName+species~Sample, value.var = "norm")
      within_Rep$AB <- within_Rep$A-within_Rep$B
      within_Rep$BC <- within_Rep$B-within_Rep$C
      within_Rep$AC <- within_Rep$A-within_Rep$C

      for (cmp in names(comparisons)){
        cc <- comparisons[[cmp]]
        save_quant(quant(within_Rep, numerator = cc$num, denominator = cc$den, exp_FC = cc$exp_FC,
                         title_actual = paste0("LF NoPlex: 8 ng HY... Rep ", p), alphaScaler = alphascaler),
                   paste0(subDir,"/LF_NoPlex_Quant_Rep",uni_runs[p],"_",cmp,"_bestchannel_",use_best_channel,"_",quant_col,".pdf"))
      }

    }

  }
  if (is.null(fig4a)) warning("Fig. 4a not made: no LF timePlex run '", fig_run, "' with comparison ", fig_comparison)
  invisible(list(fig4a = fig4a))
}


## ---- Ratio densities per species, faceted (Fig. 4b, Fig. 5d) ----------------

# Densities of log2 ratios per species with the expected ratios (dashed) and
# each distribution's median (solid line drawn up to the density curve)
ratio_density_plot <- function(d, facet, exp_FC, ylab, legend = FALSE, x_breaks = waiver()){
  med <- d %>% dplyr::group_by(.data[[facet]], species) %>%
    dplyr::summarise(median_value = median(value), dens = density_at(value, median(value)),
                     .groups = "drop")
  ggplot(d, aes(y = value, fill = species)) +
    geom_density(alpha = 0.75, linewidth = 0.25) +
    expected_lines(exp_FC) +
    geom_segment(data = med, aes(x = 0, xend = dens, y = median_value, yend = median_value),
                 inherit.aes = FALSE, linewidth = 0.3) +
    facet_grid(cols = vars(.data[[facet]])) +
    scale_fill_manual(values = paper_cols$species[c("Human", "Yeast")], name = NULL) +
    scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0, 0.05))) +
    labs(x = "density", y = ylab) +
    guides(fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.25))) +
    paper_theme_box() +
    theme(legend.position = if (legend) "inside" else "none",
          legend.position.inside = c(0.95, 0.95), legend.justification = c(1, 1))
}


## ---- Separate vs single emitter implementations (Fig. 4b) -------------------

quant_multi_single_emitter <- function(IO, uPAC, grouper="seqcharge", datatype = "precursor", min_num_reps=2, alphascaler=1,
                                       out_dir = "Figs/Fig4"){

  IO$implment <- "Multi"
  uPAC$implment <- "Single"

  normed_quant <- rbind.fill(IO, uPAC)
  normed_quant <- normed_quant %>% dplyr::mutate("Sample" = ifelse(Sample %in% c("A","B","C"), Sample,
                                                                   ifelse(grepl("C1|D1",Sample),"A",
                                                                          ifelse(grepl("C2|D2",Sample),"B","C"))))

  normed_quant <- normed_quant %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tp_pd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))

  normed_quant <- normed_quant[normed_quant$Repeatability==F,]
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",ProteinName),"Human",
                                                                    ifelse(grepl("YEAST",ProteinName),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))

  tp <- normed_quant[normed_quant$Type=="tp",] %>% dplyr::filter(Rep == 2) #get the middle set of triplicates

  joined <- tp

  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()

  joined <- joined %>%
    left_join(human_medians, by = "run_chan") %>%
    dplyr::mutate(norm = norm-median_Human) %>%
    dplyr::select(-median_Human)


  joined_quant <- joined %>% dplyr::group_by(Type, Sample,.data[[grouper]],species,implment) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_num_reps-1)) %>% #require atleast 2/3
    dplyr::summarize(med_quant = median(norm,na.rm=T))

  dynamic_column <- grouper

  formula_str <- paste(dynamic_column, "species", sep = "+")
  joined_quant_d <- dcast(
    joined_quant,
    formula = as.formula(paste(formula_str, "~ Type+Sample+implment")),
    value.var = "med_quant"
  )


  joined_quant_d$tp_BC_single <- joined_quant_d$tp_B_Single-joined_quant_d$tp_C_Single
  joined_quant_d$tp_BC_multi <- joined_quant_d$tp_B_Multi-joined_quant_d$tp_C_Multi

  joined_quant_d$tp_AC_single <- joined_quant_d$tp_A_Single-joined_quant_d$tp_C_Single
  joined_quant_d$tp_AC_multi <- joined_quant_d$tp_A_Multi-joined_quant_d$tp_C_Multi

  joined_quant_d$tp_AB_single <- joined_quant_d$tp_A_Single-joined_quant_d$tp_B_Single
  joined_quant_d$tp_AB_multi <- joined_quant_d$tp_A_Multi-joined_quant_d$tp_B_Multi


  joined_quant_m <- joined_quant_d %>% dplyr::select(.data[[grouper]],species,
                                                     tp_BC_single, tp_BC_multi,
                                                     tp_AC_single, tp_AC_multi,
                                                     tp_AB_single, tp_AB_multi
  ) %>% melt()


  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()

  joined_quant_m <- joined_quant_m %>%
    left_join(human_medians, by = "variable") %>%
    dplyr::mutate(value = value-median_Human) %>%
    dplyr::select(-median_Human) %>% na.omit()

  ## one density plot per comparison, precursors intersected across the two
  ## implementations; B/C is the manuscript panel, A/B and A/C are checks
  emitter_labs <- c(multi = "Separate\nemitters", single = "Single\nemitter")
  comparisons <- list(BC = list(num = "B", den = "C", exp_FC = 3),
                      AB = list(num = "A", den = "B", exp_FC = -2),
                      AC = list(num = "A", den = "C", exp_FC = 1))
  fig4b <- NULL
  for (cmp in names(comparisons)){
    cc <- comparisons[[cmp]]
    joined_cmp <- joined_quant_m[grepl(cmp, joined_quant_m$variable),]
    int_cmp <- joined_cmp %>% na.omit() %>% dplyr::count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
    joined_cmp_int <- joined_cmp %>% dplyr::filter(.data[[grouper]]%in%int_cmp[[grouper]]) %>% na.omit()
    joined_cmp_int$emitter <- factor(emitter_labs[sub(".*_", "", joined_cmp_int$variable)], levels = emitter_labs)

    p <- ratio_density_plot(joined_cmp_int, facet = "emitter", exp_FC = cc$exp_FC, x_breaks = c(0, 0.5),
                            ylab = bquote(Log[2]*", "*.(paste0(cc$num, "/", cc$den)))) +
      scale_y_continuous(breaks = paper_int_breaks(2)) +
      labs(title = "LF 3-timePlex",
           subtitle = bquote(italic(n)*" = "*.(scales::comma(nrow(int_cmp)))*" precursors"))

    if (cmp == "BC"){
      fig4b <- paper_panel(p, file.path(out_dir, "Fig4b.pdf"), width = 1.64, height = 2.45)
    } else {
      paper_save(p, file.path(out_dir, paste0("Fig4b_", cmp, ".pdf")), width = 1.64, height = 2.45)   # other comparisons
    }
  }

  invisible(list(fig4b = fig4b))
}


## ---- LF no-plex, 3-plexDIA, 3-timePlex and 3x3 (Fig. 4c-e) ------------------

quant_all_comparisons_v2 <- function(normed_quant, grouper="seqcharge", extra="", datatype = "precursor", alphascaler=1, min_num_reps=2,
                                     out_dir = "Figs/Fig4"){

  normed_quant <- normed_quant %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tp_pd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))

  normed_quant <- normed_quant[normed_quant$Repeatability==F,] #remove Runs that are just for computing CVs
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
  dynamic_column <- grouper

  formula_str <- paste(dynamic_column, "species", sep = "+")
  joined_quant_d <- dcast(
    joined_quant,
    formula = as.formula(paste(formula_str, "~ Type+Sample")),
    value.var = "med_quant"
  )
  joined_quant_d$np_AB <- joined_quant_d$np_A-joined_quant_d$np_B
  joined_quant_d$np_AC <- joined_quant_d$np_A-joined_quant_d$np_C
  joined_quant_d$np_BC <- joined_quant_d$np_B-joined_quant_d$np_C

  joined_quant_d$tp_AB <- joined_quant_d$tp_A-joined_quant_d$tp_B
  joined_quant_d$tp_AC <- joined_quant_d$tp_A-joined_quant_d$tp_C
  joined_quant_d$tp_BC <- joined_quant_d$tp_B-joined_quant_d$tp_C

  joined_quant_d$pd_AB <- joined_quant_d$pd_A-joined_quant_d$pd_B
  joined_quant_d$pd_AC <- joined_quant_d$pd_A-joined_quant_d$pd_C
  joined_quant_d$pd_BC <- joined_quant_d$pd_B-joined_quant_d$pd_C

  joined_quant_d$tp_pd_AB <- joined_quant_d$tp_pd_A-joined_quant_d$tp_pd_B
  joined_quant_d$tp_pd_AC <- joined_quant_d$tp_pd_A-joined_quant_d$tp_pd_C
  joined_quant_d$tp_pd_BC <- joined_quant_d$tp_pd_B-joined_quant_d$tp_pd_C

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

  allcounts  <- joined_quant_m %>% na.omit() %>% dplyr::count(variable)

  ## ---- one row per comparison: ratio counts | scatter per data type | box plots
  ## c = A/B, d = A/C, e = B/C. Ratios are intersected across the four data types.
  comparisons <- list(AB = list(num = "A", den = "B", exp_FC = -2, panel = "c"),
                      AC = list(num = "A", den = "C", exp_FC = 1,  panel = "d"),
                      BC = list(num = "B", den = "C", exp_FC = 3,  panel = "e"))
  type_labs  <- c(np = "LF no-plex", pd = "3-plexDIA", tp = "3-timePlex", tp_pd = "3x3 plexDIA & timePlex")
  short_labs <- c(np = "LF np", pd = "3-pD", tp = "3-tP", tp_pd = "3x3", Intersected = "Intersect")
  level_lab  <- ifelse(grouper == "seqcharge", "Precursor", "Protein")
  sp <- paper_cols$species[c("Human", "Yeast")]

  rows <- lapply(names(comparisons), function(cmp){
    joined_cmp <- joined_quant_m[grepl(cmp, joined_quant_m$variable),]
    counts_cmp <- allcounts[grepl(cmp, allcounts$variable),]
    int_cmp <- joined_cmp %>% na.omit() %>% dplyr::count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
    counts_cmp <- rbind(counts_cmp, data.frame(variable = "Intersected", n = nrow(int_cmp)))
    counts_cmp$type <- factor(sub(paste0("_", cmp, "$"), "", counts_cmp$variable), levels = names(short_labs))
    joined_cmp_int <- joined_cmp %>% dplyr::filter(.data[[grouper]]%in%int_cmp[[grouper]]) %>% na.omit()
    joined_cmp_int$type <- factor(sub(paste0("_", cmp, "$"), "", joined_cmp_int$variable), levels = names(type_labs))
    joined_cmp_int$Alpha <- ifelse(grepl("Human", joined_cmp_int$species), 0.04, 0.12) * alphascaler
    list(counts = counts_cmp, int = joined_cmp_int,
         x_lim = as.numeric(quantile(sort(log2(joined_cmp_int$quant_avg)), probs = c(.005, .995))),
         y_lim = as.numeric(quantile(sort(joined_cmp_int$value), probs = c(.005, .995))))
  })
  names(rows) <- names(comparisons)

  ## the same x range and ticks for every scatter plot in c-e; the rows share
  ## their x axes, so tick labels and titles appear only under e and the facet
  ## strips only above c
  x_lim <- range(sapply(rows, `[[`, "x_lim"))
  
  ce_grid <- function(letters = FALSE){
    plots <- list()
    for (k in seq_along(comparisons)){
      cmp <- names(comparisons)[k]; cc <- comparisons[[cmp]]; r <- rows[[cmp]]
      first <- k == 1; last <- k == length(comparisons)
      ratio_lab <- paste0(cc$num, "/", cc$den)
      
      counts_plot <- ggplot(r$counts, aes(x = type, y = n, fill = type)) +
        geom_col(colour = "black", alpha = 0.7, linewidth = 0.25, width = 0.9) +
        scale_fill_manual(values = c("grey10", "grey30", "grey50", "grey70", "grey90")) +
        scale_x_discrete(labels = short_labs) +
        scale_y_continuous(labels = scales::comma, breaks = scales::breaks_pretty(n = 3),
                           expand = expansion(mult = c(0, 0.14))) +
        labs(x = NULL, y = paste0("# ratios (", ratio_lab, ")")) +
        paper_theme() +
        theme(legend.position = "none", axis.ticks.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              plot.margin = margin(1, 0, 1, 0))
      
      scatter <- ggplot(r$int, aes(x = log2(quant_avg), y = value, colour = species, alpha = Alpha)) +
        expected_lines(cc$exp_FC) +
        paper_raster(geom_point(size = 0.8, stroke = 0)) +
        facet_grid(cols = vars(type), labeller = as_labeller(type_labs)) +
        scale_colour_manual(values = sp) +
        scale_alpha_identity() +
        scale_x_continuous(breaks = paper_int_breaks(2), labels = paper_int_labels) +
        scale_y_continuous(breaks = paper_int_breaks(1)) +
        coord_cartesian(xlim = x_lim, ylim = r$y_lim) +
        labs(x = bquote(Log[2]*", denominator"), y = bquote(Log[2]*", "*.(ratio_lab))) +
        paper_theme_box() +
        theme(legend.position = "none", strip.text = element_text(size = paper_size("text")),
              panel.spacing.x = unit(2, "pt"), plot.margin = margin(1, 0, 1, 0))
      
      box <- ggplot(r$int, aes(x = type, y = value, fill = species)) +
        expected_lines(cc$exp_FC) +
        geom_boxplot(alpha = 0.85, outlier.shape = NA, linewidth = 0.25) +
        paper_box_outliers(alpha = alphascaler * 0.01) +
        scale_fill_manual(values = sp) +
        scale_x_discrete(labels = short_labs) +
        scale_y_continuous(breaks = paper_int_breaks(1)) +
        coord_cartesian(ylim = r$y_lim) +
        paper_theme_box() +
        theme(axis.title = element_blank(), axis.text.y = element_blank(), legend.position = "none",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              plot.margin = margin(1, 1, 1, 2))
      
      if (!first) scatter <- scatter + theme(strip.text = element_blank(), strip.background = element_blank())
      if (!last) {
        counts_plot <- counts_plot + paper_drop_x()
        scatter <- scatter + paper_drop_x()
        box <- box + paper_drop_x()
      } else {
        # rotated tick labels of the bar and box plots hang into the bottom
        # margin instead of pushing the scatter plots' x-axis title down
        counts_plot <- patchwork::free(counts_plot, type = "space", side = "b")
        box <- patchwork::free(box, type = "space", side = "b")
      }
      if (letters) {
        tag <- ggplot() +
          annotate("text", x = 0, y = 1, label = cc$panel, hjust = 0, vjust = if (first) 0.1 else 1,
                   fontface = "bold", size = paper_geom_text(paper_size("tag"))) +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
          theme_void()
        plots <- c(plots, list(tag))
      }
      plots <- c(plots, list(counts_plot, scatter, box))
    }
    widths <- if (letters) grid::unit(c(0.16, 0.43, 3.62, 0.62), c("in", "null", "null", "null"))
              else c(0.43, 3.62, 0.62)
    patchwork::wrap_plots(plots, ncol = length(widths), widths = widths) +
      patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 1, 14, 0)))
  }
  
  # the individual PDF has no panel letters; the assembled figure's copy has c, d, e
  ce_file <- file.path(out_dir, "Fig4c-e.pdf")
  paper_save(ce_grid(), ce_file, width = 5.44, height = 2.68)
  fig4ce <- list(grob = paper_grob(ce_grid(letters = TRUE)), width = 5.6, height = 2.68, file = ce_file)
  
  invisible(list(fig4ce = fig4ce))
}


## ---- MS1 vs MS2 accuracy (Fig. 4f) ------------------------------------------

compare_quants <- function(dat1=normed_quant_IO_MS1_bestfit_relaxed_SR1,
                           dat2=normed_quant_IO_relaxed_coeff,
                           quant1="MS1_bestfit",
                           quant2="coeff", min_num_reps=2, grouper="seqcharge",
                           quant_labels = c("MS1 quantitation", "MS2 quantitation"),
                           outfile="Figs/Fig4/Fig4f.pdf"){

  dat1$quant <- "quant1"
  dat2$quant <- "quant2"
  normed_quant <- rbind.fill(dat1,dat2)

  normed_quant <- normed_quant %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tppd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))

  if (!"ProteinName" %in% colnames(normed_quant) || all(is.na(normed_quant$ProteinName))) {
    normed_quant$ProteinName <- normed_quant$prot
  }
  normed_quant <- normed_quant[normed_quant$Repeatability==F,] #remove Runs that are just for computing CVs
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",ProteinName),"Human",
                                                                    ifelse(grepl("YEAST",ProteinName),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))

  normed_quant <- normed_quant %>% dplyr::mutate("Sample" = ifelse(grepl("C1|D1",Sample),"A",
                                                                   ifelse(grepl("C2|D2",Sample),"B","C")))

  np <- normed_quant[normed_quant$Type=="np",]

  #take middle replicate for all so that number of samples matches number of samples for LF no-plex
  pD <- normed_quant[normed_quant$Type=="pd",] %>% dplyr::filter(grepl("JD0425|JD0426|JD0427",run_chan))
  tp_pd <- normed_quant[normed_quant$Type=="tppd",] %>% dplyr::filter(Rep == 2)
  tp <- normed_quant[normed_quant$Type=="tp",] %>% dplyr::filter(Rep == 2)

  joined <- rbind(np, pD, tp_pd, tp)

  # per run/channel (and per quant) human-median centering, as in quant_all_comparisons_v2
  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan, quant) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop")

  joined <- joined %>%
    left_join(human_medians, by = c("run_chan", "quant")) %>%
    dplyr::mutate(norm = norm-median_Human) %>%
    dplyr::select(-median_Human)

  joined_quant <- joined %>% dplyr::group_by(Type, Sample,.data[[grouper]],species,quant) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_num_reps-1)) %>% #require atleast 1/3
    dplyr::summarize(med_quant = median(norm,na.rm=T), .groups = "drop")

  formula_str <- paste(grouper, "species", sep = "+")
  joined_quant_d <- dcast(
    joined_quant,
    formula = as.formula(paste(formula_str, "~ Type+Sample+quant")),
    value.var = "med_quant"
  )

  # same 12 differences as quant_all_comparisons_v2, once per quant
  prefixes    <- c("np", "pd", "tp", "tppd")
  comparisons <- c("AB", "AC", "BC")
  suffixes    <- c("quant1", "quant2")
  comparison_cols <- c()
  for (prefix in prefixes) for (suffix in suffixes) for (cmp in comparisons) {
    col1 <- paste0(prefix, "_", substr(cmp, 1, 1), "_", suffix)
    col2 <- paste0(prefix, "_", substr(cmp, 2, 2), "_", suffix)
    if (!all(c(col1, col2) %in% colnames(joined_quant_d))) {
      warning("missing dcast column(s) for ", prefix, " ", cmp, " ", suffix, " - skipped")
      next
    }
    colname <- paste0(prefix, "_", cmp, "_", suffix)
    joined_quant_d[[colname]] <- joined_quant_d[[col1]] - joined_quant_d[[col2]]
    comparison_cols <- c(comparison_cols, colname)
  }
  if (length(comparison_cols) == 0) stop("No Type/Sample/quant combinations had data. dcast columns: ",
                                         paste(colnames(joined_quant_d), collapse = ", "))

  joined_quant_m <- joined_quant_d %>%
    dplyr::select(all_of(c(grouper, "species", comparison_cols))) %>%
    melt(id.vars = c(grouper, "species"))

  # center each ratio column on its human median, as in quant_all_comparisons_v2
  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop")

  joined_quant_m <- joined_quant_m %>%
    left_join(human_medians, by = "variable") %>%
    dplyr::mutate(value = value-median_Human) %>%
    dplyr::select(-median_Human) %>% na.omit()

  joined_quant_m$variable   <- as.character(joined_quant_m$variable)
  joined_quant_m$data_type  <- sub("_.*", "", joined_quant_m$variable)
  joined_quant_m$comparison <- sub("^[^_]+_([A-C]{2})_.*$", "\\1", joined_quant_m$variable)
  joined_quant_m$quant_type <- sub(".*_", "", joined_quant_m$variable)
  joined_quant_m <- joined_quant_m %>% dplyr::mutate(theoretical = dplyr::case_when(
    species == "Yeast" & comparison == "AB" ~ -2,
    species == "Yeast" & comparison == "BC" ~  3,
    species == "Yeast" & comparison == "AC" ~  1,
    TRUE ~ 0))

  # ---- intersect PER COMPARISON (AB, BC, AC separately), as quant_all_comparisons_v2 does ----
  # a precursor is kept for a comparison if it has that ratio in every data type x quant present
  int_list <- list()
  count_rows <- list()
  for (cmp in comparisons) {
    joined_cmp <- joined_quant_m[joined_quant_m$comparison == cmp,] %>% na.omit()
    if (nrow(joined_cmp) == 0) next
    allcounts <- joined_cmp %>% dplyr::count(variable)
    int_cmp <- joined_cmp %>% dplyr::count(.data[[grouper]]) %>% dplyr::filter(n == max(n))
    count_rows[[cmp]] <- rbind(allcounts, data.frame(variable = paste0("Intersected_", cmp), n = nrow(int_cmp)))
    int_list[[cmp]] <- joined_cmp %>% dplyr::filter(.data[[grouper]] %in% int_cmp[[grouper]])
  }
  counts_tbl <- do.call(rbind, count_rows)
  message("compare_quants: yeast precursor ratios per column and intersected per comparison:")
  print(counts_tbl, row.names = FALSE)

  joined_quant_m_int <- do.call(rbind, int_list)
  joined_quant_m_int$residual <- abs(joined_quant_m_int$value - joined_quant_m_int$theoretical)

  y_d0d4 <- joined_quant_m_int$residual %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(0,.95)))

  type_labs <- c(np = "LF-np", pd = "3-pD", tp = "3-tP", tppd = "3x3 pD\n& tP")
  joined_quant_m_int$data_type  <- factor(joined_quant_m_int$data_type, levels = prefixes)
  joined_quant_m_int$quant_type <- factor(joined_quant_m_int$quant_type, levels = suffixes,
                                          labels = quant_labels)
  # ratios per ridge (the same for every data type x quant after intersecting)
  n_ridge <- max(table(joined_quant_m_int$data_type, joined_quant_m_int$quant_type))

  p <- ggplot(joined_quant_m_int, aes(x = residual, y = fct_rev(data_type))) +
    facet_grid(rows = vars(quant_type)) +
    geom_density_ridges2(
      aes(fill = data_type),
      alpha = 0.75,
      scale = 1.0,
      color = "black",
      linewidth = 0.25
    ) +
    stat_summary(
      aes(y = fct_rev(data_type)),
      fun = median,
      geom = "crossbar",
      width = 0.33,
      color = "black",
      linewidth = 0.25,
      middle.linewidth = 0.5
    ) +
    coord_cartesian(xlim = y_d0d4_val) +
    scale_y_discrete(labels = type_labs) +
    scale_fill_manual(values = c(np = "grey10", pd = "grey30", tp = "grey50", tppd = "grey70")) +
    labs(subtitle = bquote(italic(n)*" = "*.(scales::comma(n_ridge))*" precursor ratios"),
         x = bquote("|"*Log[2]*", precursor ratio error|"), y = NULL) +
    paper_theme_box() +
    theme(legend.position = "none", plot.subtitle = element_text(hjust = 1),
          axis.title.x = element_text(hjust = 1), plot.margin = margin(1, 1, 1, 1))

  fig4f <- paper_panel(p, outfile, width = 1.49, height = 2.52)
  invisible(list(fig4f = fig4f, data = joined_quant_m_int, counts = counts_tbl))
}


## ---- 9-plexDIA x 3-timePlex ratios (Fig. 5d) --------------------------------

level_quant <- function(normed_quant, grouper="prot", datatype = "protein",min_reps=2, extra="p",exclude_LF = F,
                        out_file = "Figs/Fig5/Fig5d.pdf"){
  paper_figure(out_file)
  
  normed_quant <- normed_quant %>%
    mutate(Type = case_when(
      timePlex == TRUE & Label != "LF" ~ "tp_pd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))
  
  if (!"ProteinName" %in% names(normed_quant)) {
    normed_quant$ProteinName <- normed_quant$prot
  }
  
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",ProteinName),"Human",
                                                                    ifelse(grepl("YEAST",ProteinName),"Yeast","remove"))) %>% filter(!grepl("remove",species))
  
  pD <- normed_quant[normed_quant$Type=="pd",] %>% filter(Rep == 2)
  tp_pd <- normed_quant[normed_quant$Type=="tp_pd",] %>% filter(Rep == 2)

  if (exclude_LF==F){
    np <- normed_quant[normed_quant$Type=="np",]
    joined <- rbind(tp_pd, pD, np)
  } else{
    joined <- rbind(tp_pd, pD)
  }

  normed_quant_avg_ints <- joined %>% dplyr::group_by(Type, .data[[grouper]], Sample) %>%
    dplyr::summarize("quant_avg" = mean(quant_val,na.rm=T)) %>% dplyr::rename("denom_sample"="Sample")
  
  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop") %>% ungroup()
  
  joined <- joined %>%
    dplyr::left_join(human_medians, by = "run_chan") %>%  
    dplyr::mutate(norm = norm-median_Human) %>%  
    dplyr::select(-median_Human) 

  joined_quant <- joined %>% group_by(Type, Sample,.data[[grouper]],species) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_reps-1)) %>% #require atleast certain number of reps
    dplyr::summarize(med_quant = median(norm,na.rm=T))

  dynamic_column <- grouper  
  
  formula_str <- paste(dynamic_column, "species", sep = "+")
  joined_quant_d <- dcast(
    joined_quant, 
    formula = as.formula(paste(formula_str, "~ Type+Sample")), 
    value.var = "med_quant"
  )
  
  joined_quant_d$pd_JM <- joined_quant_d$pd_J-joined_quant_d$pd_M
  
  joined_quant_d$tp_pd_JM <- joined_quant_d$tp_pd_J-joined_quant_d$tp_pd_M
  
  if (exclude_LF==F){
    joined_quant_d$np_JM <- joined_quant_d$np_J-joined_quant_d$np_M
    joined_quant_m <- joined_quant_d %>% dplyr::select(.data[[grouper]],species, 
                                                       np_JM, pd_JM, tp_pd_JM) %>% melt()
  } else{
    joined_quant_m <- joined_quant_d %>% dplyr::select(.data[[grouper]],species, 
                                                       pd_JM, tp_pd_JM) %>% melt()
  }
  
  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined_quant_m <- joined_quant_m %>%
    dplyr::left_join(human_medians, by = "variable") %>%  
    dplyr::mutate(value = value-median_Human) %>%  
    dplyr::select(-median_Human) %>% na.omit()
  
  joined_quant_m <- joined_quant_m %>% 
    dplyr::mutate(denom_sample = ifelse(grepl("J",variable),"J",
                                        ifelse(grepl("M",variable),"M", "K")))
  joined_quant_m$Type <- sub("_[^_]*$", "", joined_quant_m$variable)
  
  joined_quant_m <- joined_quant_m %>% dplyr::left_join(normed_quant_avg_ints, by =c("Type","denom_sample",paste0(grouper)))
  
  allcounts  <- joined_quant_m %>% dplyr::count(variable)

  ##### AB
  joined_AB <- joined_quant_m[grepl("JM", joined_quant_m$variable),]

  counts_AB <- allcounts[grepl("JM",allcounts$variable),]
  
  int_AB <- joined_AB %>% na.omit() %>% dplyr::count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  
  df_int <- data.frame(variable = "Intersected", n=nrow(int_AB))
  counts_AB <- rbind(counts_AB, df_int)
  
  # Fig. 5d: number of ratios per data type and the densities of the
  # intersected ratios
  keep_types <- if (exclude_LF) c("pd_JM", "tp_pd_JM") else c("np_JM", "pd_JM", "tp_pd_JM")
  counts_AB$variable <- factor(counts_AB$variable, levels = c(keep_types, "Intersected"))
  joined_AB_int <- joined_AB %>% dplyr::filter(.data[[grouper]]%in%int_AB[[grouper]])
  joined_AB_int$variable <- factor(joined_AB_int$variable, levels = c(keep_types, "Intersected"))
  
  bar_labs  <- c(np_JM = "LF np", pd_JM = "9-pD", tp_pd_JM = "9-pD &\n3-tP", Intersected = "Intersect")
  bar_fills <- c(np_JM = "grey45", pd_JM = "#A6A6A6", tp_pd_JM = "#D9D9D9", Intersected = "#F0F0F0")
  type_labs <- c(np_JM = "LF no-plex", pd_JM = "9-plexDIA", tp_pd_JM = "9-plexDIA\n3-timePlex")
  level_lab <- ifelse(grouper == "seqcharge", "Precursor", "Protein")
  
  counts_AB_plot <- ggplot(counts_AB, aes(x = variable, y = n, fill = variable)) +
    geom_col(colour = "black", linewidth = 0.25, width = 0.9) +
    geom_text(aes(label = scales::comma(n)), vjust = -0.3, size = paper_geom_text(6)) +
    scale_fill_manual(values = bar_fills) +
    scale_x_discrete(labels = bar_labs) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = paste0(level_lab, " ratios (A/D)")) +
    coord_cartesian(clip = "off") +              # counts wider than their bar stay whole
    paper_theme() +
    theme(legend.position = "none", axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  joined_AB_int$type <- factor(type_labs[as.character(joined_AB_int$variable)], levels = type_labs[keep_types])
  
  # Compute the 0.25% and 99.75% quantiles of the plotted variable
  y_limits <- quantile(joined_AB_int$value, probs = c(0.0025, 0.9975), na.rm = TRUE)
  
  denseplot <- ratio_density_plot(joined_AB_int, facet = "type", exp_FC = -2,
                                  ylab = bquote(Log[2]*", A/D"), legend = TRUE) +
    scale_y_continuous(breaks = seq(-6, 6, 3)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(y_limits[1], y_limits[2])) +
    labs(x = "Density") +
    theme(legend.position.inside = c(0.48, 0.97),        # small legend and 7-pt strips, as in the manuscript
          legend.text = element_text(size = 6), legend.key.size = unit(6.5, "pt"),
          strip.text = element_text(size = 7))
  
  # the rotated bar labels hang into the bottom margin instead of pushing the
  # density axis title down
  p <- patchwork::wrap_plots(patchwork::free(counts_AB_plot, type = "space", side = "b"), denseplot,
                             ncol = 2, widths = c(0.6, 1.6)) +
    patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 16, 0)))
  fig5d <- paper_panel(p, out_file, width = 3.08, height = 2.05)

  invisible(list(fig5d = fig5d))
}

