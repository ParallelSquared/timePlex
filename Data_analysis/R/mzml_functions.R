################################################################################
# Raw data from mzML files: extracted and total ion chromatograms
################################################################################


## ---- XICs (Fig. 3b: 3-timePlex; Fig. 5b: 9-plexDIA x 3-timePlex) ------------

extra_MS1_mzml_tp <- function(mzml_file, target_mz_first = 534.3383, z=2, 
                              mztol=0.0075, iso=14, rt_min=31.6, rt_max=40.6,
                              tp1=31.93026,
                              tp2=36.03937,
                              tp3=40.28072,
                              target_mz_1st = 435.2601,
                              ppm=10,
                              extra="",
                              title = NULL, peptide = NULL, line_color = "black",
                              tp_labels = c("\u03940 min", "\u03944 min", "\u03948 min"),
                              out_file = paste0("Smoothed_Summed_Isotope_Signals_Bulk_TP_", extra, ".pdf")){
  paper_figure(out_file)
  
  mzdata <- openMSfile(mzml_file)
  hdr <- header(mzdata)
  ms1_scans <- hdr[hdr$msLevel == 1, ]
  
  mztol <- target_mz_1st*ppm/1000000

  target_mz_values <- c(target_mz_1st)
  mz_tolerance <- target_mz_1st*ppm/1000000
  
  charge <- z
  
  compute_isotopes <- function(mono_mz, charge, num_isotopes = 1) {
    return((mono_mz) + (0:num_isotopes * 1.00335 / charge))
  }
  
  isotope_mz_list <- lapply(target_mz_values, compute_isotopes, charge = charge)
  
  extract_binned_signal <- function(scan, isotope_mz_set, mzdata, mz_tolerance) {
    peaks <- peaks(mzdata, scan)
    mz_values <- peaks[, 1]
    intensities <- peaks[, 2]
    
    summed_intensity <- sum(sapply(isotope_mz_set, function(mz_center) {
      sum(intensities[abs(mz_values - mz_center) <= mz_tolerance])
    }))
    
    return(summed_intensity)
  }
  
  binned_data <- data.frame(
    scan = ms1_scans$seqNum,
    retention_time = ms1_scans$retentionTime / 60  
  )
  
  for (i in seq_along(target_mz_values)) {
    target_mass <- target_mz_values[i]
    isotopes <- isotope_mz_list[[i]]
    
    summed_intensities <- sapply(ms1_scans$seqNum, function(scan) {
      extract_binned_signal(scan, isotopes, mzdata, mz_tolerance)
    })
    
    binned_data[[paste0("sum_", round(target_mass, 4))]] <- summed_intensities
  }
  
  binned_long <- melt(binned_data, id.vars = c("retention_time", "scan"),
                      variable.name = "mass_category", value.name = "Intensity")

  binned_data <- data.frame(
    scan = ms1_scans$seqNum,
    retention_time = ms1_scans$retentionTime / 60 
  )
  
  for (i in seq_along(target_mz_values)) {
    target_mass <- target_mz_values[i]
    isotopes <- isotope_mz_list[[i]]
    
    summed_intensities <- sapply(ms1_scans$seqNum, function(scan) {
      extract_binned_signal(scan, isotopes, mzdata, mz_tolerance)
    })
    
    binned_data[[paste0("sum_", round(target_mass, 4))]] <- summed_intensities
  }

  binned_long <- melt(binned_data, id.vars = c("retention_time", "scan"),
                      variable.name = "mass_category", value.name = "Intensity")

  binned_long <- binned_long %>% dplyr::mutate("timePlex" = ifelse(retention_time<tp1+0.33 & retention_time>tp1-0.33, "1",
                                                                   ifelse(retention_time<tp2+0.33 & retention_time>tp2-0.33, "2",
                                                                          ifelse(retention_time<tp3+0.33 & retention_time>tp3-0.33, "3","none")))) %>%
    dplyr::filter(timePlex!="none")

  binned_long$timePlex <- factor(binned_long$timePlex, levels = c("3", "2", "1"))
  
  binned_long <- binned_long %>%
    dplyr::mutate(interpolated = FALSE)

  smooth_intensity <- function(df) {
    df <- df %>% dplyr::arrange(retention_time)  
    if (nrow(df) > 3) {  
      spline_fit <- spline(df$retention_time, df$Intensity, n = 50)  
      df_interp <- data.frame(
        retention_time = spline_fit$x, 
        y_interpolated = spline_fit$y, 
        mass_category = unique(df$mass_category), 
        timePlex = unique(df$timePlex),
        interpolated = TRUE  
      )
      return(df_interp)
    } else {
      return(df %>% dplyr::mutate(y_interpolated = Intensity, interpolated = FALSE))  
    }
  }
  
  #smoothing per mass_category & celltype
  binned_long_smooth <- binned_long %>% 
    dplyr::group_by(mass_category, timePlex) %>%  
    group_split() %>% 
    lapply(smooth_intensity) %>% 
    bind_rows()
  
  binned_long <- binned_long %>% dplyr::rename(Intensity_real = Intensity)
  
  binned_long_smooth <- full_join(binned_long_smooth, binned_long, 
                                  by = c("retention_time", "mass_category", "timePlex", "interpolated")) %>%
    dplyr::mutate(y_interpolated = ifelse(is.na(y_interpolated), Intensity_real, y_interpolated))
  
  binned_long_smooth <- binned_long_smooth %>%
    
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4)  

  binned_long_smooth <- binned_long_smooth %>%
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4)  
  binned_long_smooth$y_interpolated[binned_long_smooth$y_interpolated<0] <- 0
  binned_long_smooth$timePlex <- factor(binned_long_smooth$timePlex, levels = c("1", "2", "3"))
  binned_long_smooth$tp_lab <- factor(tp_labels[as.integer(as.character(binned_long_smooth$timePlex))],
                                      levels = tp_labels)
  
  # Fig. 3b: one XIC per timePlex channel
  p <- ggplot(binned_long_smooth, 
              aes(x = retention_time, y = y_interpolated, group = mass_category))+
    geom_hline(yintercept = 0, linewidth = 0.35) +
    facet_grid(cols = vars(tp_lab), scales = "free") +
    geom_point(data = binned_long_smooth %>% dplyr::filter(interpolated == FALSE), 
               colour = "black", size = 1.4, alpha = 0.9) +
    geom_line(colour = line_color, linewidth = 0.5) +
    scale_x_continuous(breaks = xic_breaks, labels = function(x) sprintf("%.1f", x)) +
    labs(x = "Retention Time (min)", y = "MS1 intensity", title = title,
         subtitle = if (!is.null(peptide)) bquote(.(peptide)^.(paste0("+", z)))) +
    paper_theme_box(base_size = 7, frame = 0.12) +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(0, 0, 1, 0)),
          plot.subtitle = element_text(size = 9, hjust = 0.5, margin = margin(0, 0, 3, 0)),
          panel.spacing.x = unit(0.15, "in"))
  
  invisible(list(xic = paper_panel(p, out_file, width = 2.95, height = 1.45)))
}

xic_breaks <- function(lim){
  b <- pretty(lim, n = 3)
  b <- b[b >= lim[1] & b <= lim[2]]
  b[c(1, length(b))]
}

# Two x breaks at ~20% and ~80% of a narrow XIC window, on 0.05 min steps, so
# the labels of neighbouring panels cannot run into each other (Fig. 5b)
xic_breaks_inner <- function(lim){
  round((lim[1] + c(0.2, 0.8) * diff(lim)) / 0.05) * 0.05
}

extra_MS1_mzml_v2 <- function(mzml_file, target_mz_first = 534.3383, z=2, 
                              mztol=0.0075, iso=14, rt_min=31.6, rt_max=40.6,
                              tp1=31.93026,
                              tp2=36.03937,
                              tp3=40.28072,
                              target_mz_1st = 534.3383,
                              target_mz_2nd = 536.3418,
                              target_mz_3rd = 538.3453,
                              target_mz_4th = NA,
                              target_mz_5th = NA,
                              target_mz_6th = NA,
                              target_mz_7th = NA,
                              target_mz_8th = NA,
                              target_mz_9th = NA,
                              pD_plex = 3,
                              ppm=10,
                              extra="",mycolor="#e8a792",
                              title = NULL, peptide = NULL,
                              tp_labels = c("\u03940 min", "\u03944 min", "\u03948 min"),
                              channel_labels = paste0("\u0394", seq(0, 16, 2), " Da"),
                              channel_samples = c("D", "A", "B", "C", "D", "C", "A", "B", "C"),
                              channel_amounts = c("4x", "1x", "2x", "3x", "4x", "3x", "1x", "2x", "3x"),
                              out_file = paste0("Smoothed_Summed_Isotope_Signals_normalized", extra, ".pdf")){
  paper_figure(out_file)
  
  mzdata <- openMSfile(mzml_file)
  hdr <- header(mzdata)
  ms1_scans <- hdr[hdr$msLevel == 1, ]
  
  mztol <- target_mz_2nd*ppm/1000000

  if(pD_plex==3){
    target_mz_values <- c(target_mz_1st, target_mz_2nd, target_mz_3rd)
  } else{
    target_mz_values <- c(target_mz_1st, target_mz_2nd, target_mz_3rd, target_mz_4th,
                          target_mz_5th,
                          target_mz_6th,
                          target_mz_7th,
                          target_mz_8th,
                          target_mz_9th)
  }

  mz_tolerance <- target_mz_2nd*ppm/1000000
  
  charge <- z
  
  compute_isotopes <- function(mono_mz, charge, num_isotopes = 1) {
    return((mono_mz) + (0:num_isotopes * 1.00335 / charge))
  }
  
  isotope_mz_list <- lapply(target_mz_values, compute_isotopes, charge = charge,num_isotopes = 1)
  
  extract_binned_signal <- function(scan, isotope_mz_set, mzdata, mz_tolerance) {
    peaks <- peaks(mzdata, scan)
    mz_values <- peaks[, 1]
    intensities <- peaks[, 2]
    
    summed_intensity <- sum(sapply(isotope_mz_set, function(mz_center) {
      sum(intensities[abs(mz_values - mz_center) <= mz_tolerance])
    }))
    
    return(summed_intensity)
  }
  
  binned_data <- data.frame(
    scan = ms1_scans$seqNum,
    retention_time = ms1_scans$retentionTime / 60  # Convert to minutes
  )
  
  for (i in seq_along(target_mz_values)) {
    target_mass <- target_mz_values[i]
    isotopes <- isotope_mz_list[[i]]
    
    summed_intensities <- sapply(ms1_scans$seqNum, function(scan) {
      extract_binned_signal(scan, isotopes, mzdata, mz_tolerance)
    })
    
    binned_data[[paste0("sum_", round(target_mass, 4))]] <- summed_intensities
  }

  binned_long <- melt(binned_data, id.vars = c("retention_time", "scan"),
                      variable.name = "mass_category", value.name = "Intensity")

  binned_long <- binned_long %>% dplyr::mutate("timePlex" = ifelse(retention_time<tp1+0.13 & retention_time>tp1-0.13, "1",
                                                                   ifelse(retention_time<tp2+0.13 & retention_time>tp2-0.13, "2",
                                                                          ifelse(retention_time<tp3+0.13 & retention_time>tp3-0.13, "3","none")))) %>%
    dplyr::filter(timePlex!="none")

  binned_long$timePlex <- factor(binned_long$timePlex, levels = c("3", "2", "1"))
  
  if(pD_plex==3){
    df_map <- data.frame("mapper" = c(paste0("sum_",target_mz_1st,"_1"), paste0("sum_",target_mz_2nd,"_1"), paste0("sum_",target_mz_3rd,"_1"),
                                      paste0("sum_",target_mz_1st,"_2"), paste0("sum_",target_mz_2nd,"_2"), paste0("sum_",target_mz_3rd,"_2"),
                                      paste0("sum_",target_mz_1st,"_3"), paste0("sum_",target_mz_2nd,"_3"), paste0("sum_",target_mz_3rd,"_3")),
                         "celltype" = c("U937","K562","K562","U937","K562","U937","K562","U937","U937"))
   } else{
    df_map <- data.frame("mapper" = c(paste0("sum_",target_mz_1st,"_1"), paste0("sum_",target_mz_2nd,"_1"), paste0("sum_",target_mz_3rd,"_1"),
                                      paste0("sum_",target_mz_4th,"_1"), paste0("sum_",target_mz_5th,"_1"), paste0("sum_",target_mz_6th,"_1"),
                                      paste0("sum_",target_mz_7th,"_1"), paste0("sum_",target_mz_8th,"_1"), paste0("sum_",target_mz_9th,"_1"),
                                      
                                      paste0("sum_",target_mz_1st,"_2"), paste0("sum_",target_mz_2nd,"_2"), paste0("sum_",target_mz_3rd,"_2"),
                                      paste0("sum_",target_mz_4th,"_2"), paste0("sum_",target_mz_5th,"_2"), paste0("sum_",target_mz_6th,"_2"),
                                      paste0("sum_",target_mz_7th,"_2"), paste0("sum_",target_mz_8th,"_2"), paste0("sum_",target_mz_9th,"_2"),
                                      
                                      paste0("sum_",target_mz_1st,"_3"), paste0("sum_",target_mz_2nd,"_3"), paste0("sum_",target_mz_3rd,"_3"),
                                      paste0("sum_",target_mz_4th,"_3"), paste0("sum_",target_mz_5th,"_3"), paste0("sum_",target_mz_6th,"_3"),
                                      paste0("sum_",target_mz_7th,"_3"), paste0("sum_",target_mz_8th,"_3"), paste0("sum_",target_mz_9th,"_3")),
                         "celltype" = c("4x","1x","2x","3x","4x","3x","1x","2x","3x",
                                        "4x","1x","2x","3x","4x","3x","1x","2x","3x",
                                        "4x","1x","2x","3x","4x","3x","1x","2x","3x"))
   }
  
  binned_long$mapper <- paste0(binned_long$mass_category,"_",binned_long$timePlex)
  binned_long <- binned_long %>% left_join(df_map, by =c("mapper"="mapper"))
  
  binned_long <- binned_long %>%
    dplyr::mutate(interpolated = FALSE)  
  
  smooth_intensity <- function(df) {
    df <- df %>% dplyr::arrange(retention_time)  
    if (nrow(df) > 3) {  
      spline_fit <- spline(df$retention_time, df$Intensity, n = 50)  
      df_interp <- data.frame(
        retention_time = spline_fit$x, 
        y_interpolated = spline_fit$y, 
        mass_category = unique(df$mass_category), 
        celltype = unique(df$celltype),
        timePlex = unique(df$timePlex),
        interpolated = TRUE  
      )
      return(df_interp)
    } else {
      return(df %>% dplyr::mutate(y_interpolated = Intensity, interpolated = FALSE))  
    }
  }
  
  binned_long_smooth <- binned_long %>% 
    dplyr::group_by(mass_category, celltype, timePlex) %>%  
    group_split() %>% 
    lapply(smooth_intensity) %>% 
    bind_rows()
  
  binned_long <- binned_long %>% dplyr::rename(Intensity_real = Intensity)
  binned_long_smooth <- full_join(binned_long_smooth, binned_long, 
                                  by = c("retention_time", "mass_category", "celltype", "timePlex", "interpolated")) %>%
    dplyr::mutate(y_interpolated = ifelse(is.na(y_interpolated), Intensity_real, y_interpolated))
  
  binned_long_smooth <- binned_long_smooth %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4)  

  binned_long_smooth <- binned_long_smooth %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4,  
           color_group = case_when(  
             celltype == "K562" ~ norm_int,   
             celltype == "U937" ~ -norm_int   
           ))  
  binned_long_smooth$y_interpolated[binned_long_smooth$y_interpolated<0] <- 0
  binned_long_smooth$timePlex <- factor(binned_long_smooth$timePlex, levels = c("1", "2", "3"))

  binned_long_smooth_tp_norm <- binned_long_smooth %>% dplyr::group_by(timePlex) %>%
    dplyr::mutate("tp_norm" = y_interpolated/max(y_interpolated,na.rm=T)) %>% 
    dplyr::group_by(timePlex,mass_category) %>%
    dplyr::mutate("AUC" = sum(tp_norm,na.rm=T)) %>% ungroup() %>%
    dplyr::group_by(timePlex) %>%
    dplyr::mutate(AUC_rel = AUC/min(AUC)) %>%
    ungroup()
  
  binned_long_smooth_tp_norm_dis <- binned_long_smooth_tp_norm %>% dplyr::distinct(timePlex, mass_category,.keep_all = T)
  
  # Fig. 5b: rows = plexDIA channels, columns = timePlex channels. Each row's
  # sample and yeast amount ("D (4x)") is written right of the strips, the
  # letter in the sample color: the black label carries a phantom letter and
  # the colored letter a phantom remainder, so the two overlap exactly.
  d <- binned_long_smooth_tp_norm
  ch <- as.integer(d$mass_category)          # channel order = target m/z order
  d$offset_lab <- factor(channel_labels[ch], levels = channel_labels)
  d$tp_lab     <- factor(tp_labels[as.integer(as.character(d$timePlex))], levels = tp_labels)
  
  side <- data.frame(offset_lab = factor(channel_labels, levels = channel_labels),
                     tp_lab     = factor(tp_labels[length(tp_labels)], levels = tp_labels),
                     colour     = paper_cols$samples4[channel_samples],
                     black      = sprintf('phantom(bold("%s"))*" (%s)"', channel_samples, channel_amounts),
                     letter     = sprintf('bold("%s")*phantom(" (%s)")', channel_samples, channel_amounts))
  side_text <- function(label, colour) {
    geom_text(data = side, aes(x = Inf, y = 0.5, label = .data[[label]]), colour = colour, inherit.aes = FALSE,
              parse = TRUE, angle = -90, vjust = -1.9, size = paper_geom_text(7))
  }
  
  p <- ggplot(d, aes(x = retention_time, y = tp_norm, group = mass_category)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    facet_grid(rows = vars(offset_lab), cols = vars(tp_lab), scales = "free_x") +
    geom_line(linewidth = 0.6, colour = mycolor) +
    geom_point(data = d %>% dplyr::filter(interpolated == FALSE),
               colour = "black", size = 0.6, alpha = 0.9) +
    side_text("black", "black") +
    geom_text(data = side, aes(x = Inf, y = 0.5, label = letter, colour = colour), inherit.aes = FALSE,
              parse = TRUE, angle = -90, vjust = -1.9, size = paper_geom_text(7)) +
    scale_colour_identity() +
    scale_x_continuous(breaks = xic_breaks_inner, labels = function(x) sprintf("%.2f", x)) +
    coord_cartesian(clip = "off") +
    labs(x = "Retention Time (min)", y = "MS1 intensity (timePlex normalized)", title = title,
         subtitle = if (!is.null(peptide)) bquote(.(peptide)^.(paste0("+", z)))) +
    paper_theme_box(base_size = 7) +              
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 8),
          strip.text.y = element_text(angle = -90),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          panel.spacing.x = unit(4, "pt"), panel.spacing.y = unit(2, "pt"),
          plot.margin = margin(3, 15, 3, 3))
  
  invisible(list(xic = paper_panel(p, out_file, width = 1.92, height = 3.95)))
}

## ---- TICs (Figs. 3b and 5b) -------------------------------------------------

jusTIC <- function(mzml_file, tp1=31.93026,
                   tp2=36.03937,
                   tp3=40.28072,
                   tp1_2=31.93026,
                   tp2_2=36.03937,
                   tp3_2=40.28072,
                   extra="",
                   band_colors = "mediumpurple1",
                   out_file = paste0(extra, "_TIC_xaxis.pdf"), width = 6.7, height = 1.45){
  paper_figure(out_file)
  
  mzdata <- openMSfile(mzml_file)
  hdr <- header(mzdata)
  ms1_scans <- hdr[hdr$msLevel == 1, ]

  if (nrow(ms1_scans) > 0) {
    ms1_intensity_data <- data.frame(RT = ms1_scans$retentionTime, TotalIntensity = NA)
    
    for (i in seq_len(nrow(ms1_scans))) {
      peaks_data <- peaks(mzdata, ms1_scans$acquisitionNum[i])  
      ms1_intensity_data$TotalIntensity[i] <- sum(peaks_data[, 2], na.rm = TRUE) 
    }

    rect_data <- data.frame(
      ymin = -Inf, ymax = Inf,  
      xmin = c(tp1 - 0.13, tp2 - 0.13, tp3 - 0.13, tp1_2 - 0.13, tp2_2 - 0.13, tp3_2 - 0.13),
      xmax = c(tp1 + 0.13, tp2 + 0.13, tp3 + 0.13, tp1_2 + 0.13, tp2_2 + 0.13, tp3_2 + 0.13)
    )

    # Figs. 3b and 5b: TIC with the XIC windows marked. The panel keeps its
    # x range and band centres so the figure can draw the zoom wedges.
    rect_data$fill <- rep_len(band_colors, nrow(rect_data))
    p <- ggplot(ms1_intensity_data, aes(y = (TotalIntensity), x = RT/60)) +  
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
                alpha = 1, inherit.aes = FALSE) +  
      geom_line(color = "black", linewidth = 0.13) +
      scale_fill_identity() +
      labs(
        y = "TIC",
        x = "Retention Time (min)") +
      paper_theme(base_size = 8, frame = 0.12) +
      theme(axis.text = element_text(size = 7))
    tic <- paper_panel(p, out_file, width = width, height = height)
    tic$x_range <- ggplot_build(p)$layout$panel_params[[1]]$x.range
    tic$bands   <- (rect_data$xmin + rect_data$xmax) / 2
    return(invisible(list(tic = tic)))
  }
}

