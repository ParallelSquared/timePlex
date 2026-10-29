# Fig. 3: timePlex coverage, plus Extended Data Figs. 3-4

list(

  ## ---- Fig. 3a: missing-species run searched with a human, yeast and arabidopsis library ----
  tar_target(MissingSpec, Jmod_read_folder(missingspec_fpath, lib_LF_arabidopsis_path) %>% Jmod_add_library(lib_LF_arabidopsis_path), format="qs"),
  tar_target(MissingSpec_rm_mixedSpecies, rm_overlaping_sequences(MissingSpec, species_overlap_fpath),format = "qs"),
  tar_target(HistQuant_MissingSpec, Missing_species_hist_quant(MissingSpec_rm_mixedSpecies, extra="", use_best_channel=T, filter_level=0.5), format="qs"),

  ## ---- Fig. 3b: TIC and XICs of the missing-species run -----------------------
  tar_target(H_XICs, extra_MS1_mzml_tp(missingspec_mzml_fpath, target_mz_first = 478.7857, z=2, 
                                                 mztol=0.0075, iso=14, rt_min=14.8, rt_max=23.2,
                                       tp1=29.54530,
                                       tp2=33.76385,
                                       tp3=38.18150,
                                       target_mz_1st = 435.2601,
                                                 ppm=10,extra="H",
                                       title = "GAPDH (Human)", peptide = "VIPELNGK", line_color = "#5BBFD3",
                                       out_file = "Figs/Fig3/Fig3b_human_XIC.pdf"), format="qs"),
  tar_target(Y_XICs, extra_MS1_mzml_tp(missingspec_mzml_fpath, target_mz_first = 574.8343, z=2, 
                                       mztol=0.0075, iso=14, rt_min=14.8, rt_max=23.2,
                                       tp1=53.71556,
                                       tp2=57.91969,
                                       tp3=62.71254,
                                       target_mz_1st = 574.8343,
                                       ppm=10,extra="Y",
                                       title = "TDH1 (Yeast)", peptide = "VVDLIEYVAK", line_color = "#E8836B",
                                       out_file = "Figs/Fig3/Fig3b_yeast_XIC.pdf"), format="qs"),
  tar_target(TIC_HY, jusTIC(missingspec_mzml_fpath,                                         
                                tp1=29.60982,
                                tp2=33.86009,
                                tp3=38.25386,
                                tp1_2=53.71556,
                                tp2_2=57.91969,
                                tp3_2=62.71254,
                                extra="HY_bulk", out_file = "Figs/Fig3/Fig3b_TIC.pdf"),format = "qs"),

  ## ---- Fig. 3d: 'separate emitters' -------------------------------------------
  tar_target(coverage_IO_bestchannel, tp_coverage(benchmarking_dat_IO$Amount_8, extra="IO_", use_best_channel = T, filter_level=0.5),format = "qs"),

  ## ---- Fig. 3e, Extended Data Fig. 3: 'single emitter, merged flow' -----------
  tar_target(coverage_uPAC_bestchannel, tp_coverage(benchmarking_dat_uPAC$Amount_8, extra="uPAC_", use_best_channel = T, filter_level=0.5,
                                                     out_files = c(prec = "Figs/ExtDataFig3/ExtDataFig3a.pdf",
                                                                   prot = "Figs/ExtDataFig3/ExtDataFig3b.pdf"),
                                                     width = 3.545, height = 2.8),format = "qs"),

  ## ---- Extended Data Fig. 4: optimized (mTRAQ) vs unoptimized (LF) library ----
  tar_target(searched_mTRAQ_lib_dat, Jmod_read_folder(dat_mTRAQ_lib, lib_path) %>% Jmod_add_ProteinName(lib_path), format="qs"),
  tar_target(searched_mTRAQ_lib_dat2, pD_get_data_Jmod(searched_mTRAQ_lib_dat, meta_path, infer_sample=T, use_best_channel = T, filter_level=0.5),format = "qs"),
  tar_target(searched_LF_lib_dat, Jmod_read_folder(dat_LF_lib, lib_path) %>% Jmod_add_ProteinName(lib_path), format="qs"),
  tar_target(searched_LF_lib_dat2, pD_get_data_Jmod(searched_LF_lib_dat, meta_path, infer_sample=T, use_best_channel = T, filter_level=0.5),format = "qs"),
  tar_target(compare_coverage_libs, compare_coverage(searched_mTRAQ_lib_dat2$Amount_8, searched_LF_lib_dat2$Amount_8),format = "qs"),

  ## ---- assembled figures (c is a table; the zoom wedges in b are drawn by hand) ----
  tar_target(fig3_full, paper_fig3(HistQuant_MissingSpec$fig3a, H_XICs$xic, Y_XICs$xic, TIC_HY$tic,
                                   coverage_IO_bestchannel$prec, coverage_IO_bestchannel$prot), format = "file"),
  tar_target(ed3_full, paper_ed3(coverage_uPAC_bestchannel$prec, coverage_uPAC_bestchannel$prot), format = "file"),
  tar_target(ed4_full, paper_ed4(compare_coverage_libs$ed4a, compare_coverage_libs$ed4b), format = "file")
)

