# Implementing timePlex: Fig. 2f and Extended Data Figs. 1-2

list(

  ## ---- Fig. 2f: library vs fine-tuned RTs -------------------------------------
  tar_target(plot_JD0319_RTs, pD_plot_RTs_TP(dat = benchmarking_dat_uPAC$Amount_8, run="JD0319",
                                             filter_level=0.5, use_best_channel=T, lib_fpath=lib_path_RT,
                                             tp_labels = c("\u03940 min", "\u03945 min", "\u039411 min")),format = "qs"),

  ## ---- Extended Data Fig. 1: carryover between timePlex channels --------------
  tar_target(dat_all_carryover, Jmod_read_folder(carryover_fpath, lib_path,distinct=F), format="qs"),
  tar_target(carryover_result, compare_carryover(dat_all_carryover, use_best_channel=T, filter_level=0.5, meta=meta_path), format="qs"),

  ## ---- Extended Data Fig. 2: 3-column vs single-column spray ------------------
  tar_target(dat_all_IO_singlecolumn, Jmod_read_folder(single_column_fpath, lib_path) %>% Jmod_add_ProteinName(lib_path), format="qs"),
  tar_target(benchmarking_dat_IO_singlecolumn, pD_get_data_Jmod(dat_all_IO_singlecolumn, meta_path, infer_sample=T, use_best_channel=T, filter_level=0.5),format = "qs"),
  tar_target(coverage_IO_imp, tp_coverage_flow(benchmarking_dat_IO_singlecolumn$Amount_8, extra="", filter_level=0.5),format = "qs"),

  ## ---- assembled figures (Figs/<figure>/<figure>_full.pdf) -----------------------
  tar_target(fig2_full, paper_fig2(plot_JD0319_RTs$fig2f), format = "file"),
  tar_target(ed1_full, paper_ed1(carryover_result[grep("^ed1", names(carryover_result))]), format = "file"),
  tar_target(ed2_full, paper_ed2(coverage_IO_imp$ed2a, coverage_IO_imp$ed2b, coverage_IO_imp$ed2c), format = "file")
)

