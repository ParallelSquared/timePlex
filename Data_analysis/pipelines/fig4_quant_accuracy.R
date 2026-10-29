# Fig. 4: quantitative accuracy, plus Extended Data Figs. 5-6
# (the 'Investigating and correcting biases' section of the Results)

list(

  ## ---- Extended Data Fig. 5: column-specific peak widths and coverage ---------
  tar_target(FWHM_implementations, Jmod_FWHMs(uPAC=benchmarking_dat_uPAC$Amount_8, IO=benchmarking_dat_IO$Amount_8), format = "qs"),
  tar_target(coverage_IO_columnspecific, columns_specific_coverage(benchmarking_dat_IO$Amount_8, 
                                                                   extra="IO", use_best_channel = T, filter_level=0.5,
                                                                   title = "3-timePlex, separate emitters\n(IonOpticks Aurora)",
                                                                   out_file = "Figs/ExtDataFig5/ExtDataFig5b.pdf"), format = "qs"),
  tar_target(coverage_uPAC_columnspecific, columns_specific_coverage(benchmarking_dat_uPAC$Amount_8, 
                                                                     extra="uPAC", use_best_channel = T, filter_level=0.5,
                                                                     title = "3-timePlex, single emitter\n(Thermo \u00b5PAC Neo)",
                                                                     out_file = "Figs/ExtDataFig5/ExtDataFig5c.pdf"), format = "qs"),

  ## ---- Extended Data Fig. 6a-c: timePlex vs no-plex ratios against TIC and RT ----
  tar_target(jd0422_mzml, "/Users/jasonderks/Benchmarking/All_07202026/JD0422.mzML", format = "file"),
  tar_target(jd0431_mzml, "/Users/jasonderks/Benchmarking/All_07202026/JD0431.mzML", format = "file"),
  tar_target(tic_jd0422, read_ms1_tic(jd0422_mzml)),
  tar_target(tic_jd0431, read_ms1_tic(jd0431_mzml)),
  tar_target(ids_jd0422_tc2, read_jmod_ids(file.path(all_fpath_IO, "JD0422_results"), time_channel = 2)),
  tar_target(ids_jd0431,     read_jmod_ids(file.path(all_fpath_IO, "JD0431_results"))),
  tar_target(paired_jd0422_tc2_jd0431, pair_tp_noplex(ids_jd0422_tc2, ids_jd0431, tic_jd0422, tic_jd0431)),
  tar_target(fig_tp_noplex, plot_tp_noplex_figure(paired_jd0422_tc2_jd0431)),
  tar_target(fig_tp_noplex_pdf,
             save_tp_noplex_figure(fig_tp_noplex, "Figs/ExtDataFig6/ExtDataFig6_full.pdf"),
             format = "file"),

  ## ---- Fig. 4a: single run before vs after RT correction ----------------------
  tar_target(plot_quants_IO_MS1, single_run_quant(normed_quant_plexArea, use_best_channel = F,
                                                                      quant_col="plex_Area",
                                                  fig_run = "JD0422", fig_chan = "JD0422_2_LF",
                                                  fig_comparison = "BC"),format = "qs"),

  ## ---- Fig. 4b: separate vs single emitter ------------------------------------
  tar_target(plot_single_multiEmitter, quant_multi_single_emitter(normed_quant_plexArea,
                                                                  normed_quant_plexArea_uPAC,
                                                                  min_num_reps=2, alphascaler=1.5),format = "qs"),

  ## ---- Fig. 4c-e: LF no-plex, 3-plexDIA, 3-timePlex and 3x3 -------------------
  tar_target(quant_benchmarking_plexArea, quant_all_comparisons_v2(normed_quant_plexArea,
                                                                                 datatype = "precursor_plexArea",
                                                                                 alphascaler=1.5, min_num_reps=2), format = "qs"),

  ## ---- Fig. 4f: MS1 vs MS2 ----------------------------------------------------
  tar_target(compared_MS1_MS2, compare_quants(dat1=normed_quant_plexArea,
                                              dat2=normed_quant_CU_Area,
                                              quant1="plex_Area_rtnorm",
                                              quant2="CU_Area_rtnorm", min_num_reps=2, grouper="seqcharge",
                                              outfile="Figs/Fig4/Fig4f.pdf"), format="qs"),

  ## ---- assembled figures (Figs/<figure>/<figure>_full.pdf) -----------------------
  tar_target(fig4_full, paper_fig4(plot_quants_IO_MS1$fig4a, plot_single_multiEmitter$fig4b,
                                   quant_benchmarking_plexArea$fig4ce, compared_MS1_MS2$fig4f), format = "file"),
  tar_target(ed5_full, paper_ed5(FWHM_implementations$ed5a, coverage_IO_columnspecific$panel,
                                 coverage_uPAC_columnspecific$panel), format = "file"),

  ## ---- Fig. 4 v2: samples per run above c-e, f stretched to match (R/fig4_v2_functions.R) ----
  tar_target(quant_benchmarking_plexArea_v2, quant_all_comparisons_v3(normed_quant_plexArea,
                                                                      alphascaler = 1.5, min_num_reps = 2), format = "qs"),
  tar_target(fig4_full_v2, paper_fig4_v2(plot_quants_IO_MS1$fig4a, plot_single_multiEmitter$fig4b,
                                         quant_benchmarking_plexArea_v2$fig4ce, compared_MS1_MS2$fig4f), format = "file"),
  # Extended Data Fig. 7: the precursor-ratio counts the v2 c-e no longer shows
  tar_target(ed7_full, paper_ed7(quant_benchmarking_plexArea_v2$counts), format = "file"),

  ## ---- Extended Data Fig. 8: B/C accuracy across the gradient, per column, one 3x3 run (R/ed8_gradient_functions.R) ----
  tar_target(ed8_ratios_rt, ed8_ratios(normed_quant_plexArea, run = "JD0413"), format = "qs"),
  tar_target(ed8_full, paper_ed8(ed8_ratios_rt), format = "file")
)

