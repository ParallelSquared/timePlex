
# Here we analyze the benchmarking data


list(

  #######################################################
  ###############  Conceptualization:  ###########
  #######################################################
  
  tar_target(cost, pD_cost(d0=778.64, d4=2348.04, d8=2343.62, numcolumns=3, pack10columns=10160,runspercolumn=1000), format="qs"),
  tar_target(cost_pd_only, pD_cost_plexDIA(d0=778.64, allothers=2348.04,midplex=10,highplex=100), format="qs"),
  
  #######################################################
  ##########  Human yeast quant evaluations:  ###########
  #######################################################
 
  ########### investigating carryover during timePlex sample loading
  tar_target(dat_all_carryover, Jmod_read_folder(carryover_fpath, lib_path,distinct=F, filtered=F), format="qs"),
  tar_target(carryover_result, compare_carryover(dat_all_carryover), format="qs"),
  
  ############ uPAC columns
  tar_target(dat_all_uPAC, Jmod_read_folder(all_fpath_uPAC, lib_path, distinct=F, filtered=F), format="qs"),
  tar_target(benchmarking_dat_uPAC, pD_get_data_Jmod(dat_all_uPAC, meta_path, infer_sample=F),format = "qs"), 
  tar_target(plot_JD0319_RTs, pD_plot_RTs_TP(dat = benchmarking_dat_uPAC$Amount_8, libpath = lib_path, run="JD0319"),format = "qs"), 
  tar_target(coverage_uPAC_bestchannel, tp_coverage(benchmarking_dat_uPAC$Amount_8, extra="uPAC",use_best_channel = T, filter_level=0.01),format = "qs"), 
  tar_target(removed_overlapping_species_uPAC, rm_overlaping_sequences(benchmarking_dat_uPAC$Amount_8, species_overlap_fpath),format = "qs"), 
  
  # prepare quant
  tar_target(quantified_SR1_uPAC_2, pD_ms1_quant(removed_overlapping_species_uPAC, adjacent_scans=1, dontapplyLF=T),format = "qs"), 
  tar_target(normed_quant_IO_MS1_bestfit_relaxed_SR1_uPAC_2, prec_normalizeQuant(quantified_SR1_uPAC_2, quant_col="MS1_bestfit", 
                                                                                 use_best_channel = T, filter_level=0.01),format = "qs"), 
  
  ############ IO columns
  tar_target(dat_all_IO, Jmod_read_folder(all_fpath_IO, lib_path,distinct=F, filtered=F), format="qs"),
  tar_target(benchmarking_dat_IO, pD_get_data_Jmod(dat_all_IO, meta_path, infer_sample=T),format = "qs"), 
  tar_target(coverage_IO_bestchannel, tp_coverage(benchmarking_dat_IO$Amount_8, extra="IO", use_best_channel = T, filter_level=0.01),format = "qs"), 
  tar_target(removed_overlapping_species, rm_overlaping_sequences(benchmarking_dat_IO$Amount_8, species_overlap_fpath),format = "qs"), 
  
  # prepare quant
  tar_target(quantified_SR1_2, pD_ms1_quant(removed_overlapping_species, adjacent_scans=1, dontapplyLF=T),format = "qs"), 
  tar_target(normed_quant_IO_MS1_bestfit_relaxed_SR1_2, prec_normalizeQuant(quantified_SR1_2, quant_col="MS1_bestfit", 
                                                                            use_best_channel = T, filter_level=0.01),format = "qs"), 
  
  # compare single vs multi-emitter approaches
  tar_target(plot_single_multiEmitter, quant_multi_single_emitter(normed_quant_IO_MS1_bestfit_relaxed_SR1_2, 
                                                                  normed_quant_IO_MS1_bestfit_relaxed_SR1_uPAC_2, 
                                                                  min_num_reps=2, alphascaler=1.5),format = "qs"), 
  
  tar_target(plot_quants_IO_MS1_bestfit_relaxed_SR1, single_run_quant(normed_quant_IO_MS1_bestfit_relaxed_SR1_2, use_best_channel = T, 
                                                                      quant_col="MS1_bestfit"),format = "qs"), 
  
  tar_target(quant_int_benchmarking_relaxed_SR1_MNR2_2, quant_all_comparisons_v2(normed_quant_IO_MS1_bestfit_relaxed_SR1_2, 
                                                                                 datatype = "precursor_relaxed_scanrad1_minNumReps2",
                                                                                 alphascaler=1.5, min_num_reps=2), format = "qs"), 
   
  #### compare MS1 and MS2
  tar_target(normed_quant_IO_relaxed_coeff, prec_normalizeQuant(quantified_SR1_2, quant_col="coeff", use_best_channel = T, filter_level=0.01),format = "qs"), 
  tar_target(compared_MS1_MS2, compare_quants(dat1=normed_quant_IO_MS1_bestfit_relaxed_SR1_2, 
                                              dat2=normed_quant_IO_relaxed_coeff,
                                              quant1="MS1_bestfit", 
                                              quant2="coeff", min_num_reps=2, grouper="seqcharge"), format="qs"),
  

  ############ IO columns... single column
  tar_target(dat_all_IO_singlecolumn, Jmod_read_folder(single_column_fpath, lib_path, distinct=F, filtered=F), format="qs"),
  tar_target(benchmarking_dat_IO_singlecolumn, pD_get_data_Jmod(dat_all_IO_singlecolumn, meta_path, infer_sample=T),format = "qs"), 
  tar_target(joined_implementations, join_implementations(benchmarking_dat_IO$Amount_8, benchmarking_dat_IO_singlecolumn$Amount_8),format = "qs"), 
  tar_target(coverage_IO_imp, tp_coverage_compare_implementation(joined_implementations, extra="IO_joined_implementations", use_best_channel = T, filter_level=0.01),format = "qs"), 
  
  ############ IO columns
  tar_target(normed_quant_IO, prec_normalizeQuant(removed_overlapping_species, quant_col="plexfitMS1", use_best_channel = F, filter_level=0.01),format = "qs"), 
  tar_target(plot_quants_IO_MS1, single_run_quant(normed_quant_IO, use_best_channel = F, quant_col="plexfitMS1"),format = "qs"), 
 
  tar_target(normed_quant_IO_coeff, prec_normalizeQuant(removed_overlapping_species, quant_col="coeff", use_best_channel = F, filter_level=0.01),format = "qs"), 
  tar_target(plot_quants_IO_coeff, single_run_quant(normed_quant_IO_coeff, use_best_channel = T, quant_col="coeff"),format = "qs"), 
  tar_target(plot_quants_IO_MS1_BestFit, single_run_quant(normed_quant_IO_coeff, use_best_channel = T, quant_col="MS1_bestfit"),format = "qs"), 
  tar_target(plot_quants_IO_MS1plexfit, single_run_quant(normed_quant_IO_coeff, use_best_channel = T, quant_col="MS1_bestfit"),format = "qs"), 
 
 
  # Chromatography
  tar_target(FWHM_implementations, Jmod_FWHMs(uPAC=benchmarking_dat_uPAC$Amount_8, IO=benchmarking_dat_IO$Amount_8),format = "qs"), 
  tar_target(coverage_uPAC_columnspecific, columns_specific_coverage(benchmarking_dat_uPAC$Amount_8, extra="uPAC_relaxed",use_best_channel = T, filter_level=0.01),format = "qs"), 
  tar_target(coverage_IO_columnspecific, columns_specific_coverage(benchmarking_dat_IO$Amount_8, extra="IO_relaxed",use_best_channel = T, filter_level=0.01),format = "qs"), 
 
 
  ############ Investigate how library affects high-plex data:
  # Searched with optimized library
  tar_target(searched_mTRAQ_lib_dat, Jmod_read_folder(dat_mTRAQ_lib, lib_path, distinct=F,filtered=F), format="qs"),
  tar_target(searched_mTRAQ_lib_dat2, pD_get_data_Jmod(searched_mTRAQ_lib_dat, meta_path, infer_sample=T),format = "qs"), 
 
  # Searched with unoptimized library
  tar_target(searched_LF_lib_dat, Jmod_read_folder(dat_LF_lib, lib_path, distinct=F,filtered=F), format="qs"),
  tar_target(searched_LF_lib_dat2, pD_get_data_Jmod(searched_LF_lib_dat, meta_path, infer_sample=T),format = "qs"), 
 
  # compare
  tar_target(compare_coverage_libs, compare_coverage(searched_mTRAQ_lib_dat2$Amount_8, searched_LF_lib_dat2$Amount_8),format = "qs"), 
 
 
  #######################################################
  ################   Missing species:   #################
  #######################################################
  tar_target(MissingSpec, Jmod_read_folder(missingspec_fpath, lib_LF_arabidopsis_path, distinct=T), format="qs"),
  tar_target(MissingSpec_rm_mixedSpecies, rm_overlaping_sequences(MissingSpec, species_overlap_fpath),format = "qs"), 

  tar_target(HistQuant_MissingSpec, Missing_species_hist_quant(MissingSpec_rm_mixedSpecies), format="qs"),

  tar_target(H_XICs, extra_MS1_mzml_tp(missingspec_mzml_fpath, target_mz_first = 478.7857, z=2, 
                                                 mztol=0.0075, iso=14, rt_min=14.8, rt_max=23.2,
                                       tp1=29.54530,
                                       tp2=33.76385,
                                       tp3=38.18150,
                                       target_mz_1st = 435.2601,
                                                 ppm=10,extra="H"), format="qs"),

  tar_target(Y_XICs, extra_MS1_mzml_tp(missingspec_mzml_fpath, target_mz_first = 574.8343, z=2, 
                                       mztol=0.0075, iso=14, rt_min=14.8, rt_max=23.2,
                                       tp1=53.71556,
                                       tp2=57.91969,
                                       tp3=62.71254,
                                       target_mz_1st = 574.8343,
                                       ppm=10,extra="Y"), format="qs"),
  
  tar_target(TIC_HY, jusTIC(missingspec_mzml_fpath,                                         
                                tp1=29.60982,
                                tp2=33.86009,
                                tp3=38.25386,
                                tp1_2=53.71556,
                                tp2_2=57.91969,
                                tp3_2=62.71254,
                                extra="HY_bulk"),format = "qs")
   
)

