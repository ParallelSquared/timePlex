
# analyzing the 9-plexDIA & 3-timePlex (27-plex) data.

list(
  
  # read data
  tar_target(dat_all_3x9_IO, Jmod_read_folder(all_3x9_fpath_IO, lib_path_PSMtag_1, distinct=T, filtered=F), format="qs"),
  tar_target(benchmarking_dat_IO_3x9, pD_get_data_Jmod(dat_all_3x9_IO, meta_path, infer_sample=F),format = "qs"), 
  
  # coverage
  tar_target(coverage_IO_bestchannel_3x9, tp_coverage_3x9(benchmarking_dat_IO_3x9$Amount_1, extra="3x9_IO", use_best_channel = T, filter_level=0.01, plex3x9=T),format = "qs"), 
  tar_target(removed_overlapping_species_3x9, rm_overlaping_sequences(benchmarking_dat_IO_3x9$Amount_1, species_overlap_fpath),format = "qs"), 
  
  # prepare quant
  tar_target(quantified_SR0_2_3x9_new, pD_ms1_quant(removed_overlapping_species_3x9, adjacent_scans=0, dontapplyLF=T, strict=F),format = "qs"), 
  tar_target(normed_quant_IO_MS1_bestfit_SR0_3x9_coeff, prec_normalizeQuant(quantified_SR0_2_3x9_new, quant_col="coeff", 
                                                                              use_best_channel = T, filter_level=0.01),format = "qs"), 
  
  # protein level, require at least 2 b-ions per precursor if it is non-K C-terminating.
  tar_target(maxLFQ_prots_strict_coeff, pD_jmod_maxlfq(normed_quant_IO_MS1_bestfit_SR0_3x9_coeff,
                                                 group.header="prot", id.header = "seqcharge", quantity.header = "quant_val",
                                                 minimum_b_count=2),format = "qs"),
  tar_target(maxLFQ_prots_strict_m_coeff, pD_melt_jmod_MaxLFQ(maxLFQ_prots_strict_coeff, meta_path),format = "qs"),
  tar_target(prot_normed_strict_coeff, prot_normalizeQuant(maxLFQ_prots_strict_m_coeff,quant_col="value"),format = "qs"),
  tar_target(protein_level_MS2Quant_3x9, Protein_level_quant(prot_normed_strict_coeff, grouper="prot", datatype = "protein", min_reps=2),format = "qs")

  # plotting XICs
  tar_target(PSMtag_XICs_abundant_yeast, extra_MS1_mzml_v2(mzML_fpath_3tP_9pD, target_mz_first = 771.4280, z=3, 
                                                 mztol=0.0075, iso=14, rt_min=14.8, rt_max=23.2,
                                                 tp1=16.10419,
                                                 tp2=20.66682,
                                                 tp3=24.71370,
                                                 target_mz_1st = 664.3248,
                                                 target_mz_2nd = 665.6626,
                                                 target_mz_3rd = 666.9987,                              
                                                 target_mz_4th = 668.3382,
                                                 target_mz_5th = 669.6760,
                                                 target_mz_6th = 671.0122,
                                                 target_mz_7th = 672.3474,
                                                 target_mz_8th = 673.6852,
                                                 target_mz_9th = 675.0172,
                                                 pD_plex = 9,
                                                 ppm=15,extra="PSMtag_yeast", mycolor="#e8a792",SC=F), format = "qs"),
  
  tar_target(PSMtag_XICs_abundant_human, extra_MS1_mzml_v2(mzML_fpath_3tP_9pD, target_mz_first = 771.4280, z=3, 
                                                           mztol=0.0075, iso=14, rt_min=14.8, rt_max=23.2,
                                                           tp1=23.74140,
                                                           tp2=28.47098,
                                                           tp3=32.12,
                                                           target_mz_1st = 593.2900,
                                                           target_mz_2nd = 594.6279,
                                                           target_mz_3rd = 595.9641,                              
                                                           target_mz_4th = 597.3035,
                                                           target_mz_5th = 598.6413,
                                                           target_mz_6th = 599.9775,
                                                           target_mz_7th = 601.3127,
                                                           target_mz_8th = 602.6505,
                                                           target_mz_9th = 603.9824,
                                                           pD_plex = 9,
                                                           ppm=15,extra="PSMtag_human",mycolor="#8dbfcd",SC=F), format = "qs"),
  
  tar_target(justTIC_PSMtag, jusTIC(mzML_fpath_3tP_9pD,                                         
                                tp1=16.10419,
                                tp2=20.66682,
                                tp3=24.71370,
                                tp1_2=23.74140,
                                tp2_2=28.47098,
                                tp3_2=32.12,
                                extra="abundant",log_trans=F), format = "qs")
  
)
