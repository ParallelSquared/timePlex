
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
  
)
