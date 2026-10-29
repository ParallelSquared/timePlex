# Mixed-species (human/yeast) 3-plex benchmarks: shared inputs
# Used by Figs. 2-4 and Extended Data Figs. 3, 5 and 6.

list(

  ## ---- 'separate emitters' (IonOpticks columns) -------------------------------
  tar_target(dat_all_IO, Jmod_read_folder(all_fpath_IO, lib_path) %>% Jmod_add_ProteinName(lib_path), format="qs"),
  tar_target(benchmarking_dat_IO, pD_get_data_Jmod(dat_all_IO, meta_path, infer_sample=T, use_best_channel=T, filter_level=0.5),format = "qs"),
  tar_target(removed_overlapping_species, rm_overlaping_sequences(benchmarking_dat_IO$Amount_8, species_overlap_fpath),format = "qs"),
  tar_target(normed_quant, prec_normalizeQuant(removed_overlapping_species, quant_col="coeff_ChannelUnique_CommonApex_rtnorm",
                                                                            use_best_channel = T, filter_level=0.5,
                                               compute_unadjusted_quant = F),format = "qs"),
  tar_target(normed_quant_CU_Area, prec_normalizeQuant(removed_overlapping_species, quant_col="coeff_ChannelUnique_Area_rtnorm",
                                               use_best_channel = T, filter_level=0.5,
                                               compute_unadjusted_quant = F),format = "qs"),
  tar_target(normed_quant_plexArea,prec_normalizeQuant(removed_overlapping_species, quant_col="plex_Area_rtnorm", 
                                               use_best_channel = T, filter_level=0.5,
                                               compute_unadjusted_quant = T),format = "qs"),

  ## ---- 'single emitter, merged flow' (uPAC columns) ---------------------------
  tar_target(dat_all_uPAC, Jmod_read_folder(all_fpath_uPAC, lib_path) %>% Jmod_add_ProteinName(lib_path), format="qs"),
  tar_target(benchmarking_dat_uPAC, pD_get_data_Jmod(dat_all_uPAC, meta_path, infer_sample=T, use_best_channel=T, filter_level=0.5),format = "qs"),
  tar_target(removed_overlapping_species_uPAC, rm_overlaping_sequences(benchmarking_dat_uPAC$Amount_8, species_overlap_fpath),format = "qs"),
  tar_target(normed_quant_plexArea_uPAC, prec_normalizeQuant(removed_overlapping_species_uPAC, quant_col="plex_Area_rtnorm", 
                                                       use_best_channel = T, filter_level=0.5),format = "qs")
)

