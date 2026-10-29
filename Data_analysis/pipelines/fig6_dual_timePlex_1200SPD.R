# Fig. 6 and Extended Data Fig. 10: dual-timePlex at >1,200 samples/day

list(

  ## ---- data: run order and dual-plex A/B setup from the mzML headers, JMod results ----
  tar_target(meta_DualTP, DualTP_build_metadata(astralDualTP_fpath_raw, plex_period = 2), format = "qs"),
  tar_target(dat_all_DualTP, Jmod_read_folder(astralDualTP_fpath, lib_DualTP, distinct=T), format = "qs"),
  tar_target(dat_DualTP,
             DualTP_annotate(dat_all_DualTP, meta_DualTP, sample_map = astralDualTP_samplemap, use_best_channel=T, filter_level=0.5) %>%
               DualTP_add_ProteinName(lib_DualTP), format = "qs"),

  ## ---- Fig. 6c-f: QC multipanel -----------------------------------------------
  tar_target(qc_multipanel_DualTP_CUCA_4each,
             DualTP_QC_multipanel(dat_DualTP,
                                  runs = c("JD1189","JD1190","JD1191","JD1192",
                                           "JD1231","JD1232","JD1233","JD1234"),
                                  binwidth = 1, extra ="MS2CommonApex_allRatios_dualTP_4each_60secbins",
                                  quant_col = "coeff_ChannelUnique_CommonApex_rtnorm"), format = "qs"),

  ## ---- Extended Data Fig. 10: FWHM per column, IDs per channel and set --------
  tar_target(counts_FWHM_dual_tP, prec_prot_counts(dat_DualTP,
                                                   runs = c("JD1189","JD1190","JD1191","JD1192",
                                                            "JD1231","JD1232","JD1233","JD1234")), format = "qs"),

  ## ---- assembled figures (Fig. 6a is a schematic, 6b a photo) ---------------------
  tar_target(fig6_full, paper_fig6(qc_multipanel_DualTP_CUCA_4each$fig6cf), format = "file"),
  tar_target(ed10_full, paper_ed10(counts_FWHM_dual_tP$ed10a, counts_FWHM_dual_tP$ed10b, counts_FWHM_dual_tP$ed10c),
             format = "file")
)

