# Fig. 5e-h: single cells by 9x3 plexDIA & timePlex

list(

  ## ---- data -------------------------------------------------------------------
  tar_target(sc_dat, Jmod_read_folder(sc_fpath, lib_path) %>% Jmod_add_ProteinName(lib_path), format="qs"),
  tar_target(sc_dat2, pD_get_data_Jmod(sc_dat, meta_SC, infer_sample=F, filter=T, use_best_channel=T, filter_level=0.5, non_uniform_name=T), format = "qs"),

  ## ---- Fig. 5e-f: quality control and IDs per cell ----------------------------
  tar_target(sc_filt, sc_quality_control(sc_dat2$Amount_SC, number_sd = 2), format = "qs"),
  tar_target(sc_qc_plot, plot_sc_qc(sc_dat2$Amount_SC, number_sd = 2), format = "qs"),
  tar_target(sc_plot, plot_counts(sc_dat2$Amount_SC %>% dplyr::filter(run_chan%in%sc_filt$run_chan)), format = "qs"),

  ## ---- protein quantification: LC drift, maxLFQ, normalization, imputation, batch correction ----
  tar_target(joined_SC_bulk, pD_join_and_intersect(sc_dat2$Amount_SC, sc_dat2$Amount_bulk,intersect=F), format = "qs"),
  tar_target(joined_SC_bulk_corrected, pD_spline_correct_forLC(joined_SC_bulk), format = "qs"),
  tar_target(SC_bulk_maxLFQ, diann_maxlfq(joined_SC_bulk_corrected %>%
                                            dplyr::filter(!grepl("Neg", Sample) & plex_Area_rtnorm > 1 & Protein_Qvalue < 0.01 &
                                                            (run_chan %in% sc_filt$run_chan | grepl("ulk", run_chan))) %>%
                                            pD_prec_norm(Quant = "plex_Area_rtnorm"),
                                        group.header="prot", id.header = "seqcharge", quantity.header = "norm_prec_lin"), format = "qs"),
  tar_target(SC_bulk_maxLFQ_m, pD_melt_MaxLFQ(SC_bulk_maxLFQ, meta_SC), format = "qs"),
  tar_target(SC_bulk_norm_MLFQ, pD_prot_norm(SC_bulk_maxLFQ_m, Quant="norm_prot"), format = "qs"),
  tar_target(SC_bulk_imputed_MLFQ_log2, pD_impute_log2(SC_bulk_norm_MLFQ, missing.prot.frac = 0.9, missing.cell.frac = 0.9, k=10), format = "qs"),
  tar_target(BC_MLFQ, pD_batchCorrect_labs_LC(SC_bulk_imputed_MLFQ_log2, meta_SC), format = "qs"),

  ## ---- Fig. 5g: single cell vs bulk fold changes ------------------------------
  tar_target(SC_bulk_cors, cor_bulk_SC(BC_MLFQ, meta_SC, min_frac_quant=0.5), format = "qs"),

  ## ---- Fig. 5h: PCA -----------------------------------------------------------
  tar_target(SC_bulk_PCAs, pD_weighted_PCA_timePlex(BC_MLFQ, min_sc_frac = 0.01), format = "qs"),

  ## ---- data tables for sharing: proteins/precursors x samples and sample metadata ----
  tar_target(sc_data_tables, sc_export_tables(sc_dat2, sc_filt, joined_SC_bulk_corrected, SC_bulk_maxLFQ, BC_MLFQ,
                                              meta_SC, "timeplex_SC_acquisition_10272026.tsv"), format = "file"),

  ## ---- assembled Fig. 5 (a is a table; the zoom wedges in b are drawn by hand) ----
  tar_target(fig5_full, paper_fig5(PSMtag_XICs_abundant_yeast$xic, PSMtag_XICs_abundant_human$xic, justTIC_PSMtag$tic,
                                   coverage_IO_bestchannel_3x9$prec, coverage_IO_bestchannel_3x9$prot,
                                   prec_level_CUArea_3x9$fig5d, sc_qc_plot$fig5e, sc_plot$fig5f,
                                   SC_bulk_cors$fig5g, SC_bulk_PCAs$fig5h), format = "file")
)

