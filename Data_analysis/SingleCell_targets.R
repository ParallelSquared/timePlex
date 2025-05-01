
# Here we analyze the single-cell data

list(
  
  #######################################################
  #######  Single cell and bulk K and U cells:  #########
  #######################################################

  tar_target(throughputplot, throughput_plot(), format="qs"),
  tar_target(sc_dat, Jmod_read_folder(sc_fpath, lib_path, distinct=F, filtered=F), format="qs"),
  
  tar_target(sc_dat2, pD_get_data_Jmod(sc_dat, meta_SC, infer_sample=F), format = "qs"),
  tar_target(bulk, compute_prec_prot_qvalues_v2(sc_dat2$Amount_100x, toPlot=F), format = "qs"),
  tar_target(sc, compute_prec_prot_qvalues_v2(sc_dat2$Amount_SC, toPlot=T), format = "qs"),

  tar_target(bulk_ms1, pD_ms1_quantINT(bulk$df, adjacent_scans=1), format = "qs"), 
  tar_target(sc_ms1, pD_ms1_quant(sc$df, adjacent_scans=1), format = "qs"), 
  tar_target(sc_ms1_filt, sc_quality_control(sc_ms1, min_counts = 250), format = "qs"),
  tar_target(sc_plot, plot_counts(sc, sc_ms1_filt, filt=F), format = "qs"),
  
  tar_target(bulk_normedQuant, prec_normalizeQuant(bulk_ms1, quant_col="MS1_bestfit", use_best_channel = T, filter_level=0.01), format = "qs"), 
  tar_target(sc_normedQuant, prec_normalizeQuant(sc_ms1_filt, quant_col="MS1_Int", use_best_channel = T, filter_level=0.01), format = "qs"), 

  tar_target(output_MSEmpPrep, 
             pD_MsEmpire_prep_jmod(bulk_normedQuant, BR=c(1), conds=c("K562","U937")), format = "qs"),
  
  tar_target(MsEmpire_bulk_KU,  pD_MsEmpire_run_specific_jmod(output_MSEmpPrep), format = "qs"), 
  
######
  tar_target(SC_maxLFQ, diann_maxlfq(sc_normedQuant, group.header="prot", id.header = "seqcharge", quantity.header = "quant_val"), format = "qs"), 
  tar_target(Bulk_maxLFQ, diann_maxlfq(bulk_normedQuant, group.header="prot", id.header = "seqcharge", quantity.header = "quant_val"), format = "qs"), 

  tar_target(dat_SC_maxLFQ_m, pD_melt_MaxLFQ(SC_maxLFQ, meta_SC), format = "qs"),
  tar_target(dat_Bulk_maxLFQ_m, pD_melt_MaxLFQ(Bulk_maxLFQ, meta_SC), format = "qs"),

  # normalize, impute, batch correct
  tar_target(SC_norm_MLFQ, pD_prot_norm(dat_SC_maxLFQ_m, Quant="norm_prot"), format = "qs"), 
  tar_target(SC_imputed_MLFQ_log2, pD_impute_log2(SC_norm_MLFQ, missing.prot.frac = 0.95, missing.cell.frac = 0.95, k=5), format = "qs"),

  # normalize, impute, batch correct
  tar_target(Bulk_norm_MLFQ, pD_prot_norm(dat_Bulk_maxLFQ_m, Quant="norm_prot"), format = "qs"), 
  tar_target(Bulk_imputed_MLFQ_log2, pD_impute_log2(Bulk_norm_MLFQ, missing.prot.frac = 0.95, missing.cell.frac = 0.95, k=2), format = "qs"),

  tar_target(joined_bulk_SC, join_post_imp(SC_imputed_MLFQ_log2, Bulk_imputed_MLFQ_log2), format = "qs"),

  tar_target(BC_MLFQ, pD_batchCorrect_labs_LC(joined_bulk_SC, meta_SC), format = "qs"),
  tar_target(SC_bulk_cors, cor_bulk_SC(BC_MLFQ, meta_SC), format = "qs"),

  tar_target(SC_XICs, extra_MS1_mzml_v2(SC_mzml_fpath, target_mz_first = 534.3383, z=2, 
                                        mztol=0.0075, iso=14, rt_min=31.6, rt_max=40.6,
                                        tp1=31.93026,
                                        tp2=36.03937,
                                        tp3=40.28072,
                                        target_mz_1st = 534.3383,
                                        target_mz_2nd = 536.3418,
                                        target_mz_3rd = 538.3453,
                                        ppm=10,extra="dif"),format = "qs"),


  tar_target(SC_XICs_abundant, extra_MS1_mzml_v2(SC_mzml_fpath, target_mz_first = 478.7857, z=2, 
                                        mztol=0.0075, iso=14, rt_min=14.8, rt_max=23.2,
                                        tp1=15.03870,
                                        tp2=19.12919,
                                        tp3=22.94236,
                                        target_mz_1st = 478.7857,
                                        target_mz_2nd = 480.7896,
                                        target_mz_3rd = 482.7928,
                                        ppm=10,extra="abundant"), format = "qs"),

  tar_target(SC_justTIC, jusTIC(SC_mzml_fpath,                                         
                                tp1=31.93026,
                                tp2=36.03937,
                                tp3=40.28072,
                                                 tp1_2=15.03870,
                                                 tp2_2=19.12919,
                                                 tp3_2=22.94236,
                                                 extra="abundant",log_trans=T), format = "qs"),


  tar_target(SC_bulk_PCAs, pD_weighted_PCA_timePlex(BC_MLFQ), format = "qs")

)
