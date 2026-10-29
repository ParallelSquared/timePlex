################################################################################
# timePlex (Derks et al. 2026) analysis pipeline. Run with targets::tar_make().
#
#   R/          functions
#   pipelines/  targets, one file per manuscript figure
################################################################################

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
invisible(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))

print(Sys.time())

pacman::p_load(eulerr, ggalt, ggbeeswarm, OrgMassSpecR, robustbase, plyr,
               stringi, stringr, tidyr, patchwork, ggplot2, reshape2,
               Cairo, gridGraphics, data.table, devtools, diann, gridExtra, ggpubr,
               ggpointdensity, viridis, scales, matrixStats, dplyr, targets,
               MASS, flux, psych, sva, ggExtra, gprofiler2, NCmisc, ggrepel, ggridges, rstatix, dtplyr,
               lawstat,lme4, lsmeans, boot,msEmpiRe,rsvd,enviPat,qs,PupillometryR,gtools,stats, plotly,
               htmlwidgets,htmltools,gghalves,bayestestR,mzR,xcms,progressr,purrr,ggh4x,qs2,
               MSnbase,sn,forcats, quantreg, splines
)

tar_source("R")
tar_option_set(error = "trim",   # a failed target only skips what depends on it; everything else still builds
               memory = "transient", garbage_collection = TRUE)   # drop inputs after each target (36 GB RAM)


## ---- Inputs shared by the bulk benchmarks ------------------------------------
meta_path <- "timeplex_meta_10082026.tsv"
lib_path <- "ProteinGroup_ProteinName_map_10232026.tsv"  # ProteinGroup -> ProteinName, merged from the search libraries
lib_path_RT <- "singlefragment_HY_03132025.tsv"  # only for the original library RT (Tr_recalibrated) in Fig. 2f
species_overlap_fpath <- "/Users/jasonderks/Overlapped_sequences.tsv"

## JMod-timePlex_beta re-search, one folder per pipeline input below.
new_search <- "/Users/jasonderks/Benchmarking/JMod_timePlex_beta_10232026"

## ---- benchmark_data.R --------------------------------------------------------
all_fpath_IO <- file.path(new_search, "IO")
all_fpath_uPAC <- file.path(new_search, "uPAC")

## ---- fig2_implementation.R ---------------------------------------------------
carryover_fpath <- file.path(new_search, "New_Carryover")
single_column_fpath <- file.path(new_search, "New_tP_implementation")

## ---- fig3_coverage.R ---------------------------------------------------------
missingspec_fpath <- file.path(new_search, "MissingSpecies")
missingspec_mzml_fpath <- "/Users/jasonderks/Downloads/MissingSpecies/JD0483.mzML"
lib_LF_arabidopsis_path <- "/Users/jasonderks/Downloads/LF_HY_Arab_top12_noX_v1p9.tsv"
dat_mTRAQ_lib <- file.path(new_search, "mTRAQ_LF_int")
dat_LF_lib <- file.path(new_search, "LF_mTRAQ_int")

## ---- fig5_9x3_plexDIA_timePlex.R ---------------------------------------------
all_3x9_fpath_IO <- file.path(new_search, "IO_3x9")
mzML_fpath_3tP_9pD <- "/Users/jasonderks/Downloads/raw_3x9_benchmarking/JD0582_re.mzML"

## ---- fig5_single_cell.R ------------------------------------------------------
sc_fpath <- file.path(new_search, "SC_original_optimum")
meta_SC <- "timeplex_SCmeta_allRuns_08262026.tsv"

## ---- fig6_dual_timePlex_1200SPD.R --------------------------------------------
astralDualTP_fpath <- file.path(new_search, "Astral_1300SPD")
astralDualTP_fpath_raw <- "/Volumes/Lab/JD/Raw_data/timePlex_revisions/Astral_1300SPD_raw"
lib_DualTP <- "/Users/jasonderks/Benchmarking/58kPrecursor_HY_lib_PSMTag_07152026.tsv"
astralDualTP_samplemap <- data.frame(plex_channel = seq(0, 16, 2),
                                     Sample = c("A","B","C","B","C","A","C","A","B"))

## ---- Targets, one file per manuscript figure ---------------------------------
list(
  source("pipelines/benchmark_data.R")$value,
  source("pipelines/fig1_concept.R")$value,
  source("pipelines/fig2_implementation.R")$value,
  source("pipelines/fig3_coverage.R")$value,
  source("pipelines/fig4_quant_accuracy.R")$value,
  source("pipelines/fig5_9x3_plexDIA_timePlex.R")$value,
  source("pipelines/fig5_single_cell.R")$value,
  source("pipelines/fig6_dual_timePlex_1200SPD.R")$value
)
