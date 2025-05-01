
source("Functions.R")


library(pacman)
print(Sys.time())

pacman::p_load(eulerr, ggalt, ggbeeswarm, OrgMassSpecR, robustbase, plyr,
               stringi, stringr, tidyr, patchwork, ggplot2, reshape2,
               Cairo, gridGraphics, data.table, devtools, diann, gridExtra, ggpubr,
               ggpointdensity, viridis, scales, matrixStats, dplyr, targets,
               MASS, flux, psych, sva, ggExtra, gprofiler2, NCmisc, ggrepel, ggridges, rstatix, dtplyr,
               lawstat,lme4, lsmeans, boot,msEmpiRe,rsvd,enviPat,qs,PupillometryR,gtools,stats, plotly,
               htmlwidgets,htmltools,gghalves,bayestestR,distributions,mzR,openMSfile,xcms,progressr,purrr,ggh4x,qs2,
               MSnbase,sn
)


############ BULK
all_fpath_IO <- "timePlex_IO_03262025"

all_fpath_uPAC <- "timePlex_uPAC_03262025"

lib_path <- "singlefragment_HY_03132025.tsv"

############ Single cell
sc_fpath <- "SC_03262025"

SC_mzml_fpath <- "JD0439.mzML"

############ Carryover during timePlex loadings
carryover_fpath <- "Checking_columnCarryover"

############ META
meta_SC <- "timeplex_SCmeta_03192025.tsv"
meta_path <- "timeplex_meta_03152025.tsv"

IO_features_fpath <- "Features_IO"
uPAC_features_fpath <- "Features_uPAC"
  
############ LIB
lib_LF_arabidopsis_path <- "LF_HY_Arab_top12_v1p9.tsv"


######### missing species data (H, HY, HY)
missingspec_fpath <- "MissingSpec_03262025"
missingspec_mzml_fpath <- "JD0483.mzML"

######### overlapping species sequence
species_overlap_fpath <-  "Overlapped_sequences.tsv"

####### single column (no timplex implementation)
single_column_fpath <- "Single_column_03262025"

##### comparing impact of spec lib on high-plex data
dat_mTRAQ_lib <- "mTRAQ_LF_int"
dat_LF_lib <- "LF_mTRAQ_int"



source_targets <- function(path) {
  source(path)$value
}

targets_list <- c(
  source_targets("Benchmarking_targets.R"),
  source_targets("SingleCell_targets.R")
)

targets_list <- unlist(targets_list, recursive = FALSE)

targets_list

