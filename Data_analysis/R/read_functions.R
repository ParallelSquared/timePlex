################################################################################
# Reading and annotating JMod search results and run metadata
################################################################################
# Jmod_read_folder() reads <folder>/<Run>_*/filtered_IDs.csv for every run,
# drops decoys and, with distinct = TRUE, joins ProteinName etc. from the
# spectral library on untag_prec. time_channel indexes the timePlex column
# (0/1/2; 1 = middle).


## ---- JMod search results ----------------------------------------------------

Jmod_read_folder <- function(folder_path, lib_fpath,distinct=F){
  print("starting")
  file_list <- list.files(folder_path)
  data_list <- list()
  counter = 0
  for (file_name in file_list) {
    counter = counter + 1 
    file_path <- file.path(folder_path, file_name,'filtered_IDs.csv')
    data <- data.frame(fread(file_path, 
                             select = c("seq","z","is_decoy","time_channel","rt","lib_rt","stripped_seq","MS1_Int","untag_prec","channel","coeff", "protein",
                                        "plexfittrace_ps_all","plexfittrace_all", "MS1_Area", "all_ms1_iso0vals","mz",
                                        "Qvalue","Protein_Qvalue","Best_pDChannel_Qvalue","plexfitMS1","plexfitMS1_p", "Best_tP_specific_channelQvalue",
                                        "plex_Area",
                                        "plex_Area_rtnorm","MS1ConsensusApexAligned","coeff_ChannelUnique", "MS1_Area_rtnorm",
                                        "MS1ConsensusApexAligned_rtnorm", "coeff_rtnorm", "coeff_ChannelUnique_rtnorm", 
                                        "coeff_ChannelUnique_Area_rtnorm",
                                        "coeff_ChannelUnique_CommonApex_rtnorm", "MS1_FWHM",
                                        
                                        "rtnorm_factor_MS1", "rtnorm_factor_MS2"
                                        )))
    data <- data[data$is_decoy==F,]
    temp <- str_extract(file_name, ".+?(?=_re|_\\d{14})")
    
    data$Run <- temp
    data_list[[counter]] <- data
    print(paste0("Loaded ",counter,"/",length(file_list)," files"))
  }
  df <- do.call("rbind.fill", data_list)
  if(distinct==T){
    lib <- fread(lib_fpath)
    pepcol <- if ("PeptideSequence" %in% names(lib)) "PeptideSequence" else "StrippedPeptide"
    if (!pepcol %in% names(lib)) stop("Neither 'PeptideSequence' nor 'StrippedPeptide' found in lib")
    
    lib$seqcharge <- paste0(lib[[pepcol]], "_", lib$PrecursorCharge)
    lib <- lib %>% distinct(transition_group_id, ProteinName, .keep_all=T)
    df <- df %>% left_join(lib, by =c("untag_prec" = "seqcharge"))
  }

  return(df)
}

# ProteinName (for species) from a ProteinGroup -> ProteinName map, matched on the JMod protein column
Jmod_add_ProteinName <- function(df, map_fpath){
  map <- fread(map_fpath, select = c("ProteinGroup", "ProteinName"))
  df <- df %>% left_join(map, by = c("protein" = "ProteinGroup"))
  return(df)
}

# ProteinName and Proteotypic from a DIA-NN library, matched on the modified
# precursor as JMod writes it in untag_prec ("TLETANC(UniMod:4)MSSQTK_2"); the
# stripped-sequence join of Jmod_read_folder(distinct = TRUE) misses every modified precursor
Jmod_add_library <- function(df, lib_fpath){
  lib <- fread(lib_fpath, select = c("FullUniModPeptideName", "PrecursorCharge", "ProteinName", "Proteotypic", "decoy"))
  lib <- unique(lib[decoy == 0])
  lib$untag_prec <- paste0(lib$FullUniModPeptideName, "_", lib$PrecursorCharge)
  df %>% left_join(lib[, c("untag_prec", "ProteinName", "Proteotypic")], by = "untag_prec")
}

pD_filter <- function(dat, use_best_channel=T, filter_level=0.5){
  if (isTRUE(use_best_channel)) {
    # NOTE: BestChannel cut kept at the original hard-coded 0.01, not filter_level
    dat <- dat %>%
      dplyr::filter(Best_pDChannel_Qvalue < 0.01) %>%
      dplyr::filter(Qvalue < filter_level)
  } else {
    dat <- dat %>% dplyr::filter(Qvalue < 0.01)
  }
}

pD_get_data_Jmod <- function(dat, meta_path,infer_sample,filter=T,use_best_channel=T, filter_level=0.01, non_uniform_name=F){
  m <- fread(meta_path)
  dat <- dat %>% dplyr::rename("seqcharge" = "untag_prec")
  dat <- dat %>% dplyr::rename("prot" = "protein")
  if (!"time_channel" %in% colnames(dat)) {
    dat$time_channel <- NA
  }
  
  dat$channel[is.na(dat$time_channel)&(!grepl("mTRAQ|PSMtag_",dat$seq))] <- "NoPlex_LF"
  dat$channel[is.na(dat$channel) & is.na(dat$time_channel) & grepl("mTRAQ|PSMtag_", dat$seq)] <- ""
  
  dat$time_channel[is.na(dat$time_channel)] <- 0
  dat <- dat %>% dplyr::mutate("channel" = ifelse(grepl("_", channel),channel,
                                                  ifelse(grepl("mTRAQ|PSMtag_",seq), paste0("_",channel),paste0(channel,"_LF"))))
  
  dat <- dat %>% dplyr::group_by(Run, time_channel, seqcharge) %>% 
    dplyr::mutate("Best_pDChannel_Qvalue" = min(Qvalue,na.rm=T)) %>% ungroup()
  
  dat$run_chan <- paste0(dat$Run, "_", dat$channel)
  m$Column <- as.numeric(m$Column)
  
  if(non_uniform_name==F){
    m$Run <- sub("_.*", "", m$Run)
  }
  
  m <- m %>% dplyr::mutate("TP" = ifelse(timePlex==T,Column-1,
                                         ifelse(Label=="LF","NoPlex","")))
  m$run_chan <- paste0(m$Run, "_", m$TP,"_",m$Label)
  m$batch <- paste0(m$Label,"_",m$Column)
  m <- m %>% dplyr::select(-c("Run"))
  
  if (infer_sample == TRUE) {
    m <- m %>%
      dplyr::mutate(
        Sample = as.character(Sample),
        Sample = case_when(
          Sample == "" & V9 == "C1" & Label == "LF" ~ "A",
          Sample == "" & V9 == "C2" & Label == "LF" ~ "B",
          Sample == "" & V9 == "C3" & Label == "LF" ~ "C",
          
          Sample == "" & V9 == "D1" & Label == "0"  ~ "A",
          Sample == "" & V9 == "D1" & Label == "4"  ~ "B",
          Sample == "" & V9 == "D1" & Label == "8"  ~ "C",
          
          Sample == "" & V9 == "D2" & Label == "0"  ~ "C",
          Sample == "" & V9 == "D2" & Label == "4"  ~ "A",
          Sample == "" & V9 == "D2" & Label == "8"  ~ "B",
          
          Sample == "" & V9 == "D3" & Label == "0"  ~ "B",
          Sample == "" & V9 == "D3" & Label == "4"  ~ "C",
          Sample == "" & V9 == "D3" & Label == "8"  ~ "A",
          
          TRUE ~ Sample  
        )
      )
  }
  dat<-dat %>% inner_join(m, by =c("run_chan"="run_chan"))
  dat$File.Name <- dat$run_chan
  uniruns <- dat %>% distinct(File.Name,Sample)
  print(uniruns)
  unique_loadings <- unique(dat$Amount)
  if (filter == TRUE) {
    dat <- pD_filter(dat, use_best_channel=use_best_channel, filter_level=filter_level)
  }
  
  fin <- list()
  for(i in 1:length(unique_loadings)){
    temp <- dat[dat$Amount==unique_loadings[i],]
    dats <- paste0("Amount_",unique_loadings[i])
    fin[[dats]] <- temp
  }
  return(fin)
}

rm_overlaping_sequences <- function(dat, overlap_path){
  overlaps <- fread(overlap_path)
  print(nrow(dat))
  dat <- dat[dat$stripped_seq%!in%overlaps$PeptideSequence,]
  print(nrow(dat))
  
  return(dat)
}


## ---- Dual-timePlex run metadata (Fig. 6) ------------------------------------
## mzR cannot open Orbitrap Astral mzML files (its ProteoWizard rejects cvParam
## MS:1003378), so acquisition start times are read straight from the mzML
## header.

#' Acquisition start timestamp from the mzML header.
DualTP_mzml_start_time <- function(path, nchars = 2e6, max_nchars = 32e6) {
  pat <- 'startTimeStamp="([^"]+)"'
  repeat {
    txt <- tryCatch(readChar(path, nchars = nchars, useBytes = TRUE),
                    error = function(e) NA_character_)
    if (is.na(txt)) return(NA_character_)
    m <- regmatches(txt, regexpr(pat, txt, perl = TRUE))
    if (length(m)) return(sub(pat, "\\1", m, perl = TRUE))
    if (nchar(txt, type = "bytes") < nchars || nchars >= max_nchars)
      return(NA_character_)
    nchars <- nchars * 4
  }
}

#' ISO-8601 -> POSIXct UTC, tolerating "Z", "+HH:MM" and bare forms.
DualTP_parse_iso8601 <- function(x) {
  x <- trimws(x)
  x <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", x)
  out <- as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) out[bad] <- as.POSIXct(x[bad], format = "%Y-%m-%dT%H:%M:%OS%z", tz = "UTC")
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) out[bad] <- as.POSIXct(x[bad], format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  out
}

#' Scan every .mzML header for its acquisition start timestamp; order runs,
#' zero the clock, and assign the alternating dual-timePlex setup.
#'
#' @param plex_period 2 = strict A,B,A,B...   4 = A,A,B,B...
#' @param plex_csv    optional TSV/CSV with columns file, plex to override the
#'                    automatic alternation for any run (hand-edited exceptions).
#' @return data.table: file, path, Run, start_time (UTC), injection_index,
#'         t_start_hr, plex
DualTP_build_metadata <- function(mzml_dir, plex_period = 2L,
                                  plex_labels = c("timePlex A", "timePlex B"),
                                  plex_csv = NULL) {
  files <- list.files(mzml_dir, pattern = "\\.mzML$", full.names = TRUE,
                      ignore.case = TRUE)
  stopifnot(length(files) > 0)
  message(sprintf("Found %d mzML files in %s", length(files), mzml_dir))

  ts <- vapply(seq_along(files), function(i) {
    if (i %% 100 == 0) message(sprintf("  ...%d/%d", i, length(files)))
    DualTP_mzml_start_time(files[i])
  }, character(1))

  meta <- data.table(
    file       = basename(files),
    path       = files,
    Run        = str_extract(basename(files), "[^_.]+"),
    start_time = DualTP_parse_iso8601(ts),
    file_mtime = file.info(files)$mtime
  )
  n_bad <- sum(is.na(meta$start_time))
  if (n_bad > 0) {
    message(sprintf("WARNING: %d/%d files have no parseable startTimeStamp; using file mtime.",
                    n_bad, nrow(meta)))
    meta[is.na(start_time), start_time := as.POSIXct(file_mtime, tz = "UTC")]
  }

  setorder(meta, start_time, file)
  meta[, injection_index := seq_len(.N)]
  t0 <- min(meta$start_time)
  meta[, t_start_hr := as.numeric(difftime(start_time, t0, units = "hours"))]

  meta[, plex := plex_labels[((injection_index - 1L) %/% (plex_period %/% 2L)) %%
                               length(plex_labels) + 1L]]
  if (!is.null(plex_csv) && file.exists(plex_csv)) {
    ov <- fread(plex_csv)
    meta[ov, plex := i.plex, on = "file"]
    message("Applied plex overrides from ", plex_csv)
  }

  meta[, file_mtime := NULL]
  message(sprintf("%d runs | t0 = %s UTC | span = %.2f h | median gap = %.1f s",
                  nrow(meta), format(t0), max(meta$t_start_hr),
                  median(diff(as.numeric(meta$start_time)))))
  print(meta[, .N, by = plex])
  meta[]
}


## ---- Dual-timePlex annotation (Fig. 6) --------------------------------------

#' Rename to downstream conventions, attach acquisition metadata, and assign
#' Sample from the plexDIA mass-tag channel. `channel` strings like "d4",
#' "mTRAQ-4" or "4" all reduce to the number 4 before the join.
DualTP_annotate <- function(dat, meta, sample_map = NULL, use_best_channel=T, filter_level=0.01) {
  dat <- dat %>% dplyr::rename("seqcharge" = "untag_prec", "prot" = "protein")
  dat <- dat %>% dplyr::group_by(Run, time_channel, seqcharge) %>%
    dplyr::mutate("Best_pDChannel_Qvalue" = min(Qvalue,na.rm=T)) %>% ungroup()
  
  dat$run_chan <- paste0(dat$Run, "_", dat$channel)
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(Best_pDChannel_Qvalue<0.01) %>% dplyr::filter(Qvalue < filter_level)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<0.01)
  }
  
  meta <- as.data.frame(meta)[, c("Run", "plex", "injection_index", "t_start_hr")]
  dat <- dat %>% inner_join(meta, by = "Run")
  if (nrow(dat) == 0) stop("No overlap between result Run names and mzML Run names")
  
  ## channel strings carry the timePlex column too (e.g. "2_d16"), so grab only
  ## the mass-tag number: digits right after a "d", else the LAST digit run
  chan_str <- as.character(dat$channel)
  tag <- stringr::str_extract(chan_str, "(?i)(?<=d)\\d+")
  no_d <- is.na(tag)
  if (any(no_d)) tag[no_d] <- stringr::str_extract(chan_str[no_d], "\\d+(?=\\D*$)")
  dat$plex_channel <- suppressWarnings(as.numeric(tag))
  if (all(is.na(dat$plex_channel)))
    stop("Could not parse a mass-tag number out of `channel`; unique values were: ",
         paste(head(unique(dat$channel), 20), collapse = ", "))
  
  if (is.null(sample_map)) {
    sample_map <- data.frame(plex_channel = seq(0, 16, 2),
                             Sample = c("A","B","C","B","A","C","C","A","B"))
  }
  sample_map <- as.data.frame(sample_map)
  if (!"plex_channel" %in% colnames(sample_map) && "channel" %in% colnames(sample_map)) {
    sample_map$plex_channel <-
      suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(sample_map$channel))))
    sample_map$channel <- NULL
  }
  by_cols <- intersect(c("plex", "time_channel", "plex_channel"), colnames(sample_map))
  dat <- dat %>% left_join(sample_map[, c(by_cols, "Sample")], by = by_cols)
  
  if (anyNA(dat$Sample)) {
    bad <- unique(dat$plex_channel[is.na(dat$Sample)])
    stop("Channels with no Sample in sample_map: ", paste(bad, collapse = ", "),
         " -- sample_map covers: ", paste(sort(unique(sample_map$plex_channel)), collapse = ", "))
  }
  print(dat %>% dplyr::distinct(channel, time_channel, plex_channel, Sample) %>%
          dplyr::arrange(time_channel, plex_channel), n = 100)
  return(dat)
}

DualTP_add_ProteinName <- function(dat, lib_path) {
  lib_uni <- unique(fread(lib_path, select = c("ProteinGroup", "ProteinName")))
  lib_uni <- lib_uni[!is.na(ProteinGroup) & ProteinGroup != "" &
                       !is.na(ProteinName)  & ProteinName  != ""]
  
  ## exact group-string map (a group with conflicting names keeps them all)
  map_g <- lib_uni[, .(pn = paste(unique(ProteinName), collapse = ";")),
                   by = ProteinGroup]
  
  ## per-accession map from exploded groups
  map_a <- lib_uni[, .(acc = unlist(strsplit(ProteinGroup, ";", fixed = TRUE))),
                   by = ProteinName]
  map_a <- map_a[, .(pn = paste(unique(ProteinName), collapse = ";")), by = acc]
  
  pn_new <- map_g$pn[match(dat$prot, map_g$ProteinGroup)]
  miss <- is.na(pn_new)
  if (any(miss)) {
    first_acc <- sub(";.*", "", dat$prot[miss])
    pn_new[miss] <- map_a$pn[match(first_acc, map_a$acc)]
  }
  
  had <- if ("ProteinName" %in% colnames(dat)) dat$ProteinName else NA_character_
  dat$ProteinName <- ifelse(is.na(pn_new), had, pn_new)
  
  n <- length(dat$ProteinName)
  tagged <- grepl("HUMAN|YEAST", dat$ProteinName)
  message(sprintf("ProteinName from library: %s of %s rows matched (%.1f%%); %s rows (%.1f%%) now carry a HUMAN/YEAST tag",
                  format(sum(!is.na(pn_new)), big.mark = ","), format(n, big.mark = ","),
                  100 * mean(!is.na(pn_new)),
                  format(sum(tagged, na.rm = TRUE), big.mark = ","),
                  100 * mean(tagged, na.rm = TRUE)))
  un <- unique(dat$prot[is.na(dat$ProteinName) | dat$ProteinName == ""])
  if (length(un))
    message("Still unannotated: ", length(un), " protein groups, e.g. ",
            paste(head(un, 5), collapse = ", "))
  dat
}

#' Tag rows Human/Yeast from ProteinName (falls back to `prot`). Groups whose
#' name carries BOTH suffixes (mixed-species groups) are removed as ambiguous.
DualTP_assign_species <- function(dat) {
  src <- if ("ProteinName" %in% colnames(dat)) dat$ProteinName else dat$prot
  is_h <- grepl("HUMAN", src); is_y <- grepl("YEAST", src)
  dat$species <- ifelse(is_h & !is_y, "Human",
                        ifelse(is_y & !is_h, "Yeast", "remove"))
  tab <- table(dat$species)
  message("species assignment: ", paste(names(tab), format(as.vector(tab), big.mark = ","),
                                        sep = "=", collapse = ", "))
  dat <- dat %>% dplyr::filter(species != "remove")
  if (nrow(dat) == 0)
    stop("No Human/Yeast rows after species assignment; values looked like: ",
         paste(head(unique(src), 5), collapse = ", "))
  dat
}

