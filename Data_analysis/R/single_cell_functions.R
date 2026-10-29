################################################################################
# Fig. 5e-h: single cells and 10-cell bulk standards (9x3 plexDIA & timePlex)
################################################################################


## ---- Quality control and coverage (Fig. 5e-f) -------------------------------

# Cell types in ggplot's default hues, negative controls grey (Fig. 5e, f, h)
sc_cell_cols <- function(){
  c(stats::setNames(scales::hue_pal()(3), c("HEK293", "K562", "Raji")), "Negative control" = "grey60")
}

# Summed MS1 intensity (plex_Area_rtnorm, 1% FDR) per cell, and rel_MS1: its log10
# relative to the median of the other cells in the same run and timePlex channel
# (the 9 labels that went through one LC column together). This removes run and
# column differences, if any exist.
sc_summed_ms1 <- function(dat){
  d <- dat %>% dplyr::filter(Qvalue < 0.01) %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarize(summed_MS1 = sum(plex_Area_rtnorm),
                     Sample = dplyr::first(Sample)) %>% ungroup()
  set  <- sub("_[^_]*$", "", d$run_chan)   # run_chan = <Run>_<timePlex channel>_<label>
  cell <- !grepl("Neg", d$Sample)
  d$rel_MS1 <- sapply(seq_len(nrow(d)), function(i){
    others <- set == set[i] & cell & seq_len(nrow(d)) != i
    log10(d$summed_MS1[i]) - median(log10(d$summed_MS1[others]))
  })
  d
}

# Threshold from the negative controls: mean + number_sd SD of their rel_MS1
sc_neg_threshold <- function(d, number_sd = 2){
  neg <- d %>% dplyr::filter(grepl("Neg", Sample))
  if(nrow(neg) < 2){
    stop("Fewer than 2 Neg samples found - cannot compute mean + sd threshold.")
  }
  neg_thresh <- mean(neg$rel_MS1, na.rm = TRUE) + number_sd*sd(neg$rel_MS1, na.rm = TRUE)
  message(sprintf("Neg-based threshold: %.3f (log10 relative summed MS1; n = %d Neg samples)", neg_thresh, nrow(neg)))
  neg_thresh
}

# Cells above the negative-control threshold (Neg samples dropped), for Fig. 5f-h
sc_quality_control <- function(dat, number_sd = 2){
  d <- sc_summed_ms1(dat)
  neg_thresh <- sc_neg_threshold(d, number_sd)
  d %>% dplyr::filter(rel_MS1 >= neg_thresh, !grepl("Neg", Sample))
}

# Fig. 5e: relative summed MS1 per cell with the negative-control threshold
plot_sc_qc <- function(dat, number_sd = 2, out_file = "Figs/Fig5/Fig5e.pdf"){
  paper_figure(out_file)
  d <- sc_summed_ms1(dat)
  neg_thresh <- sc_neg_threshold(d, number_sd)
  n_keep  <- sum(d$rel_MS1 >= neg_thresh & !grepl("Neg", d$Sample))
  n_total <- sum(!grepl("Neg", d$Sample))

  d$Sample <- factor(ifelse(grepl("Neg", d$Sample), "Negative control", d$Sample),
                     levels = c("HEK293", "K562", "Raji", "Negative control"))

  p <- ggplot(d, aes(x = rel_MS1, fill = Sample)) +
    geom_density(alpha = 0.5, colour = "black", linewidth = 0.25) +
    geom_vline(xintercept = neg_thresh, linetype = "dashed", linewidth = 0.4) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, lineheight = 0.9, size = paper_geom_text(),
             label = sprintf("%d single cells\nretained\n(%.1f%%)", n_keep, 100 * n_keep / n_total)) +
    scale_fill_manual(values = sc_cell_cols(), name = NULL, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
    labs(x = "Log<sub>10</sub>, summed MS1 int<br>relative to median cell", y = "Density") +
    guides(fill = guide_legend(nrow = 2, byrow = FALSE,
                               override.aes = list(colour = "black", linewidth = 0.25))) +
    paper_theme() +
    theme(legend.position = "top", legend.justification = "left",
          legend.key.size = unit(7, "pt"), legend.key.spacing.y = unit(1.5, "pt"),
          axis.title.x = ggtext::element_markdown(lineheight = 1.1))
  
  invisible(list(fig5e = paper_panel(p, out_file, width = 1.56, height = 2.0)))
}

plot_counts <- function(sc, out_file = "Figs/Fig5/Fig5f.pdf"){
  paper_figure(out_file)

  prec <- sc  %>%
    dplyr::select("run_chan","prot","Sample") %>%
    dplyr::add_count(run_chan) %>% distinct(run_chan,.keep_all=T)
  prec$type <- "precursors"
  
  prot <- sc %>% dplyr::filter(Protein_Qvalue<0.01) %>%
    dplyr::distinct(run_chan,prot,.keep_all=T) %>%
    dplyr::select("run_chan","prot","Sample") %>% 
    dplyr::add_count(run_chan) %>% distinct(run_chan,.keep_all=T)
  prot$type <- "proteins"
  
  all <- rbind.fill(prot,prec)
  
  sum <- all %>% dplyr::group_by(type,Sample) %>% dplyr::summarize(med_n=median(n,na.rm=T))
  sum <- all %>% dplyr::group_by(type) %>% dplyr::summarize(max_n=max(n,na.rm=T))
  
  p <- ggplot(all %>% dplyr::filter(grepl("K562|Raji|HEK293",Sample) ), 
              aes(x=Sample,y=n,fill=Sample)) +
    geom_beeswarm(alpha = 0.5, size = 0.5, cex = 1.5) +
    geom_boxplot(alpha=0.75,outliers = F,width=0.5, linewidth = 0.25) + 
    facet_grid(~type) +
    scale_fill_manual(values = sc_cell_cols()) +
    scale_y_continuous(labels = scales::comma, limits = c(0, max(all$n))) +
    labs(x = NULL, y = "Count") +
    paper_theme_box() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 7))
  
  invisible(list(fig5f = paper_panel(p, out_file, width = 1.55, height = 2.0)))   
}


## ---- Join single cells with bulk, correct LC drift across runs --------------

pD_join_and_intersect <- function(data1, data2, intersect=T){
  dd1 <- unique(data1$seqcharge)
  dd2 <- unique(data2$seqcharge)
  keep <- dd1[dd1%in%dd2]
  all <- rbind.fill(data1,data2)
  if(intersect==T){
    all <- all %>% dplyr::filter(seqcharge%in%keep)
  }
  
  return(all)
}

pD_spline_correct_forLC <- function(SN_maxLFQprep_select){
  SN_maxLFQprep_select_1 <- SN_maxLFQprep_select %>% dplyr::mutate(log2val = log2(plex_Area_rtnorm)) %>% 
    dplyr::group_by(run_chan) %>% dplyr::mutate("med_log2val" = median(log2val[is.finite(log2val)])) %>%   # zeros (log2 = -Inf) must not set the median
    dplyr::mutate("log2val"=log2val-med_log2val) %>% ungroup() %>% dplyr::filter(is.finite(log2val))
  SN_maxLFQprep_select_1$run <- as.numeric(as.character(SN_maxLFQprep_select_1$RunOrder))
  
  results <- SN_maxLFQprep_select_1 %>%
    group_by(seqcharge, time_channel) %>%
    mutate(n = n(),
           n_runs = n_distinct(run)) %>%
    filter(n >= 50, n_runs >= 8) %>%
    do({
      fit <- smooth.spline(.$run, .$log2val)
      predicted <- predict(fit, .$run)
      residuals <- .$log2val - predicted$y
      RSS <- sum(residuals^2)
      TSS <- sum((.$log2val - mean(.$log2val))^2)
      R_squared <- 1 - (RSS / TSS)
      
      data.frame(run_chan = .$run_chan, time_channel = .$time_channel, run = .$run,
                 adjusted_log2val = residuals, log2val = .$log2val,
                 degrees = fit$df, R_squared = R_squared)
    }) %>%
    ungroup()
  results_keep <- results %>% filter(R_squared>0.1) %>%
    dplyr::select("run_chan","time_channel","seqcharge","adjusted_log2val")
  
  SN_maxLFQprep_select_1 <- SN_maxLFQprep_select_1 %>%
    left_join(results_keep, by = c("run_chan"="run_chan","time_channel"="time_channel","seqcharge"="seqcharge")) %>%
    dplyr::mutate("adj_val" = 2^(adjusted_log2val+med_log2val)) %>%
    dplyr::mutate("adj_val" = ifelse(is.na(adj_val)==F,adj_val,plex_Area_rtnorm))
  
  return(SN_maxLFQprep_select_1)
}


## ---- Protein quantification: normalization, imputation, batch correction ----

pD_prec_norm <- function(df, Quant = "plex_Area_rtnorm"){
  
  
  df$Quantnew <- log2(df[[Quant]])
  
  dat_fin <- df %>%
    # per-seqcharge mean BEFORE normalization, for restoring the scale later
    dplyr::group_by(seqcharge) %>%
    dplyr::mutate(prec_mean = mean(Quantnew, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # 1) sample loading: center each run_chan at its median
    dplyr::group_by(run_chan) %>%
    dplyr::mutate("norm_prec" = Quantnew - median(Quantnew, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # 2) within-run/column batch: center each seqcharge within Run x Column
    dplyr::group_by(seqcharge, Run, Column) %>%
    dplyr::mutate("norm_prec" = norm_prec - mean(norm_prec, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # 3) global: center each seqcharge across everything
    dplyr::group_by(seqcharge) %>%
    dplyr::mutate("norm_prec" = norm_prec - mean(norm_prec, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    # restore original scale: add back the pre-normalization seqcharge mean
    dplyr::mutate("norm_prec_log2" = norm_prec + prec_mean,
                  "norm_prec_lin"  = 2^norm_prec_log2)
  
  return(dat_fin)
}

pD_melt_MaxLFQ <- function(dat,meta_bulk_path){
  pD_1 <- dat
  print("new")
  pD_3_m <- melt(as.matrix(pD_1))
  m <- fread(meta_bulk_path)
  m$run_chan <- paste0(m$Run,"_",(m$Column-1),"_",m$Label)
  m$run_chan <- gsub('\\-', '.', m$run_chan)
  pD_3_m <- pD_3_m %>% left_join(m, by=c("Var2"="run_chan"))
  pD_3_m <- pD_3_m %>% dplyr::rename("prot"="Var1", "run_chan"="Var2")
  pD_3_m$norm_prot <- log2(as.numeric(pD_3_m$value))
  return(pD_3_m)
}

pD_prot_norm <- function(df, Quant = "prot_norm_adj"){

  df$Quantnew <- df[[Quant]]
  df$Genes <- df$prot
  dat_fin <- df %>% dplyr::group_by(run_chan) %>%
    dplyr::mutate("norm_prot" = Quantnew - median(Quantnew, na.rm=T)) %>%
    ungroup() %>% dplyr::group_by(Genes,Run,Column) %>%
    dplyr::mutate("norm_prot" = norm_prot - mean(norm_prot, na.rm=T)) %>% 
    ungroup() %>% dplyr::group_by(Genes) %>%
    dplyr::mutate("norm_prot" = norm_prot - mean(norm_prot, na.rm=T)) %>% ungroup()
  return(dat_fin)
}

pD_impute_log2 <- function(df, missing.prot.frac=0.7, missing.cell.frac=0.95, k=10){
  dat_fin1 <- df %>% dplyr::select("Genes", "run_chan", "norm_prot")
  dat_fin1$Genes <- as.factor(dat_fin1$Genes)
  dat_fin1$norm_prot <- as.numeric(dat_fin1$norm_prot)
  dat_fin1$run_chan <- gsub('\\.', '-', dat_fin1$run_chan)
  dat_fin1$run_chan <- as.factor(dat_fin1$run_chan)
  dat_fin2<-dcast(dat_fin1, Genes~run_chan, value.var = "norm_prot", fill=NA)
  dat3<-as.matrix(dat_fin2[,-1]); row.names(dat3)<-dat_fin2[,1]
  
  pre_impute <-(filt.mat.rc(dat3,missing.prot.frac,missing.cell.frac))
  sc.imp<-cr_norm_log(pre_impute)
  sc.imp <- hknn_weighted_EucNorm(sc.imp, k)
  sc.imp<-cr_norm_log((sc.imp))
  
  dat <- list(pre_impute = pre_impute, sc.imp = sc.imp)
  
  return(dat)
}

#from HS SCoPE2
filt.mat.rc<-function(mat, pct.r,pct.c){
  
  kr<-c()
  for(k in 1:nrow(mat)){
    
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    
  }
  
  mat<-mat[kr,]
  
  kc<-c()
  for(k in 1:ncol(mat)){
    
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    
  }
  
  mat<-mat[,kc]
  
  return(mat)
  
}

hknn_weighted_EucNorm <-function(dat, k){
  
  dat.imp<-dat
  
  #similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  counts <- pairwiseCount(dat)
  dist.mat <- dist.mat/counts
  
  cnames<-colnames(dist.mat)
  
  for(X in cnames){
    
    distances<-dist.mat[, X]
    
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    dat.reordered<-dat[ , names(distances.ordered ) ]
    vec<-dat[, X]
    na.index<-which( is.na(vec) )
    for(i in na.index){
      
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      distances.ordered_mat <- as.matrix(distances.ordered)
      
      distances.ordered_mat <- distances.ordered_mat[rownames(distances.ordered_mat)%in%closest.columns,] 
      closest_vec <- distances.ordered_mat[1:(k)]
      closest_vec_inv <- 1/closest_vec #inverted euclidian distances (small distances = more similar = should have greater weight)
      total <- sum(closest_vec_inv)
      
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the weighted average of the same row in k of the most similar columns
        vec[i]<-sum(closest_vec_inv*as.matrix(dat[ i, closest.columns[1:k] ] ))/total
        
      }
      
      
      if( length(closest.columns)<=k ){
        j <- length(closest.columns)
        closest_vec <- distances.ordered_mat[1:(j)]
        closest_vec_inv <- 1/closest_vec
        total <- sum(closest_vec_inv)
        vec[i]<-sum(closest_vec_inv*as.matrix(dat[ i, closest.columns[1:j] ] ))/total
        
      }
      
      
    }
    
    dat.imp[,X]<-vec
    
  }
  
  return(dat.imp)
  
}

cr_norm_log<-function(dat){
  
  for(k in 1:ncol(dat)){
    
    dat[,k]<-dat[,k]-median(dat[,k], na.rm = T)
    
    
  }
  
  
  for(k in 1:nrow(dat)){
    
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
    
  }
  
  return(dat)
}

pD_batchCorrect_labs_LC <- function(dat, meta_fpath){
  sc.imp <- as.matrix(dat$sc.imp)
  sc.imp <- cr_norm_log(sc.imp)
  m <- fread(meta_fpath)
  m$BioRep <- "1"
  m$run_chan <- paste0(m$Run,"_",(m$Column-1),"_",m$Label)
  m$Celltype <- m$Sample
  m$LC_batch <- m$Column-1
  no_neg <- m[!grepl("neg|Carrier|Neg",m$Celltype),]
  uni_BR <- unique(no_neg$BioRep)
  temp_list <- list()

  # Correction sequence: LC_batch, then Label
  cov_sequence <- c("LC_batch","Label")

  for(i in 1:length(uni_BR)){
    no_neg_temp <- no_neg[no_neg$BioRep==uni_BR[i],]
    temp <- sc.imp[,colnames(sc.imp)%in%no_neg_temp$run_chan]

    for(cv in cov_sequence){

      batch.covs <- as.character(no_neg_temp[[cv]][match(colnames(temp), no_neg_temp$run_chan)])

      # Filter groups with >= 3 cells for this covariate
      batch_counts <- table(batch.covs)
      valid_batches <- names(batch_counts[batch_counts >= 3])
      keep <- batch.covs %in% valid_batches
      temp <- temp[, keep]
      no_neg_temp <- no_neg_temp[no_neg_temp$run_chan %in% colnames(temp),]

      temp <- cr_norm_log(temp)
      batch.covs <- as.character(no_neg_temp[[cv]][match(colnames(temp), no_neg_temp$run_chan)])

      # Biological covariate to preserve: Sample
      sample.covs <- as.character(no_neg_temp$Sample[match(colnames(temp), no_neg_temp$run_chan)])

      # Only usable if Sample varies AND is not confounded with the batch
      # (confounded = every batch contains only one Sample level)
      use_mod <- length(unique(sample.covs)) > 1 &&
        !all(rowSums(table(batch.covs, sample.covs) > 0) <= 1) &&
        qr(model.matrix(~ sample.covs))$rank == ncol(model.matrix(~ sample.covs))

      # Skip ComBat if only one batch level remains
      if(length(unique(batch.covs)) > 1){
        if(use_mod){
          mod <- model.matrix(~ sample.covs)
          temp <- ComBat(temp, batch = batch.covs, mod = mod)
        } else {
          temp <- ComBat(temp, batch = batch.covs)
        }
      }
      temp <- cr_norm_log(temp)
    }

    temp_list[[i]] <- temp
  }

  matrix.sc.batch <- do.call("cbind", temp_list)
  imputed.BC <- cr_norm_log(matrix.sc.batch)
  imputed.BC_updated <- imputed.BC

  mat.sc.imp_NA <- imputed.BC_updated
  unimp <- dat$pre_impute
  unimp <- unimp[,colnames(unimp)%in%no_neg$run_chan]
  unimp <- unimp[,colnames(unimp)%in%colnames(matrix.sc.batch)]
  mat.sc.imp_NA[is.na(unimp)==T] <- NA
  mat.sc.imp_NA <- cr_norm_log(mat.sc.imp_NA)
  unimputed.BC <- melt(mat.sc.imp_NA); colnames(unimputed.BC)<-c("Genes","run_chan","norm_prot")

  dat_fin <- list(unimputed.BC = unimputed.BC, imputed.BC_updated = imputed.BC_updated, no_neg=no_neg)
  return(dat_fin)
}


## ---- Single cell vs bulk fold changes (Fig. 5g) and PCA (Fig. 5h) -----------

cor_bulk_SC <- function(BC_MLFQ, meta_path,min_frac_quant=0.5, out_file = "Figs/Fig5/Fig5g.pdf"){
  paper_figure(out_file)
  
  dat <- BC_MLFQ$unimputed.BC
  dat$norm_prot <- 2^dat$norm_prot
  m <- fread(meta_path)
  m$run_chan <- paste0(m$Run,"_",(m$Column-1),"_",m$Label)
  dat <- dat %>% left_join(m,by =c("run_chan" = "run_chan"))
  dat_consensus <- dat %>% 
    dplyr::group_by(Sample, Amount, Genes) %>% dplyr::add_count(name="total_possible") %>% dplyr::select(-c("Rep")) %>%
    na.omit() %>% dplyr::add_count(name="total_quantified") %>% dplyr::mutate("frac_quant" = total_quantified/total_possible) %>%
    dplyr::filter(frac_quant>min_frac_quant) %>%
    dplyr::summarize("mean_prot" = 2^mean(log2(norm_prot), na.rm = TRUE))   # mean of the log2 values per cell type
  dat_consensus_d <- dcast(dat_consensus, Genes~Sample+Amount, value.var = "mean_prot")
  dat_consensus_d$bulk_KU <- log2(dat_consensus_d$K562_bulk)-log2(dat_consensus_d$HEK293_bulk)
  dat_consensus_d$SC_KU <- log2(dat_consensus_d$K562_SC)-log2(dat_consensus_d$HEK293_SC)
  dat_consensus_d <- dat_consensus_d %>% na.omit()
  fit <- tls_fit(dat_consensus_d$bulk_KU, dat_consensus_d$SC_KU)
  
  r <- cor(dat_consensus_d$bulk_KU, dat_consensus_d$SC_KU, method = "pearson")
  p <- ggplot(dat_consensus_d, aes(x = bulk_KU, y = SC_KU)) +
    geom_point(alpha = 1, shape = 21, size = 1.3, color = "black", stroke = 0.25) +
    geom_pointdensity(alpha = 0.5, size = 1.1, stroke = 0) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = paper_geom_text(), parse = TRUE,
             label = sprintf('paste(italic(n), " = %s proteins")', scales::comma(nrow(dat_consensus_d)))) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 2.9, size = paper_geom_text(),
             label = sprintf("r = %.2f", r)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = "10x Bulk<br>Log<sub>2</sub>, K562 / HEK293",
         y = "Averaged single-cells<br>Log<sub>2</sub>, K562 / HEK293") +
    scale_color_viridis() +
    paper_theme() +
    theme(legend.position = "none",
          axis.title.x = ggtext::element_markdown(lineheight = 1.1),
          axis.title.y = ggtext::element_markdown(lineheight = 1.1))
  
  invisible(list(fig5g = paper_panel(p, out_file, width = 1.79, height = 2.0), data = dat_consensus_d))
}

tls_fit <- function(x, y) {
  fit <- prcomp(~ x + y)
  intercept <- fit$center[2] - fit$rotation[2,1] / fit$rotation[1,1] * fit$center[1]
  slope <- fit$rotation[2,1] / fit$rotation[1,1]
  list(intercept = intercept, slope = slope)
}

pD_weighted_PCA_timePlex <- function(dat, min_sc_frac = 0, out_file = "Figs/Fig5/Fig5h.pdf"){
  paper_figure(out_file)

  # only proteins quantified (not imputed) in at least min_sc_frac of the single cells
  obs <- dat$unimputed.BC %>% dplyr::filter(!grepl("ulk", run_chan)) %>%
    dplyr::group_by(Genes) %>% dplyr::summarise(frac = mean(!is.na(norm_prot)))
  mat.sc.imp <- dat$imputed.BC_updated
  mat.sc.imp <- mat.sc.imp[rownames(mat.sc.imp) %in% obs$Genes[obs$frac >= min_sc_frac], ]
  no_neg <- dat$no_neg
  r1<-cor(t(mat.sc.imp))
  rsum<-rowSums(r1^2)
  X.m <- mat.sc.imp
  X.m <- diag(rsum) %*%  X.m
  
  pca.imp.cor <- cor(X.m, use="complete.obs")
  sc.pca<-eigen(pca.imp.cor)
  scx<-as.data.frame(sc.pca$vectors)
  colnames(scx)<-paste0("PC",1:ncol(scx))
  scx$cells<-colnames(pca.imp.cor)
  pca_var <- sc.pca$values
  percent_var<- pca_var/sum(pca_var)*100
  pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
  pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)
  pca.display$celltype<-no_neg$Celltype[match(pca.display$id, no_neg$run_chan)]
  pca.display$sample<-no_neg$Amount[match(pca.display$id, no_neg$run_chan)]
  pca.display$celltype <- factor(pca.display$celltype)
  if (mean(pca.display$PC1[pca.display$celltype == "HEK293"], na.rm = TRUE) > 0) pca.display$PC1 <- -pca.display$PC1
  if (mean(pca.display$PC2[pca.display$celltype == "K562"], na.rm = TRUE) < 0) pca.display$PC2 <- -pca.display$PC2

  PCx<-"PC1"
  PCy<-"PC2"

  pca.display$sample_celltype <- paste0(pca.display$celltype,"_",pca.display$sample)
  bcount <- pca.display[grepl("ulk",pca.display$sample_celltype),]
  
  # Fig. 5h: single cells (circles) drawn over the 10-cell bulk samples (diamonds)
  pca <- ggplot(pca.display, aes(x = PC1, y = PC2, fill = celltype)) +
    geom_point(data = pca.display %>% dplyr::filter(grepl("ulk",sample_celltype)),
               aes(shape = "Bulk"), alpha = 0.6, colour = "black", size = 2.4, stroke = 0.25) +
    geom_point(data = pca.display %>% dplyr::filter(!grepl("ulk",sample_celltype)),
               aes(shape = "Single cell"), alpha = 0.6, colour = "black", size = 1.275, stroke = 0.25) +   # small enough that the bulk diamonds show
    scale_shape_manual(values = c("Single cell" = 21, "Bulk" = 23), breaks = c("Single cell", "Bulk"),
                       name = NULL) +
    scale_fill_manual(values = sc_cell_cols(), name = NULL) +
    annotate("text", x = -Inf, y = -Inf, hjust = -0.08, vjust = -0.3, lineheight = 0.9,
             size = paper_geom_text(),
             label = paste0((nrow(pca.display)-nrow(bcount)), " single cells\n",
                            scales::comma(dim(mat.sc.imp)[1]), " proteins")) +
    labs(x = paste0(PCx,"  (", round(percent_var[1],0),"%)"), y = paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
    guides(fill  = guide_legend(order = 1, override.aes = list(shape = 22, size = 2.3, alpha = 0.6)),
           shape = guide_legend(order = 2, override.aes = list(fill = "white", size = 2, alpha = 1))) +
    paper_theme() +
    theme(legend.position = "top", legend.justification = "left", legend.direction = "vertical",
          legend.key.spacing.y = unit(0, "pt"))
  
  invisible(list(fig5h = paper_panel(pca, out_file, width = 1.93, height = 2.0)))
}


## ---- Data tables for sharing (Gatto et al. 2023, Nat. Methods) ------

# Writes the Fig. 5 data as tab-separated tables with one column per sample, named
# by SampleID (= run_chan: <raw file>_<timePlex channel 0-2>_<PSMtag label>):
#   proteins x samples, log2: unimputed and imputed, both batch corrected (the
#     unimputed table is the imputed one with the imputed values set back to NA)
#   precursors x samples: MS1 plex_Area_rtnorm of every measured sample, and the
#     normalized log2 values that went into maxLFQ
#   sample metadata: one row per measured sample, negative controls and cells
#     that failed QC included
# README.md in out_dir describes the files and their columns.
sc_export_tables <- function(sc_dat2, sc_filt, joined_corrected, maxlfq, BC_MLFQ, meta_path, acq_path,
                             out_dir = "Data_tables/Fig5_single_cell", prefix = "timePlex_SC_"){
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  prec <- dplyr::bind_rows(sc_dat2$Amount_SC, sc_dat2$Amount_bulk)
  imp <- BC_MLFQ$imputed.BC_updated

  ## sample metadata
  ms1 <- sc_summed_ms1(sc_dat2$Amount_SC)
  summed <- prec %>% dplyr::filter(Qvalue < 0.01) %>% dplyr::group_by(run_chan) %>%
    dplyr::summarise(SummedMS1 = sum(plex_Area_rtnorm), .groups = "drop")
  counts <- prec %>% dplyr::group_by(run_chan) %>%   # counted as in Fig. 5f
    dplyr::summarise(nPrecursors = dplyr::n(),
                     nProteins = dplyr::n_distinct(prot[Protein_Qvalue < 0.01]), .groups = "drop")
  cell_lines <- data.frame(CellLine = c("HEK293", "K562", "Raji"),
                           CellLineSource = c("ATCC CRL-1573", "ATCC CCL-243", "ATCC CCL-86"))

  meta <- data.table::fread(meta_path) %>%
    dplyr::mutate(SampleID = paste0(Run, "_", Column - 1, "_", Label),
                  SampleType = dplyr::case_when(Sample == "Neg" ~ "negative control",
                                                Amount == "bulk" ~ "10-cell bulk",
                                                TRUE ~ "single cell"),
                  CellLine = ifelse(Sample == "Neg", NA, Sample)) %>%
    dplyr::left_join(cell_lines, by = "CellLine") %>%
    dplyr::left_join(data.table::fread(acq_path), by = "Run") %>%
    dplyr::left_join(summed, by = c("SampleID" = "run_chan")) %>%
    dplyr::left_join(dplyr::select(ms1, run_chan, RelMS1 = rel_MS1), by = c("SampleID" = "run_chan")) %>%
    dplyr::left_join(counts, by = c("SampleID" = "run_chan")) %>%
    dplyr::mutate(Organism = ifelse(is.na(CellLine), NA, "Homo sapiens"),
                  Condition = ifelse(is.na(CellLine), NA, "untreated"),
                  SamplePrep = ifelse(Amount == "bulk", "mPOP lysate, 10 cells per channel", "nPOP"),
                  ChemicalLabel = paste0("PSMtag-d", Label),
                  PassedQC = ifelse(SampleType == "single cell", SampleID %in% sc_filt$run_chan, NA),
                  InProteinTables = SampleID %in% colnames(imp))
  used <- function(sc, qc, tab) paste(c(if (sc) "5e", if (qc) "5f", if (tab) c("5g", "5h")), collapse = ";")
  meta$UsedInFigure <- mapply(used, meta$Amount == "SC", meta$PassedQC %in% TRUE, meta$InProteinTables)
  meta <- meta %>% dplyr::arrange(RunOrder, Column, Label) %>%
    dplyr::select(SampleID, SampleType, CellLine, CellLineSource, Organism, Condition, SamplePrep,
                  PooledSet = V9, RawFile = Run, mzMLFile, AcquisitionStart, RunOrder, Instrument,
                  TimePlexChannel = Column, ChemicalLabel, SummedMS1, RelMS1, PassedQC,
                  nPrecursors, nProteins, InProteinTables, UsedInFigure)
  stopifnot(!anyDuplicated(meta$SampleID), all(colnames(imp) %in% meta$SampleID))
  ids <- meta$SampleID[meta$InProteinTables]

  ## proteins x samples
  prot_names <- dplyr::distinct(prec, prot, ProteinName)
  stopifnot(!anyDuplicated(prot_names$prot))
  u <- BC_MLFQ$unimputed.BC
  unimp <- matrix(NA_real_, nrow(imp), ncol(imp), dimnames = dimnames(imp))
  unimp[cbind(as.character(u$Genes), as.character(u$run_chan))] <- u$norm_prot
  prot_table <- function(mat){
    mat <- mat[order(rownames(mat)), ids]
    data.frame(Protein = rownames(mat), ProteinName = prot_names$ProteinName[match(rownames(mat), prot_names$prot)],
               mat, check.names = FALSE)
  }

  ## precursors x samples
  # the maxLFQ input, filtered as in the SC_bulk_maxLFQ target
  mlfq_in <- joined_corrected %>%
    dplyr::filter(!grepl("Neg", Sample) & plex_Area_rtnorm > 1 & Protein_Qvalue < 0.01 &
                    (run_chan %in% sc_filt$run_chan | grepl("ulk", run_chan))) %>%
    pD_prec_norm(Quant = "plex_Area_rtnorm")
  stopifnot(setequal(mlfq_in$prot, rownames(maxlfq)), setequal(mlfq_in$run_chan, colnames(maxlfq)))
  prec_table <- function(d, value, samples){
    w <- as.data.frame(data.table::dcast(data.table::as.data.table(d), prot + seqcharge ~ run_chan, value.var = value))
    data.frame(Precursor = w$seqcharge, Protein = w$prot,
               ProteinName = prot_names$ProteinName[match(w$prot, prot_names$prot)],
               w[, samples], check.names = FALSE)
  }

  tables <- list(metadata = meta, proteins_unimputed = prot_table(unimp), proteins_imputed = prot_table(imp),
                 precursors_MS1area = prec_table(prec, "plex_Area_rtnorm", meta$SampleID),
                 precursors_normalized = prec_table(mlfq_in, "norm_prec_log2", ids))
  files <- file.path(out_dir, paste0(prefix, names(tables), ".tsv"))
  for (i in seq_along(tables)) data.table::fwrite(tables[[i]], files[i], sep = "\t", na = "NA", quote = FALSE)
  files
}

