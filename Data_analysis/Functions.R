
pD_cost <- function(d0=778.64, d4=2348.04, d8=2343.62, numcolumns=3,pack10columns=10160,runspercolumn=2500){
  #Figure 1d:
  
  total_prot <- 50*3*0.1 #50 aliquots * 3 labels * 0.1 mg protein/aliquot
  
  total_cost <- d0+d4+d8
  
  cost_per_ug <- total_cost/(total_prot*1000)
  cost_extra_columns <- ((pack10columns/10)*(numcolumns-1))/runspercolumn #b/c using 3 extra column compared to non-multiplexed.. assume 1500 runs per column
  plexDIA_timePlex <- (100*1.5)/12
  timePlex <- (100*1.5)/numcolumns
  plexDIA <- (100*1.5)/3
  LFDIA <- (100*1.5)
  cost <- c(LFDIA, plexDIA, (cost_per_ug), timePlex, cost_extra_columns/(numcolumns-1), plexDIA_timePlex, (cost_per_ug), cost_extra_columns/3)
  Application_Reagent <- c("LC-MS/MS cost","LC-MS/MS cost","Label cost" ,"LC-MS/MS cost", "Extra columns cost","LC-MS/MS cost","Label cost","Extra columns cost")
  type <- c("LF-DIA", "3-plexDIA","3-plexDIA","3-timeplex","3-timeplex", "3-timeplex & 3-plexDIA","3-timeplex & 3-plexDIA","3-timeplex & 3-plexDIA")
  
  df <- data.frame(cost = cost, Application_Reagent = Application_Reagent,type = type)
  #df$error <- c(25,75)
  
  df$Application_Reagent <- factor(df$Application_Reagent, levels=c("Extra columns cost", "Label cost", "LC-MS/MS cost"))
  df$type <- factor(df$type, levels=c("LF-DIA", "3-plexDIA", "3-timeplex", "3-timeplex & 3-plexDIA"))
  
  total_costs_vec <- c(LFDIA, plexDIA + (cost_per_ug/3)*1.5, timePlex + cost_extra_columns/(numcolumns-1), plexDIA_timePlex + (cost_per_ug/3)*1.5 + cost_extra_columns/(numcolumns-1))
  type_unique <- c("LF-DIA", "3-plexDIA", "3-timeplex", "3-timeplex & 3-plexDIA")  
  df_unique <- data.frame(totalcosts = total_costs_vec, type=type_unique)
  
  ggplot(df, aes(x=type, y=cost, fill=type)) + 
    theme_classic() + labs(x="", y="Estimated cost (USD) per sample",subtitle = "Assuming 1 ug protein per sample") +
    geom_bar(data = df %>% dplyr::filter(type == "LF-DIA"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == "3-plexDIA"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == "3-timeplex"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == "3-timeplex & 3-plexDIA"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    scale_y_continuous(labels=scales::dollar_format()) + 
    scale_fill_manual(values = c("grey39","red","blue")) +
    theme(legend.position = c(0.7,0.8),
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          axis.text.x = element_text(size=13, face = "bold", color = "black"),
          axis.text.y = element_text(size=14, face = "bold", color = "black"),
          axis.title.x = element_text(size=12, color = "black")) + coord_flip()# +
  #geom_text(data=df_unique, aes(x=type, y=totalcosts+3, label=totalcosts))
  
  ggsave("Figs/Figure1_cost.pdf", width=6.3, height=1.9, dpi=500)
}

pD_cost_plexDIA <- function(d0=778.64, allothers=2348.04,midplex=10,highplex=100){
  #Figure 1d:
  
  ## 3-plexDIA
  total_prot <- 50*3*0.1 #50 aliquots * 3 labels * 0.1 mg protein/aliquot
  total_cost <- d0+(3-1)*allothers
  
  plex_cost_per_ug_total <- (total_cost/(total_prot*1000))#*3
  
  #################
  # 10-plexDIA
  midplex_total_prot <- 50*midplex*0.1 #50 aliquots * 3 labels * 0.1 mg protein/aliquot
  midplex_total_cost <- d0+(midplex-1)*allothers
  
  midplex_cost_per_ug_total <- (midplex_total_cost/(midplex_total_prot*1000))#*midplex
  
  #################
  # 100-plexDIA
  highplex_total_prot <- 50*highplex*0.1 #50 aliquots * 3 labels * 0.1 mg protein/aliquot
  highplex_total_cost <- d0+(highplex-1)*allothers
  
  highplex_cost_per_ug_total <- (highplex_total_cost/(highplex_total_prot*1000))#*highplex
  
  # cost_extra_columns <- ((pack10columns/10)*3)/runspercolumn #b/c using 3 extra column compared to non-multiplexed.. assume 1500 runs per column
  LFDIA <- (100*1.5)
  cost <- c(LFDIA, LFDIA/3, plex_cost_per_ug_total, LFDIA/midplex, midplex_cost_per_ug_total, LFDIA/highplex,highplex_cost_per_ug_total)
  Application_Reagent <- c("LC-MS/MS cost","LC-MS/MS cost","Label cost" ,"LC-MS/MS cost","Label cost","LC-MS/MS cost","Label cost")
  type <- c("LF-DIA", "3-plexDIA","3-plexDIA",paste0(midplex,"-plexDIA"),paste0(midplex,"-plexDIA"), 
            paste0(highplex,"-plexDIA"),paste0(highplex,"-plexDIA"))
  
  df <- data.frame(cost = cost, Application_Reagent = Application_Reagent,type = type)
  #df$error <- c(25,75)
  
  df$Application_Reagent <- factor(df$Application_Reagent, levels=c("Extra columns cost", "Label cost", "LC-MS/MS cost"))
  df$type <- factor(df$type, levels=c("LF-DIA", "3-plexDIA", paste0(midplex,"-plexDIA"), paste0(highplex,"-plexDIA")))
  
  total_costs_vec <- c(LFDIA, LFDIA/3 + plex_cost_per_ug_total, LFDIA/midplex + midplex_cost_per_ug_total, LFDIA/highplex + highplex_cost_per_ug_total)
  type_unique <- c("LF-DIA", "3-plexDIA", paste0(midplex,"-plexDIA"), paste0(highplex,"-plexDIA"))
  df_unique <- data.frame(totalcosts = total_costs_vec, type=type_unique)
  
  ggplot(df, aes(x=type, y=cost, fill=type)) + 
    theme_classic() + labs(x="", y="Estimated cost (USD) per sample",subtitle = "Assuming 1 ug protein per sample") +
    geom_bar(data = df %>% dplyr::filter(type == "LF-DIA"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == "3-plexDIA"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(midplex,"-plexDIA")),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(highplex,"-plexDIA")),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    scale_y_continuous(labels=scales::dollar_format()) + 
    scale_fill_manual(values = c("grey39","red","green")) +
    theme(legend.position = c(0.7,0.8),
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          axis.text.x = element_text(size=13, face = "bold", color = "black"),
          axis.text.y = element_text(size=14, face = "bold", color = "black"),
          axis.title.x = element_text(size=12, color = "black")) + coord_flip()# +
  #geom_text(data=df_unique, aes(x=type, y=totalcosts+3, label=totalcosts))
  
  
  cost <- c(LFDIA, LFDIA/3, plex_cost_per_ug_total/10, LFDIA/midplex, midplex_cost_per_ug_total/10, LFDIA/highplex,highplex_cost_per_ug_total/10)
  Application_Reagent <- c("LC-MS/MS cost","LC-MS/MS cost","Label cost" ,"LC-MS/MS cost","Label cost","LC-MS/MS cost","Label cost")
  type <- c("LF-DIA", "3-plexDIA","3-plexDIA",paste0(midplex,"-plexDIA"),paste0(midplex,"-plexDIA"), 
            paste0(highplex,"-plexDIA"),paste0(highplex,"-plexDIA"))
  
  df <- data.frame(cost = cost, Application_Reagent = Application_Reagent,type = type)
  #df$error <- c(25,75)
  
  df$Application_Reagent <- factor(df$Application_Reagent, levels=c("Extra columns cost", "Label cost", "LC-MS/MS cost"))
  df$type <- factor(df$type, levels=c("LF-DIA", "3-plexDIA", paste0(midplex,"-plexDIA"), paste0(highplex,"-plexDIA")))
  
  
  ggplot(df, aes(x=type, y=cost, fill=type)) + 
    theme_classic() + labs(x="", y="Estimated cost (USD) per sample",subtitle = "Assuming 100 ng protein per sample") +
    geom_bar(data = df %>% dplyr::filter(type == "LF-DIA"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == "3-plexDIA"),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(midplex,"-plexDIA")),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(highplex,"-plexDIA")),#colour="black", 
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    scale_y_continuous(labels=scales::dollar_format()) + 
    scale_fill_manual(values = c("grey39","red","blue")) +
    theme(legend.position = c(0.7,0.8),
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          axis.text.x = element_text(size=13, face = "bold", color = "black"),
          axis.text.y = element_text(size=14, face = "bold", color = "black"),
          axis.title.x = element_text(size=12, color = "black")) + coord_flip()# +
  
  ggsave("Figs/Figure1_cost_plexDIAOnly.pdf", width=4.1, height=1.9, dpi=500)
}

Jmod_read_folder <- function(folder_path, lib_fpath,distinct=F,filtered=T){
  print("starting")
  file_list <- list.files(folder_path)
  data_list <- list()
  counter = 0
  for (file_name in file_list) {
    counter = counter + 1 
    if(filtered==T){
      file_path <- file.path(folder_path, file_name,'filtered_IDs.csv')
    } else{
      file_path <- file.path(folder_path, file_name,'all_IDs.csv')
      
    }
    data <- data.frame(fread(file_path, 
                             select = c("seq","z","time_channel","rt","lib_rt","stripped_seq","MS1_Int","untag_prec","channel","coeff", "protein",
                                        "plexfittrace_ps_all","plexfittrace_all", "MS1_Area", "all_ms1_iso0vals","mz",
                                        "Qvalue","Protein_Qvalue","BestChannel_Qvalue","plexfitMS1","plexfitMS1_p","decoy")))
    data <- data[data$decoy==F,] #remove decoys
    temp <- str_extract(file_name, "[^_]+")
    data$Run <- temp
    data_list[[counter]] <- data
    print(paste0("Loaded ",counter,"/",length(file_list)," files"))
  }
  df <- do.call("rbind.fill", data_list)
  if(distinct==T){
    lib <- fread(lib_fpath)
    lib$seqcharge <- paste0(lib$PeptideSequence,"_",lib$PrecursorCharge)
    lib <- lib %>% distinct(transition_group_id, ProteinName, .keep_all=T)
    df <- df %>% left_join(lib, by =c("untag_prec" = "seqcharge"))
  }
  df$MS1_Int[df$MS1_Int < 1] <- NA #replace values <1 with NA
  return(df)
}

pD_get_data_Jmod <- function(dat, meta_path,infer_sample){
  m <- fread(meta_path)
  dat <- dat %>% dplyr::rename("seqcharge" = "untag_prec")
  dat <- dat %>% dplyr::rename("prot" = "protein")
  if (!"time_channel" %in% colnames(dat)) {
    dat$time_channel <- NA
  }
  
  dat$channel[is.na(dat$time_channel)&(!grepl("mTRAQ",dat$seq))] <- "NoPlex_LF"
  dat$channel[is.na(dat$channel) & is.na(dat$time_channel) & grepl("mTRAQ", dat$seq)] <- ""
  
  dat$time_channel[is.na(dat$time_channel)] <- 0
  dat <- dat %>% dplyr::mutate("channel" = ifelse(grepl("_", channel),channel,
                                                  ifelse(grepl("mTRAQ",seq), paste0("_",channel),paste0(channel,"_LF"))))
  
  dat$run_chan <- paste0(dat$Run, "_", dat$channel)
  m$Column <- as.numeric(m$Column)-1
  m$Run <- sub("_.*", "", m$Run)
  
  m <- m %>% dplyr::mutate("TP" = ifelse(timePlex==T,Column,
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
  fin <- list()
  for(i in 1:length(unique_loadings)){
    temp <- dat[dat$Amount==unique_loadings[i],]
    dats <- paste0("Amount_",unique_loadings[i])
    fin[[dats]] <- temp
  }
  return(fin)
}

pD_plot_RTs_TP <- function(dat=benchmarking_dat8, libpath = "", run="JD0311"){
  dat_specific <- dat[dat$Run==run,] %>% dplyr::filter(Qvalue<0.2 & BestChannel_Qvalue<0.01 & decoy==F)
  lib <- fread(libpath) %>% dplyr::distinct(PeptideSequence,.keep_all=T) %>% dplyr::select(PeptideSequence,Tr_recalibrated) 
  dat_specific <- dat_specific %>% left_join(lib, by = c("stripped_seq"="PeptideSequence"))
  
  chan0 <- dat_specific[dat_specific$time_channel==1,] %>% dplyr::select("seqcharge","lib_rt") %>%
    dplyr::rename("lib_rt_col1" = "lib_rt")
  dat_specific_unified <- dat_specific %>% inner_join(chan0, by =c("seqcharge"="seqcharge"))
  ggplot(dat_specific_unified, aes(y=lib_rt,x=Tr_recalibrated)) + geom_point(alpha=0.5,shape=21,aes(fill=as.factor(time_channel)),size=1.5) + theme_bw()+
    scale_fill_manual(values=c("#16330E","#4D9934","#E1F7DA")) + labs(fill="timePlex", x="Original library RT, min", y="Empirical RT, min")
  
  
  dat_specific_unified <- dat_specific_unified[!is.na(dat_specific_unified$Tr_recalibrated),]
  
  dat_specific_unified <- dat_specific_unified %>%
    dplyr::group_by(time_channel) %>%
    dplyr::mutate(lowess_fit = lowess(Tr_recalibrated, rt, f = 0.1)$y)  # Adjust f for smoothing
  
  dat_specific_unified$d_rt_lowess <- dat_specific_unified$rt-dat_specific_unified$lowess_fit
  #sel <- t1 %>% dplyr::select(seqcharge,Tr_recalibrated,rt)
  # Plot
  ggplot(dat_specific_unified, aes(y = rt, x = Tr_recalibrated)) + 
    geom_point(alpha = 0.5, shape = 21, aes(fill = as.factor(time_channel)), size = 2.5) +
    theme_bw() +
    scale_fill_manual(values = c("#16330E", "#4D9934", "#E1F7DA")) +
    scale_color_manual(values = c("#16330E", "#4D9934", "#E1F7DA")) +  # Match colors for LOWESS
    labs(fill = "timePlex", color = "timePlex", 
         x = "Original library RT, min", y = "Empirical RT, min",
         title = "LF 3-timePlex, 8ng",subtitle=paste0(run))
  ggsave(paste0("RTs_original",run,".pdf"),width=5,height=5)
  
  t1 <- dat_specific_unified[dat_specific_unified$time_channel==0,]
  ggplot(dat_specific_unified, aes(y = rt, x = Tr_recalibrated)) + 
    # Plot all other points first
    geom_point(
      data = dat_specific_unified %>% dplyr::filter(seqcharge != "TRPVVAAGAVGLAQR_3"),
      alpha = 0.5, shape = 21, size = 2.5,
      aes(fill = as.factor(time_channel))
    ) +
    # Plot the red points on top
    geom_smooth(aes(group=as.factor(time_channel),color=as.factor(time_channel)),se=F)+
    geom_point(
      data = dat_specific_unified %>% dplyr::filter(seqcharge == "TRPVVAAGAVGLAQR_3"),
      alpha = 1, shape = 21, size = 3, fill = "#FF13F0", color = "black" # Black outline for visibility
    ) +
    theme_classic() +
    scale_fill_manual(values = c("#16330E", "#4D9934", "#E1F7DA")) + # Original colors
    scale_color_manual(values = c("#16330E", "#4D9934", "#E1F7DA")) +  
    labs(fill = "timePlex", color = "timePlex", 
         x = "Original library RT, min", y = "Empirical RT, min",
         title = "LF 3-timePlex, 8ng",subtitle=paste0(run," TRPVVAAGAVGLAQR_3")) + xlim(-50,100)
  
  ggsave(paste0("RTs_original_outliers",run,".pdf"),width=5,height=5)
  
  ggplot(dat_specific_unified, aes(y=rt,x=lib_rt_col1)) + geom_point(alpha=0.1,shape=21,aes(fill=as.factor(time_channel)),size=2) + theme_bw()+
    scale_fill_manual(values=c("#16330E", "#4D9934", "#E1F7DA")) + labs(fill="timePlex", x="Fine-tuned library RT of timePlex-1, min", y="Empirical RT, min",
                                                                        title = "LF 3-timePlex, 8ng",subtitle=paste0(run))
  ggsave(paste0("RTs_finetuned_",run,".pdf"),width=5,height=5)
  
  
  ggplot(dat_specific_unified, aes(y=rt,x=lib_rt_col1)) + geom_point(alpha=0.5,shape=21,
                                                                     aes(fill=as.factor(time_channel)),size=2.5) + 
    geom_point(
      data = dat_specific_unified %>% dplyr::filter(seqcharge == "TRPVVAAGAVGLAQR_3"),
      alpha = 1, shape = 21, size = 3, fill = "#FF13F0", color = "black" # Black outline for visibility
    ) + theme_classic()+
    scale_fill_manual(values=c("#16330E", "#4D9934", "#E1F7DA")) + 
    labs(fill="timePlex", x="Fine-tuned library RT of timePlex-1, min", y="Empirical RT, min",
         title = "LF 3-timePlex, 8ng",subtitle=paste0(run," TRPVVAAGAVGLAQR_3"))
  ggsave(paste0("RTs_finetuned_outliers_",run,".pdf"),width=5,height=5)
  
  dat_specific_unified$d_RT <- dat_specific_unified$rt-dat_specific_unified$lib_rt_col1
  chan1 <- dat_specific_unified[dat_specific_unified$time_channel==1,] %>% dplyr::select("seqcharge","rt") %>% dplyr::rename("rt_first"="rt")
  
  dat_specific_unified2 <- dat_specific_unified %>% inner_join(chan1, by =c("seqcharge"="seqcharge"))
  dat_specific_unified2$delta_RT_first <- dat_specific_unified2$rt - dat_specific_unified2$rt_first
  
  ggplot(dat_specific_unified2, aes(x=rt_first,y=delta_RT_first)) + geom_point(alpha=0.5,shape=21,aes(fill=as.factor(time_channel.x)),size=2.5) + theme_bw()+
    scale_fill_manual(values=c("#16330E", "#4D9934", "#E1F7DA")) + labs(fill="timePlex", y=paste0("\U0394RT from timePlex-1, min"), x="RT from timePlex-1, min",
                                                                        title = "LF 3-timePlex, 8ng",subtitle=paste0(run)) #+ ylim(0,60)
  ggsave(paste0("RTs_differencefromFirst_",run,".pdf"),width=5,height=5)
  
  return(dat_specific_unified2)
}

tp_coverage <- function(dat,extra="", use_best_channel = F, filter_level=0.01){
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(BestChannel_Qvalue<filter_level) %>% dplyr::filter(Qvalue < 0.2)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<filter_level)
  }
  
  dat <- dat %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))
  dat <- dat[dat$Repeatability ==F,]
  dat_keep <- dat
  
  if ("Implementation" %in% colnames(dat)) {
    dat$Type <- paste0(dat$Type,"_",dat$Implementation)
    dat <- dat %>% dplyr::filter(Type == "3-plexDIA" | Type=="LF\nNo-plex")
  }
  
  dat <- dat %>% dplyr::mutate(Label = ifelse(Type == "3-plexDIA","mTRAQ",Label))
  
  
  
  dat_counts <- dat %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    dplyr::ungroup()
  
  unique_combinations <- dat %>%
    distinct(Sample, Column, Amount, Label)
  
  missing_3plexDIA <- unique_combinations %>%
    dplyr::mutate(Type = "3-plexDIA", n = 0, Run = "Placeholder_Run")
  
  if (!"3-plexDIA" %in% dat$Type) {
    dat_counts <- bind_rows(dat_counts, missing_3plexDIA)
  }
  
  dat_counts$Sample <- factor(dat_counts$Sample, levels = c("A", "B", "C"))
  
  
  # Compute mean and standard error
  dat_summary <- dat_counts %>%
    dplyr::group_by(Sample, Label, Amount, Type) %>%
    dplyr::summarise(
      mean_count = mean(n, na.rm = TRUE),
      se = sd(n, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% dplyr::ungroup()
  
  dat_summary <- dat_summary %>%
    dplyr::arrange(desc(Sample), Type) %>%
    
    dplyr::group_by(Type) %>%
    dplyr::mutate(
      stacked_mean = cumsum(mean_count)) %>% dplyr::ungroup() %>%
    dplyr::group_by(Label, Sample, Type) %>%
    dplyr::mutate(stacked_mean = ifelse(Type=="LF\nNo-plex", mean_count,stacked_mean))%>%
    dplyr::mutate(
      er_max = stacked_mean + se,
      er_min = stacked_mean - se
    ) %>%
    dplyr::ungroup()
  
  df_text <- dat_summary %>%
    dplyr::group_by(Type, Sample,Label) %>%
    dplyr::summarise(
      text_value = round(max(er_max) - max(se)),
      text_pos = max(er_max) + 2000
    )
  
  my_colors <- c("#7ca982","#4e598c","#fed766")
  
  ggplot(dat_summary, aes(x = Type, y = mean_count, fill = Sample)) +
    scale_x_discrete(limits = c("LF\nNo-plex",  "3-plexDIA", "LF\n3-timePlex", "3-plexDIA\n3-timePlex")) +
    
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA\n3-timePlex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\n3-timePlex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 1.17),
             color = "black", width = 1, alpha = 1) +
    
    scale_y_continuous(labels = comma) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("3-plexDIA\n3-timePlex", "LF\n3-timePlex", "3-plexDIA")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.18, size = 0.4) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.5, size = 0.4,
                  position = position_dodge(1.17)) +
    
    scale_fill_manual(values = my_colors) +
    theme_classic() +
    
    labs(y = "Precursors data points / Run", x = "", fill = "Samples") +
    
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = c(0.16, 0.9),
          legend.direction = "horizontal",
          legend.spacing.x = unit(1.0, 'line')) +
    
    guides(fill = guide_legend(override.aes = list(alpha = 0.7),
                               title = "Samples",
                               label.position = "bottom",
                               title.position = "top", title.vjust = 1, title.hjust = 0.5)) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("3-plexDIA\n3-timePlex", "LF\n3-timePlex", "3-plexDIA")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3) +
    
    geom_text(data = df_text %>% dplyr::filter(Type == "LF\nNo-plex"), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3,
              width = 0.5, size = 0.4,
              position = position_dodge(1.17))
  
  ggsave(paste0(extra,"_","Prec_coverage_","BestChannel_",use_best_channel,".pdf"),width=6,height=5)
  
  
  
  ##################
  dat <- dat_keep %>% dplyr::mutate(Label = ifelse(Type == "3-plexDIA","mTRAQ",Label))
  
  dat_counts <- dat %>% dplyr::filter(Protein_Qvalue<filter_level) %>% dplyr::distinct(Sample, Run, Column, Amount, Label, Type, prot,.keep_all=T) %>% 
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>% dplyr::ungroup()
  
  dat_summary <- dat_counts %>%
    dplyr::group_by(Sample, Label, Amount, Type) %>%
    dplyr::summarise(
      mean_count = mean(n, na.rm = TRUE),
      se = sd(n, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% dplyr::ungroup()
  
  dat_summary <- dat_summary %>%
    dplyr::arrange(desc(Sample), Type) %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(
      stacked_mean = cumsum(mean_count)) %>% ungroup() %>%
    dplyr::group_by(Label, Sample, Type) %>%
    dplyr::mutate(stacked_mean = ifelse(Type=="LF\nNo-plex", mean_count,stacked_mean))%>%
    dplyr::mutate(
      er_max = stacked_mean + se,
      er_min = stacked_mean - se
    ) %>%
    dplyr::ungroup()
  
  df_text <- dat_summary %>%
    dplyr::group_by(Type, Sample,Label) %>%
    dplyr::summarise(
      text_value = round(max(er_max) - max(se)),
      text_pos = max(er_max) + 310
    )
  
  
  ggplot(dat_summary, aes(x = Type, y = mean_count, fill = Sample)) +
    scale_x_discrete(limits = c("LF\nNo-plex", "3-plexDIA", "LF\n3-timePlex", "3-plexDIA\n3-timePlex")) +
    
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA\n3-timePlex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\n3-timePlex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 1.17),
             color = "black", width = 1, alpha = 1) +
    
    scale_y_continuous(labels = comma) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("3-plexDIA\n3-timePlex", "LF\n3-timePlex", "3-plexDIA")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.18, size = 0.4) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.5, size = 0.4,
                  position = position_dodge(1.17)) +
    
    scale_fill_manual(values = my_colors) +
    theme_classic() +
    
    labs(y = "Proteins data points / Run", x = "", fill = "Samples") +
    
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = c(0.16, 0.9),
          legend.direction = "horizontal",
          legend.spacing.x = unit(1.0, 'line')) +
    
    guides(fill = guide_legend(override.aes = list(alpha = 0.7),
                               title = "Samples",
                               label.position = "bottom",
                               title.position = "top", title.vjust = 1, title.hjust = 0.5)) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("3-plexDIA\n3-timePlex", "LF\n3-timePlex", "3-plexDIA")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3) +
    
    geom_text(data = df_text %>% dplyr::filter(Type == "LF\nNo-plex"), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3,
              width = 0.5, size = 0.4,
              position = position_dodge(1.17))
  
  ggsave(paste0(extra,"_","Prot_coverage_","BestChannel_",use_best_channel,".pdf"),width=6,height=5)
  
  return(dat_counts)
  
  
}

rm_overlaping_sequences <- function(dat, overlap_path){
  overlaps <- fread(overlap_path)
  print(nrow(dat))
  dat <- dat[dat$stripped_seq%!in%overlaps$PeptideSequence,]
  print(nrow(dat))
  
  return(dat)
}

sum_fitted_values <- function(fitted_str) {
  if (is.na(fitted_str)) return(NA)
  values <- str_extract_all(fitted_str, "-?\\d+\\.\\d+(e-?\\d+)?")[[1]] %>% as.numeric()
  sum(values, na.rm = TRUE)
}

pD_ms1_quant <- function(df, adjacent_scans = 1,dontapplyLF=T) {
  df_lim <- df %>% dplyr::filter(grepl("mTRAQ", seq)) %>% dplyr::filter(BestChannel_Qvalue < 0.01 & Qvalue<0.2)
  df_lim2 <- df_lim %>%
    dplyr::mutate(
      MS1_fitted = pmap_chr(list(plexfittrace_ps_all, plexfittrace_all), function(ps, tr) {
        extract_highest_pearson_vec(ps, tr, adjacent_scans)
      }),
      MS1_fitted_sum = map_dbl(MS1_fitted, sum_fitted_values)
    )
  
  # for LF get the apex and adj scans.
  df_complement <- df %>% dplyr::filter(!(grepl("mTRAQ", seq) & BestChannel_Qvalue < 0.01  & Qvalue<0.2))
  if(dontapplyLF==F){
    df_complement <- df_complement %>%
      dplyr::mutate(
        MS1_fitted = pmap_chr(list(all_ms1_iso0vals, all_ms1_iso0vals), function(ps, tr) {
          extract_highest_pearson_vec(ps, tr, adjacent_scans)
        }),
        MS1_fitted_sum = map_dbl(MS1_fitted, sum_fitted_values)
      )
  }
  all <- bind_rows(df_complement, df_lim2) %>%
    dplyr::mutate(MS1_bestfit = ifelse(is.na(MS1_fitted_sum), plexfitMS1, MS1_fitted_sum))
  
  return(all)
}

extract_highest_pearson_vec <- function(pearson_str, trace_str, n) {
  if (is.na(pearson_str) || is.na(trace_str)) return(NA)  # Handle NA cases
  
  # Split the semicolon-separated numeric strings and convert to numeric
  pearson_stats <- str_split(pearson_str, ";")[[1]] %>% as.numeric()
  trace_entries <- str_split(trace_str, ";")[[1]]
  
  if (length(pearson_stats) == 0 || length(trace_entries) == 0) return(NA)  
  
  max_idx <- which.max(pearson_stats)
  start_idx <- max(1, max_idx - n)
  end_idx <- min(length(trace_entries), max_idx + n)
  
  paste(trace_entries[start_idx:end_idx], collapse = ";")
}

prec_normalizeQuant <- function(dat,quant_col="MS1_Int",use_best_channel = F, filter_level=0.01){
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(BestChannel_Qvalue<filter_level) %>% dplyr::filter(Qvalue < 0.2)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<filter_level)
  }
  if (!"MS1_bestfit" %in% colnames(dat)) {
    dat$MS1_bestfit <- NA  #create column as place holder if not exist
  }
  dat <- dat %>% dplyr::mutate(MS1_bestfit = ifelse(is.na(MS1_bestfit), MS1_Area, MS1_bestfit))

  #normalization using only HUMAN proteins
  datH <- dat %>%
    dplyr::filter(!grepl("YEAST", prot)) %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(med_H = median(.data[[quant_col]], na.rm = TRUE), .groups = "drop")
  
  dat <- dat %>%
    left_join(datH, by = "run_chan") %>%
    dplyr::mutate(med_norm = log2(.data[[quant_col]]) - log2(med_H)) %>%
    dplyr::group_by(seqcharge) %>%
    dplyr::mutate(seqcharge_mean = mean(.data[[quant_col]], na.rm = TRUE)) %>%
    dplyr::mutate(norm = med_norm - mean(med_norm, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(norm) & norm != 0 & is.finite(norm))
  
  #loess only needs to be run per time channel. Its an ionization/column specific effect.
  dat$run_timechannel <- paste0(dat$Run,"_",dat$time_channel)
  unique_run_chans <- unique(dat$run_timechannel)
  
  dat$loess_adjustment <- NA_real_
  dat$before_adj <- dat$norm
  
  for (chan in unique_run_chans) {
    
    dat_human <- dat %>% dplyr::filter(run_timechannel == chan & grepl("HUMAN", prot))
    
    if (nrow(dat_human) > 50) {  
      loess_model <- loess(norm ~ rt, data = dat_human, span = 0.1, na.action = na.exclude)
      
      dat[dat$run_timechannel == chan, "loess_adjustment"] <- predict(loess_model, 
                                                                      newdata = dat[dat$run_timechannel == chan, c("rt")])
    }
  }
  
  dat <- dat %>%
    dplyr::mutate(norm_original = norm)
  dat <- dat %>%
    dplyr::mutate(norm = norm - loess_adjustment)
  dat$quant_val <- (2^dat$norm)*dat$seqcharge_mean
  dat$quant_val_beforeAdj <- (2^dat$before_adj)*dat$seqcharge_mean
  
  return(dat)
}

quant_multi_single_emitter <- function(IO, uPAC, grouper="seqcharge", datatype = "precursor", min_num_reps=2, alphascaler=1){
  
  IO$implment <- "Multi"
  uPAC$implment <- "Single"
  
  normed_quant <- rbind.fill(IO, uPAC)
  
  normed_quant <- normed_quant %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tp_pd", #tp_pd
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))
  
  normed_quant$ProteinName <- normed_quant$prot
  normed_quant <- normed_quant[normed_quant$Repeatability==F,] 
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",prot),"Human",
                                                                    ifelse(grepl("YEAST",prot),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))
  
  tp <- normed_quant[normed_quant$Type=="tp",] %>% dplyr::filter(Rep == 2) #get the middle set of triplicates
  
  joined <- tp
  normed_quant_avg_ints <- joined %>% dplyr::group_by(Type, .data[[grouper]], Sample,implment) %>%
    dplyr::summarize("quant_avg" = mean(quant_val,na.rm=T)) %>% dplyr::rename("denom_sample"="Sample")
  
  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined <- joined %>%
    left_join(human_medians, by = "run_chan") %>%  
    dplyr::mutate(norm = norm-median_Human) %>%  
    dplyr::select(-median_Human) 
  
  
  joined_quant <- joined %>% dplyr::group_by(Type, Sample,.data[[grouper]],species,implment) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_num_reps-1)) %>% #require atleast 2/3 
    dplyr::summarize(med_quant = median(norm,na.rm=T)) #usually median
  
  dynamic_column <- grouper  
  
  formula_str <- paste(dynamic_column, "species", sep = "+") 
  joined_quant_d <- dcast(
    joined_quant, 
    formula = as.formula(paste(formula_str, "~ Type+Sample+implment")), 
    value.var = "med_quant"
  )
  
  
  joined_quant_d$tp_BC_single <- joined_quant_d$tp_B_Single-joined_quant_d$tp_C_Single
  joined_quant_d$tp_BC_multi <- joined_quant_d$tp_B_Multi-joined_quant_d$tp_C_Multi
  
  joined_quant_d$tp_AC_single <- joined_quant_d$tp_A_Single-joined_quant_d$tp_C_Single
  joined_quant_d$tp_AC_multi <- joined_quant_d$tp_A_Multi-joined_quant_d$tp_C_Multi
  
  joined_quant_d$tp_AB_single <- joined_quant_d$tp_A_Single-joined_quant_d$tp_B_Single
  joined_quant_d$tp_AB_multi <- joined_quant_d$tp_A_Multi-joined_quant_d$tp_B_Multi
  
  
  joined_quant_m <- joined_quant_d %>% dplyr::select(.data[[grouper]],species, 
                                                     tp_BC_single, tp_BC_multi,
                                                     tp_AC_single, tp_AC_multi,
                                                     tp_AB_single, tp_AB_multi
  ) %>% melt()
  
  
  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined_quant_m <- joined_quant_m %>%
    left_join(human_medians, by = "variable") %>%  
    dplyr::mutate(value = value-median_Human) %>%  
    dplyr::select(-median_Human) %>% na.omit()
  
  joined_quant_m <- joined_quant_m %>% 
    dplyr::mutate(denom_sample = ifelse(grepl("AB",variable),"B","C"))
  joined_quant_m$Type <- sub("_.*", "", joined_quant_m$variable)
  
  ##### AB
  joined_BC <- joined_quant_m[grepl("BC", joined_quant_m$variable),]
  
  int_BC <- joined_BC %>% na.omit() %>% count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  
  
  joined_BC_int <- joined_BC %>% dplyr::filter(.data[[grouper]]%in%int_BC[[grouper]]) %>% na.omit()
  joined_BC_int$variable <- factor(joined_BC_int$variable, 
                                   levels=c("tp_BC_multi", "tp_BC_single"
                                            
                                   ))
  
  
  
  y_d0d4 <- joined_BC_int$value %>% sort()
  joined_BC_int <- joined_BC_int %>% dplyr::add_count(species)
  joined_BC_int <- joined_BC_int %>% dplyr::mutate("Alpha" = ifelse(grepl("Human",species),0.05,0.17))
  
  
  
  medians <- joined_BC_int %>%
    dplyr::group_by(species, variable) %>%
    dplyr::summarize(median_value = median(value), .groups = "drop")
  
  ggplot(joined_BC_int, aes(fill=species, y=value)) +
    facet_grid(~variable) +
    geom_density(alpha=0.75) + 
    geom_hline(data=medians, aes(yintercept=median_value), linetype="solid", color="black") + 
    geom_hline(yintercept = 3, linetype="dashed", color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + 
    labs(subtitle=paste0("n = ",nrow(int_BC)," precursors"))+
    theme_bw() +
    scale_fill_manual(values=c("#68aabc","#e08a6c"))+
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) #+ ylim(y_d0d4_val)
  
  ggsave(paste0("TP_single_multi_emtter_BC",datatype,".pdf"),width=2.8,height=4)
  
  
  #####
  ##### AB
  joined_BC <- joined_quant_m[grepl("AB", joined_quant_m$variable),]
  
  int_BC <- joined_BC %>% na.omit() %>% count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  
  
  joined_BC_int <- joined_BC %>% dplyr::filter(.data[[grouper]]%in%int_BC[[grouper]]) %>% na.omit()
  joined_BC_int$variable <- factor(joined_BC_int$variable, 
                                   levels=c("tp_AB_multi", "tp_AB_single"
                                            
                                   ))
  
  
  
  y_d0d4 <- joined_BC_int$value %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.0025,.9975)))
  
  joined_BC_int <- joined_BC_int %>% dplyr::add_count(species)
  joined_BC_int <- joined_BC_int %>% dplyr::mutate("Alpha" = ifelse(grepl("Human",species),0.05,0.17))
  
  
  
  medians <- joined_BC_int %>%
    dplyr::group_by(species, variable) %>%
    dplyr::summarize(median_value = median(value), .groups = "drop")
  
  ggplot(joined_BC_int, aes(fill=species, y=value)) +
    facet_grid(~variable) +
    geom_density(alpha=0.75) + 
    geom_hline(data=medians, aes(yintercept=median_value), linetype="solid", color="black") + 
    geom_hline(yintercept = -2, linetype="dashed", color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + 
    labs(subtitle=paste0("n = ",nrow(int_BC)," precursors"))+
    theme_bw() +
    scale_fill_manual(values=c("#68aabc","#e08a6c"))+
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + ylim(y_d0d4_val)
  
  ggsave(paste0("TP_single_multi_emtter_AB",datatype,".pdf"),width=2.8,height=4)
  
  
  ##### AC
  joined_BC <- joined_quant_m[grepl("AC", joined_quant_m$variable),]
  
  int_BC <- joined_BC %>% na.omit() %>% count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  
  
  joined_BC_int <- joined_BC %>% dplyr::filter(.data[[grouper]]%in%int_BC[[grouper]]) %>% na.omit()
  joined_BC_int$variable <- factor(joined_BC_int$variable, 
                                   levels=c("tp_AC_multi", "tp_AC_single"
                                            
                                   ))
  
  
  
  y_d0d4 <- joined_BC_int$value %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.0025,.9975)))
  
  joined_BC_int <- joined_BC_int %>% dplyr::add_count(species)
  joined_BC_int <- joined_BC_int %>% dplyr::mutate("Alpha" = ifelse(grepl("Human",species),0.05,0.17))
  
  
  
  medians <- joined_BC_int %>%
    dplyr::group_by(species, variable) %>%
    dplyr::summarize(median_value = median(value), .groups = "drop")
  
  ggplot(joined_BC_int, aes(fill=species, y=value)) +
    facet_grid(~variable) +
    geom_density(alpha=0.75) + 
    geom_hline(data=medians, aes(yintercept=median_value), linetype="solid", color="black") + 
    geom_hline(yintercept = 1, linetype="dashed", color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + 
    labs(subtitle=paste0("n = ",nrow(int_BC)," precursors"))+
    theme_bw() +
    scale_fill_manual(values=c("#68aabc","#e08a6c"))+
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + ylim(y_d0d4_val)
  
  ggsave(paste0("TP_single_multi_emtter_AC",datatype,".pdf"),width=2.8,height=4)
  
}

single_run_quant <- function(normed_quant, use_best_channel = F,quant_col="coeff",alphascaler=1){
  
  normed_quant$ProteinName <- normed_quant$prot
  normed_quant <- normed_quant[normed_quant$Repeatability==F,] #remove Runs that are just for computing CVs
  normed_quant1 <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",prot),"Human",
                                                                     ifelse(grepl("YEAST",prot),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))
  normed_quant1_noplex <- normed_quant1[(normed_quant1$timePlex==F) & (normed_quant1$Label=="LF"),]
  normed_quant1_plex <- normed_quant1[(normed_quant1$timePlex==T)|(normed_quant1$Label!="LF"),]
  
  uni_runs <- unique(normed_quant1_plex$Run)  
  for(i in 1:length(uni_runs)){
    subDir <- uni_runs[i]
    dir.create(file.path(subDir)) #create folder
    run_of_interest <- normed_quant1[normed_quant1$Run==uni_runs[i],]
    
    uni_runchans <- unique(run_of_interest$run_chan)
    
    for(j in 1:length(unique(uni_runchans))){
      
      dat_singleSample <- run_of_interest[run_of_interest$run_chan==uni_runchans[j],]
      ggplot(dat_singleSample, aes(x=rt, y=norm)) + geom_pointdensity(alpha=0.1) + theme_classic() + 
        labs(x="RT, min", y= "Log2, relative precursor quant", title=paste0("After adjustment (",uni_runchans[j],")")) +
        ylim(-8,8) +scale_color_viridis(option="mako")+theme(legend.position = "none")
      ggsave(paste0(subDir,"/AFTER_RT_adj_",uni_runchans[j],".pdf"),width=5.5,height=4) 
      
      ggplot(dat_singleSample, aes(x=rt, y=before_adj)) + geom_pointdensity(alpha=0.1) + theme_classic() + 
        labs(x="RT, min", y= "Log2, relative precursor quant", title=paste0("Before adjustment (",uni_runchans[j],")")) +
        geom_point(aes(x=rt,y=loess_adjustment),color="red",size=0.75,alpha=1) + ylim(-8,8) +scale_color_viridis(option="mako")+
        theme(legend.position = "none")
      ggsave(paste0(subDir,"/BEFORE_RT_adj_",uni_runchans[j],".pdf"),width=5.5,height=4)
      
    }
    uni_current <- unique(run_of_interest[1,])
    if(uni_current$timePlex==T){
      
      if(!grepl("LF",uni_current$Label)){ #plexDIA + timePlex
        
        for(k in 0:(length(unique(run_of_interest$time_channel))-1)){
          
          within_pD <- run_of_interest[run_of_interest$time_channel==k,]
          within_pD <- within_pD %>% dplyr::distinct(seqcharge, run_chan,.keep_all=T)
          within_pD_d <- dcast(within_pD, seqcharge+ProteinName+species~Sample, value.var = "norm")
          within_pD_d$AB <- within_pD_d$A-within_pD_d$B
          within_pD_d$BC <- within_pD_d$B-within_pD_d$C
          within_pD_d$AC <- within_pD_d$A-within_pD_d$C
          
          AB <- quant(within_pD_d, numerator = "A", denominator = "B", exp_FC=-2,
                      title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA, timePlex-",k), alphaScaler=alphascaler)
          ggsave(plot=AB, paste0(subDir,"/TP_PD_Quant_",uni_runs[i],"_AB_",k,"_timeplex_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
          
          AC <- quant(within_pD_d, numerator = "A", denominator = "C", exp_FC=1,
                      title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA, timePlex-",k),alphaScaler=alphascaler)
          ggsave(plot=AC, paste0(subDir,"/TP_PD_Quant_",uni_runs[i],"_","AC_",k,"_timeplex_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
          
          BC <- quant(within_pD_d, numerator = "B", denominator = "C", exp_FC=3,
                      title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA, timePlex-",k),alphaScaler=alphascaler)
          ggsave(plot=BC, paste0(subDir,"/TP_PD_Quant_",uni_runs[i],"_","BC_",k,"_timeplex_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
          
        }
        
      } else if((uni_current$timePlex==F) & (!grepl("LF",uni_current$Label))){

        within_pD <- run_of_interest#[run_of_interest$time_channel==k,]
        within_pD <- within_pD %>% dplyr::distinct(seqcharge, run_chan,.keep_all=T)
        within_pD_d <- dcast(within_pD, seqcharge+ProteinName+species~Sample, value.var = "norm")
        within_pD_d$AB <- within_pD_d$A-within_pD_d$B
        within_pD_d$BC <- within_pD_d$B-within_pD_d$C
        within_pD_d$AC <- within_pD_d$A-within_pD_d$C
        
        AB <- quant(within_pD_d, numerator = "A", denominator = "B", exp_FC=-2,
                    title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA"),alphaScaler=alphascaler)
        ggsave(plot=AB, paste0(subDir,"/PD_Quant_",uni_runs[i],"_AB_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
        
        AC <- quant(within_pD_d, numerator = "A", denominator = "C", exp_FC=1,
                    title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA"),alphaScaler=alphascaler)
        ggsave(plot=AC, paste0(subDir,"/PD_Quant_",uni_runs[i],"_","AC_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
        
        BC <- quant(within_pD_d, numerator = "B", denominator = "C", exp_FC=3,
                    title_actual = paste0("3-plexDIA & 3-timePlex: 8 ng HY... within plexDIA"),alphaScaler=alphascaler)
        ggsave(plot=BC, paste0(subDir,"/PD_Quant_",uni_runs[i],"_","BC_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
        
        # }
        
      } else{
        
        
        within_run <- run_of_interest %>% dplyr::distinct(seqcharge, run_chan,.keep_all=T)
        within_run$log2quant_val <- log2(within_run$quant_val)
        within_run_d <- dcast(within_run, seqcharge+ProteinName+species~Sample, value.var = "log2quant_val")
        
        within_run_d$AB <- within_run_d$A-within_run_d$B
        within_run_d$BC <- within_run_d$B-within_run_d$C
        within_run_d$AC <- within_run_d$A-within_run_d$C
        
        AB <- quant(within_run_d, numerator = "A", denominator = "B", exp_FC=-2,
                    title_actual = paste0("3-timePlex: 8 ng HY... within run"),alphaScaler=alphascaler)
        ggsave(plot=AB, paste0(subDir,"/TP_Quant_",uni_runs[i],"_","AB_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
        
        
        AC <- quant(within_run_d, numerator = "A", denominator = "C", exp_FC=1,
                    title_actual = paste0("3-timePlex: 8 ng HY... within run"),alphaScaler=alphascaler)
        ggsave(plot=AC, paste0(subDir,"/TP_Quant_",uni_runs[i],"_","AC_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
        
        
        BC <- quant(within_run_d, numerator = "B", denominator = "C", exp_FC=3,
                    title_actual = paste0("3-timePlex: 8 ng HY... within run"),alphaScaler=alphascaler)
        ggsave(plot=BC, paste0(subDir,"/TP_Quant_",uni_runs[i],"_","BC_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=10,height=3.3)
        
        
        #### before RT adjustment:
        within_run <- run_of_interest %>% dplyr::distinct(seqcharge, run_chan,.keep_all=T)
        within_run$log2quant_val_beforeAdj <- log2(within_run$quant_val_beforeAdj)
        within_run_d <- dcast(within_run, seqcharge+ProteinName+species~Sample, value.var = "log2quant_val_beforeAdj")
        
        within_run_d$AB <- within_run_d$A-within_run_d$B
        within_run_d$BC <- within_run_d$B-within_run_d$C
        within_run_d$AC <- within_run_d$A-within_run_d$C
        
        AB <- quant(within_run_d, numerator = "A", denominator = "B", exp_FC=-2,
                    title_actual = paste0("3-timePlex: 8 ng HY... within run"),alphaScaler=alphascaler)
        ggsave(plot=AB, paste0(subDir,"/TP_Quant_",uni_runs[i],"_","AB_bestchannel_",use_best_channel,"_beforeAdj",quant_col,".pdf"),width=10,height=3.3)
        
        
        AC <- quant(within_run_d, numerator = "A", denominator = "C", exp_FC=1,
                    title_actual = paste0("3-timePlex: 8 ng HY... within run"),alphaScaler=alphascaler)
        ggsave(plot=AC, paste0(subDir,"/TP_Quant_",uni_runs[i],"_","AC_bestchannel_",use_best_channel,"_beforeAdj",quant_col,".pdf"),width=10,height=3.3)
        
        
        BC <- quant(within_run_d, numerator = "B", denominator = "C", exp_FC=3,
                    title_actual = paste0("3-timePlex: 8 ng HY... within run"),alphaScaler=alphascaler)
        ggsave(plot=BC, paste0(subDir,"/TP_Quant_",uni_runs[i],"_","BC_bestchannel_",use_best_channel,"_beforeAdj",quant_col,".pdf"),width=10,height=3.3)
        
        
      }
      
    }}
  
  
  
  for(p in 1:(length(unique(normed_quant1_noplex$Rep)))){
    subDir <- paste0("NoPlex_LF_Rep_",p)
    dir.create(file.path(subDir)) #create folder

    reps_noplex <- normed_quant1_noplex %>% dplyr::distinct(Run,.keep_all=T) %>%
      dplyr::group_by(Rep,Sample) %>%
      dplyr::mutate(Rep = n()) %>%
      dplyr::ungroup()
    within_Rep <- normed_quant1_noplex[normed_quant1_noplex$Rep==p,]
    if(length(unique(within_Rep$Sample))==3){
      
      within_Rep <- dcast(within_Rep, seqcharge+ProteinName+species~Sample, value.var = "norm")
      within_Rep$AB <- within_Rep$A-within_Rep$B
      within_Rep$BC <- within_Rep$B-within_Rep$C
      within_Rep$AC <- within_Rep$A-within_Rep$C
      
      AB <- quant(within_Rep, numerator = "A", denominator = "B", exp_FC=-2,
                  title_actual = paste0("LF NoPlex: 8 ng HY... Rep ", p),alphaScaler=alphascaler)
      ggsave(plot=AB, paste0(subDir,"/LF_NoPlex_Quant_Rep",uni_runs[p],"_","AB_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=9,height=3.3)
      
      
      AC <- quant(within_Rep, numerator = "A", denominator = "C", exp_FC=1,
                  title_actual = paste0("LF NoPlex: 8 ng HY... Rep ", p),alphaScaler=alphascaler)
      ggsave(plot=AC, paste0(subDir,"/LF_NoPlex_Quant_Rep",uni_runs[p],"_","AC_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=9,height=3.3)
      
      
      BC <- quant(within_Rep, numerator = "B", denominator = "C", exp_FC=3,
                  title_actual = paste0("LF NoPlex: 8 ng HY... Rep ",p),alphaScaler=alphascaler)
      ggsave(plot=BC, paste0(subDir,"/LF_NoPlex_Quant_Rep",uni_runs[p],"_","BC_bestchannel_",use_best_channel,"_",quant_col,".pdf"),width=9,height=3.3)
      
    }
    
  }
}

quant_all_comparisons_v2 <- function(normed_quant, grouper="seqcharge", datatype = "precursor",alphascaler=1,min_num_reps=2){
  
  normed_quant <- normed_quant %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tp_pd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))
  
  normed_quant$ProteinName <- normed_quant$prot
  normed_quant <- normed_quant[normed_quant$Repeatability==F,] #remove Runs that are just for computing CVs
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",prot),"Human",
                                                                    ifelse(grepl("YEAST",prot),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))
  
  np <- normed_quant[normed_quant$Type=="np",]
  
  pD <- normed_quant[normed_quant$Type=="pd",] %>% dplyr::filter(grepl("JD0425|JD0426|JD0427",run_chan))
  tp_pd <- normed_quant[normed_quant$Type=="tp_pd",] %>% dplyr::filter(Rep == 2)
  tp <- normed_quant[normed_quant$Type=="tp",] %>% dplyr::filter(Rep == 2)
  
  joined <- rbind(np, pD, tp_pd, tp)
  normed_quant_avg_ints <- joined %>% dplyr::group_by(Type, .data[[grouper]], Sample) %>%
    dplyr::summarize("quant_avg" = mean(quant_val,na.rm=T)) %>% dplyr::rename("denom_sample"="Sample")
  
  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined <- joined %>%
    left_join(human_medians, by = "run_chan") %>%  
    dplyr::mutate(norm = norm-median_Human) %>%  
    dplyr::select(-median_Human) 
  
  
  joined_quant <- joined %>% dplyr::group_by(Type, Sample,.data[[grouper]],species) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_num_reps-1)) %>% #require atleast 1/3 
    dplyr::summarize(med_quant = median(norm,na.rm=T)) #usually median
  dynamic_column <- grouper  # This can change dynamically
  
  formula_str <- paste(dynamic_column, "species", sep = "+") # Keeps these in rows
  joined_quant_d <- dcast(
    joined_quant, 
    formula = as.formula(paste(formula_str, "~ Type+Sample")), 
    value.var = "med_quant"
  )
  joined_quant_d$np_AB <- joined_quant_d$np_A-joined_quant_d$np_B
  joined_quant_d$np_AC <- joined_quant_d$np_A-joined_quant_d$np_C
  joined_quant_d$np_BC <- joined_quant_d$np_B-joined_quant_d$np_C
  
  joined_quant_d$tp_AB <- joined_quant_d$tp_A-joined_quant_d$tp_B
  joined_quant_d$tp_AC <- joined_quant_d$tp_A-joined_quant_d$tp_C
  joined_quant_d$tp_BC <- joined_quant_d$tp_B-joined_quant_d$tp_C
  
  joined_quant_d$pd_AB <- joined_quant_d$pd_A-joined_quant_d$pd_B
  joined_quant_d$pd_AC <- joined_quant_d$pd_A-joined_quant_d$pd_C
  joined_quant_d$pd_BC <- joined_quant_d$pd_B-joined_quant_d$pd_C
  
  joined_quant_d$tp_pd_AB <- joined_quant_d$tp_pd_A-joined_quant_d$tp_pd_B
  joined_quant_d$tp_pd_AC <- joined_quant_d$tp_pd_A-joined_quant_d$tp_pd_C
  joined_quant_d$tp_pd_BC <- joined_quant_d$tp_pd_B-joined_quant_d$tp_pd_C
  
  joined_quant_m <- joined_quant_d %>% dplyr::select(.data[[grouper]],species, 
                                                     np_AB, np_AC, np_BC,
                                                     tp_AB, tp_AC, tp_BC,
                                                     pd_AB, pd_AC, pd_BC,
                                                     tp_pd_AB, tp_pd_AC, tp_pd_BC) %>% melt()
  
  
  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined_quant_m <- joined_quant_m %>%
    left_join(human_medians, by = "variable") %>%  
    dplyr::mutate(value = value-median_Human) %>%  
    dplyr::select(-median_Human) %>% na.omit()
  
  joined_quant_m <- joined_quant_m %>% 
    dplyr::mutate(denom_sample = ifelse(grepl("AB",variable),"B","C"))
  joined_quant_m$Type <- sub("_.*", "", joined_quant_m$variable)
  joined_quant_m <- joined_quant_m %>% left_join(normed_quant_avg_ints, by =c("Type","denom_sample",paste0(grouper)))
  
  allcounts  <- joined_quant_m %>% na.omit() %>% count(variable)
  
  
  ##### AB
  joined_AB <- joined_quant_m[grepl("AB", joined_quant_m$variable),]
  
  
  counts_AB <- allcounts[grepl("AB",allcounts$variable),]
  
  int_AB <- joined_AB %>% na.omit() %>% count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  
  df_int <- data.frame(variable = "Intersected", n=nrow(int_AB))
  counts_AB <- rbind(counts_AB, df_int)
  
  
  counts_AB$variable <- factor(counts_AB$variable, 
                               levels=c("np_AB", "pd_AB",
                                        "tp_AB", "tp_pd_AB",
                                        "Intersected"))
  
  counts_AB_plot <- ggplot(counts_AB, aes(x=as.factor(variable), y=n, fill=variable)) + 
    geom_bar(stat="identity", colour="black", alpha =0.7) +
    geom_text(data=counts_AB, aes(x=as.factor(variable), y=n+300, label=n), col='black', size=4.5) +
    theme_classic() + 
    theme(#axis.text.x = element_blank(),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14),
      axis.ticks.x = element_blank(),
      # plot.title = element_text(size=14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.title = element_text(size = 13.5, hjust = 0.5, face="bold"),
      legend.position = "none",
      strip.text.x = element_text(size=12, face = "italic"),
      legend.text = element_text(size=12))+
    labs(x = "", y = "Protein Ratios (A/B)", fill = "") +
    scale_fill_manual(values = c("#413c58", "#b6db8a", "#a3c1c1", "#d87c7c","grey40"))
  
  
  
  joined_AB_int <- joined_AB %>% dplyr::filter(.data[[grouper]]%in%int_AB[[grouper]]) %>% na.omit()
  joined_AB_int$variable <- factor(joined_AB_int$variable, 
                                   levels=c("np_AB", "pd_AB",
                                            "tp_AB", "tp_pd_AB",
                                            "Intersected"))
  
  
  x_d0d4 <- log2(joined_AB_int$quant_avg) %>% sort()
  x_d0d4_val <- as.numeric(quantile(x_d0d4,probs=c(.005,.995)))
  y_d0d4 <- joined_AB_int$value %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.005,.995)))
  
  joined_AB_int <- joined_AB_int %>% add_count(species)
  joined_AB_int <- joined_AB_int %>% dplyr::mutate("Alpha" = ifelse(grepl("Human",species),0.05,0.15))
  
  scatter <- ggplot(joined_AB_int, aes(x=log2(quant_avg), y=value, color=species, alpha=Alpha*alphascaler)) + geom_point(size=1) + 
    facet_grid(~variable)+
    geom_hline(yintercept = -2,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_bw()+
    scale_color_manual(values=c("#68aabc","#e08a6c")) +
    scale_alpha_identity() +
    labs(#title = title_actual, 
      y=paste0("Log2, A/B"),
      x=paste0("Log2, B"))+
    theme(legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(y_d0d4_val), xlim = c(x_d0d4_val))
  
  joined_AB_int$placeholder <- "placeholder"
  box <- ggplot(joined_AB_int, 
                aes(fill=species,y=value,x=variable)) +geom_boxplot(alpha=0.85,outlier.alpha=alphascaler*0.01,outlier.size = 1)+
    facet_grid(~placeholder) +
    geom_hline(yintercept = -2,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_bw()+
    scale_fill_manual(values=c("#68aabc","#e08a6c")) +
    labs(#title = title_actual, 
      #subtitle = ""
    ) +
    theme(axis.title.y = element_blank(),
          legend.position="none",
          axis.text.y = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(y_d0d4_val))
  
  quantplot_AB <- ggarrange(counts_AB_plot,scatter, box, #dense,
                            ncol = 3, nrow = 1,
                            widths = c(0.17,0.8, 0.17))
  ggsave(paste0("quantplot_AB_",datatype,".pdf"),plot=quantplot_AB,width=13,height=2.3)
  
  
  
  ########## BC
  joined_AB <- joined_quant_m[grepl("BC", joined_quant_m$variable),]
  
  
  counts_AB <- allcounts[grepl("BC",allcounts$variable),]
  
  int_AB <- joined_AB %>% na.omit() %>% count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  df_int <- data.frame(variable = "Intersected", n=nrow(int_AB))
  counts_AB <- rbind(counts_AB, df_int)
  
  
  counts_AB$variable <- factor(counts_AB$variable, 
                               levels=c("np_BC", "pd_BC",
                                        "tp_BC", "tp_pd_BC",
                                        "Intersected"))
  
  counts_AB_plot <- ggplot(counts_AB, aes(x=as.factor(variable), y=n, fill=variable)) + 
    geom_bar(stat="identity", colour="black", alpha =0.7) +
    geom_text(data=counts_AB, aes(x=as.factor(variable), y=n+300, label=n), col='black', size=4.5) +
    theme_classic() + 
    theme(#axis.text.x = element_blank(),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14),
      axis.ticks.x = element_blank(),
      # plot.title = element_text(size=14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.title = element_text(size = 13.5, hjust = 0.5, face="bold"),
      legend.position = "none",
      strip.text.x = element_text(size=12, face = "italic"),
      legend.text = element_text(size=12))+
    labs(x = "", y = "Protein Ratios (B/C)", fill = "")+
    scale_fill_manual(values = c("#413c58", "#b6db8a", "#a3c1c1", "#d87c7c","grey40"))
  
  
  joined_AB_int <- joined_AB %>% dplyr::filter(.data[[grouper]]%in%int_AB[[grouper]]) %>% na.omit()
  joined_AB_int$variable <- factor(joined_AB_int$variable, 
                                   levels=c("np_BC", "pd_BC",
                                            "tp_BC", "tp_pd_BC",
                                            "Intersected"))
  
  
  x_d0d4 <- log2(joined_AB_int$quant_avg) %>% sort()
  x_d0d4_val <- as.numeric(quantile(x_d0d4,probs=c(.005,.995)))
  y_d0d4 <- joined_AB_int$value %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.005,.995)))
  
  joined_AB_int <- joined_AB_int %>% add_count(species)
  joined_AB_int <- joined_AB_int %>% dplyr::mutate("Alpha"=0.5*(1-n/nrow(.)))
  joined_AB_int <- joined_AB_int %>% dplyr::mutate("Alpha" = ifelse(grepl("Human",species),0.05,0.15))
  
  scatter <- ggplot(joined_AB_int, aes(x=log2(quant_avg), y=value, color=species, alpha=Alpha*alphascaler)) + geom_point(size=1) + 
    facet_grid(~variable)+
    geom_hline(yintercept = 3,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_bw()+
    scale_color_manual(values=c("#68aabc","#e08a6c")) +
    scale_alpha_identity() +
    labs(#title = title_actual, 
      y=paste0("Log2, B/C"),
      x=paste0("Log2, C"))+
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(y_d0d4_val), xlim = c(x_d0d4_val))
  
  joined_AB_int$placeholder <- "placeholder"
  box <- ggplot(joined_AB_int, 
                aes(fill=species,y=value,x=variable)) +geom_boxplot(alpha=0.85,outlier.alpha=alphascaler*0.01,outlier.size = 1)+
    facet_grid(~placeholder) +
    geom_hline(yintercept = 3,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_bw()+
    scale_fill_manual(values=c("#68aabc","#e08a6c")) +
    labs(#title = title_actual, 
      
      #subtitle = ""
    ) +
    theme(axis.title.y = element_blank(),
          legend.position="none",
          axis.text.y = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(y_d0d4_val))
  
  quantplot_BC <- ggarrange(counts_AB_plot,scatter, box, #dense,
                            ncol = 3, nrow = 1,
                            widths = c(0.17,0.8, 0.17))
  ggsave(paste0("quantplot_BC_",datatype,".pdf"),plot=quantplot_BC,width=13,height=2.3)
  
  
  ########## AC
  joined_AB <- joined_quant_m[grepl("AC", joined_quant_m$variable),]
  
  
  counts_AB <- allcounts[grepl("AC",allcounts$variable),]
  
  int_AB <- joined_AB %>% na.omit() %>% count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  
  df_int <- data.frame(variable = "Intersected", n=nrow(int_AB))
  counts_AB <- rbind(counts_AB, df_int)
  
  
  counts_AB$variable <- factor(counts_AB$variable, 
                               levels=c("np_AC", "pd_AC",
                                        "tp_AC", "tp_pd_AC",
                                        "Intersected"))
  
  counts_AB_plot <- ggplot(counts_AB, aes(x=as.factor(variable), y=n, fill=variable)) + 
    geom_bar(stat="identity", colour="black", alpha =0.7) +
    geom_text(data=counts_AB, aes(x=as.factor(variable), y=n+300, label=n), col='black', size=4.5) +
    #scale_fill_manual(values = c("grey15", "grey35","grey50","grey65","grey80"))+
    theme_classic() + 
    theme(#axis.text.x = element_blank(),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14),
      axis.ticks.x = element_blank(),
      # plot.title = element_text(size=14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.title = element_text(size = 13.5, hjust = 0.5, face="bold"),
      legend.position = "none",
      strip.text.x = element_text(size=12, face = "italic"),
      legend.text = element_text(size=12))+
    labs(x = "", y = "Protein Ratios (B/C)", fill = "")+
    scale_fill_manual(values = c("#413c58", "#b6db8a", "#a3c1c1", "#d87c7c","grey40"))
  
  
  joined_AB_int <- joined_AB %>% dplyr::filter(.data[[grouper]]%in%int_AB[[grouper]]) %>% na.omit()
  joined_AB_int$variable <- factor(joined_AB_int$variable, 
                                   levels=c("np_AC", "pd_AC",
                                            "tp_AC", "tp_pd_AC",
                                            "Intersected"))
  
  
  x_d0d4 <- log2(joined_AB_int$quant_avg) %>% sort()
  x_d0d4_val <- as.numeric(quantile(x_d0d4,probs=c(.005,.995)))
  y_d0d4 <- joined_AB_int$value %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(.005,.995)))
  
  joined_AB_int <- joined_AB_int %>% add_count(species)
  joined_AB_int <- joined_AB_int %>% dplyr::mutate("Alpha"=0.5*(1-n/nrow(.)))
  joined_AB_int <- joined_AB_int %>% dplyr::mutate("Alpha" = ifelse(grepl("Human",species),0.05,0.15))
  
  scatter <- ggplot(joined_AB_int, aes(x=log2(quant_avg), y=value, color=species, alpha=Alpha*alphascaler)) + geom_point(size=1) + 
    facet_grid(~variable)+
    geom_hline(yintercept = 1,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_bw()+
    scale_color_manual(values=c("#68aabc","#e08a6c")) +
    scale_alpha_identity() +
    labs(
      y=paste0("Log2, B/C"),
      x=paste0("Log2, C"))+
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(y_d0d4_val), xlim = c(x_d0d4_val))
  
  joined_AB_int$placeholder <- "placeholder"
  box <- ggplot(joined_AB_int, 
                aes(fill=species,y=value,x=variable)) +geom_boxplot(alpha=0.85,outlier.alpha=alphascaler*0.01,outlier.size = 1)+
    facet_grid(~placeholder) +
    geom_hline(yintercept = 1,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_bw()+
    scale_fill_manual(values=c("#68aabc","#e08a6c")) +
    labs(#title = title_actual, 
      
      #subtitle = ""
    ) +
    theme(axis.title.y = element_blank(),
          legend.position="none",
          axis.text.y = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(y_d0d4_val))
  
  quantplot_AC <- ggarrange(counts_AB_plot,scatter, box, #dense,
                            ncol = 3, nrow = 1,
                            widths = c(0.17,0.8, 0.17))
  ggsave(paste0("quantplot_AC_",datatype,".pdf"),plot=quantplot_AC,width=13,height=2.3)
  
  
  
  
}


compare_quants <- function(dat1=normed_quant_IO_MS1_bestfit_relaxed_SR1, 
                           dat2=normed_quant_IO_relaxed_coeff,
                           quant1="MS1_bestfit", 
                           quant2="coeff", min_num_reps=2, grouper="seqcharge"){
  
  dat1$quant <- "quant1"
  dat2$quant <- "quant2"
  normed_quant <- rbind.fill(dat1,dat2)
  
  normed_quant <- normed_quant %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tppd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))
  
  normed_quant$ProteinName <- normed_quant$prot
  normed_quant <- normed_quant[normed_quant$Repeatability==F,] #remove Runs that are just for computing CVs
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",prot),"Human",
                                                                    ifelse(grepl("YEAST",prot),"Yeast","remove"))) %>% dplyr::filter(!grepl("remove",species))
  
  
  np <- normed_quant[normed_quant$Type=="np",]
  pD <- normed_quant[normed_quant$Type=="pd",] %>% dplyr::filter(grepl("JD0425|JD0426|JD0427",run_chan))
  tp_pd <- normed_quant[normed_quant$Type=="tppd",] %>% dplyr::filter(Rep == 2)
  tp <- normed_quant[normed_quant$Type=="tp",] %>% dplyr::filter(Rep == 2)
  
  joined <- rbind(np, pD, tp_pd, tp)
  normed_quant_avg_ints <- joined %>% dplyr::group_by(Type, .data[[grouper]], Sample, quant) %>%
    dplyr::summarize("quant_avg" = mean(quant_val,na.rm=T)) %>% dplyr::rename("denom_sample"="Sample")
  
  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan,quant) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined <- joined %>%
    left_join(human_medians, by = c("run_chan", "quant")) %>%  
    dplyr::mutate(norm = norm-median_Human) %>%  
    dplyr::select(-median_Human) 
  
  joined_quant <- joined %>% dplyr::group_by(Type, Sample,.data[[grouper]],species,quant) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_num_reps-1)) %>% #require atleast 1/3 
    dplyr::summarize(med_quant = median(norm,na.rm=T)) #usually median
  dynamic_column <- grouper  
  
  formula_str <- paste(dynamic_column, "species", sep = "+") 
  joined_quant_d <- dcast(
    joined_quant, 
    formula = as.formula(paste(formula_str, "~ Type+Sample+quant")), 
    value.var = "med_quant"
  )
  
  
  prefixes <- c("np", "tp", "pd", "tppd")
  comparisons <- list(c("A", "B"), c("A", "C"), c("B", "C"))
  suffixes <- c("quant1", "quant2")
  
  for (prefix in prefixes) {
    for (suffix in suffixes) {
      for (pair in comparisons) {
        colname <- paste0(prefix, "_", paste0(pair, collapse=""), "_", suffix)
        col1 <- paste0(prefix, "_", pair[1], "_", suffix)
        col2 <- paste0(prefix, "_", pair[2], "_", suffix)
        joined_quant_d[[colname]] <- joined_quant_d[[col1]] - joined_quant_d[[col2]]
      }
    }
  }
  
  comparison_cols <- unlist(lapply(prefixes, function(prefix) {
    unlist(lapply(suffixes, function(suffix) {
      sapply(comparisons, function(pair) {
        paste0(prefix, "_", paste0(pair, collapse=""), "_", suffix)
      })
    }))
  }))
  
  joined_quant_m <- joined_quant_d %>%
    select(all_of(c(grouper, "species", comparison_cols))) %>%
    melt()
  
  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined_quant_m <- joined_quant_m %>%
    left_join(human_medians, by = "variable") %>%  
    dplyr::mutate(value = value-median_Human) %>%  
    dplyr::select(-median_Human) %>% na.omit()
  
  joined_quant_m <- joined_quant_m %>% 
    dplyr::mutate(denom_sample = ifelse(grepl("AB",variable),"B","C"))
  joined_quant_m$Type <- sub("_.*", "", joined_quant_m$variable)
  joined_quant_m <- joined_quant_m %>% dplyr::mutate("theoretical" = 
                                                       ifelse((grepl("AB",variable) & species=="Yeast"),-2,
                                                              ifelse((grepl("BC",variable) & species=="Yeast"),3,
                                                                     ifelse((grepl("AC",variable) & species=="Yeast"),1, 0))))
  
  int_AB <- joined_quant_m %>% na.omit() %>% count(.data[[grouper]],as.factor(theoretical)) %>% dplyr::filter(n==max(n))
  
  joined_quant_m_int <- joined_quant_m[joined_quant_m$seqcharge%in%int_AB$seqcharge,] %>% na.omit()
  joined_quant_m_int$residual <- abs(joined_quant_m_int$value-joined_quant_m_int$theoretical)
  
  joined_quant_m_int$data_type <- sub("_.*", "", joined_quant_m_int$variable)
  
  joined_quant_m_int$quant_type <- sub(".*_", "", joined_quant_m_int$variable)
  
  y_d0d4 <- joined_quant_m_int$residual %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(0,.95)))
  
  extract_chars_before_underscore <- function(text) {
    matches <- str_extract_all(text, ".(?=_)")[[1]]
    return(matches)
  }
  
  joined_quant_m_int <- joined_quant_m_int %>%
    dplyr::mutate(end_aa = lapply(seqcharge, extract_chars_before_underscore))
  
  
  joined_quant_m_int <- joined_quant_m_int %>% dplyr::mutate("ends_K" = ifelse(grepl("K",end_aa),T,F))
  
  joined_quant_m_int <- joined_quant_m_int %>% 
    dplyr::mutate("quant_type_alt" = ifelse(quant_type=="quant2",paste0(quant_type,"_",ends_K),quant_type))
  

  ggplot(joined_quant_m_int, aes(x = residual, y = fct_rev(data_type))) +
    facet_grid(rows = vars(quant_type)) +
    geom_density_ridges2(
      aes(fill = data_type), 
      alpha = 0.75, 
      scale = 1.0, 
      color = "black"
    ) + labs(subtitle=paste0(nrow(joined_quant_m_int)/8," prec ratios"))+
    stat_summary(
      aes(y = fct_rev(data_type)), 
      fun = median, 
      geom = "crossbar", 
      width = 0.33, 
      color = "black", 
      fatten = 1.5
    ) +
    coord_cartesian(xlim = y_d0d4_val) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    scale_fill_manual(values = c("#413c58", "#b6db8a", "#a3c1c1", "#d87c7c"))
  
  ggsave("MS1_MS2_quant_K_R.pdf",width=2.5,height=6)
  
}

join_implementations <- function(dat1=benchmarking_dat_IO$Amount_8, dat2=benchmarking_dat_IO_singlecolumn$Amount_8){
  dat1$Implementation <- "3cols"
  dat2$Implementation <- "1col"
  
  joined <- rbind.fill(dat1,dat2)
  return(joined)
}

tp_coverage_compare_implementation <- function(dat,extra="", use_best_channel = F, filter_level=0.01){
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(BestChannel_Qvalue<filter_level) %>% filter(Qvalue < 0.2)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<filter_level)
  }
  
  dat <- dat %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))
  
  if ("Implementation" %in% colnames(dat)) {
    dat <- dat %>% dplyr::filter(Type == "3-plexDIA" | Type=="LF\nNo-plex")
    dat$Type <- paste0(dat$Type,"_",dat$Implementation)
    
  }
  
  dat <- dat[dat$Repeatability ==F,]
  dat_keep <- dat
  dat <- dat %>% dplyr::mutate(Label = ifelse(grepl("3-plexDIA",Type),"mTRAQ",Label))
  
  quant <- dat %>% dplyr::filter(grepl("LF",Type)) %>% 
    group_by(Implementation,seqcharge,Sample) %>%
    dplyr::summarize("quant" = mean(MS1_Area)) %>% ungroup()
  
  quant_d <- dcast(quant, seqcharge+Sample~Implementation, value.var = "quant")
  quant_d$ratio <- quant_d$`3cols`/quant_d$`1col`
  
  y_d0d4 <- quant_d$ratio %>% sort()
  y_d0d4_val <- as.numeric(quantile(log2(y_d0d4),probs=c(0.005,.995)))
  ggplot(quant_d, aes(y = log2(ratio))) + 
    geom_histogram(fill="grey65",color="black",size=0.25) + 
    ylim(y_d0d4_val) + 
    geom_hline(aes(yintercept = median(log2(ratio), na.rm = TRUE)), 
               color = "red", linetype = "dashed", size = 1) + 
    theme_classic()
  ggsave(paste0(extra,"_","Change_in_Intensity","BestChannel_",use_best_channel,".pdf"),width=2,height=4)
  
  
  
  dat_counts <- dat %>%
    group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    ungroup()
  
  dat_counts$Sample <- factor(dat_counts$Sample, levels = c("A", "B", "C"))
  
  
  # Compute mean and standard error
  dat_summary <- dat_counts %>%
    group_by(Sample, Label, Amount, Type) %>%
    summarise(
      mean_count = mean(n, na.rm = TRUE),
      se = sd(n, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% ungroup()
  
  dat_summary <- dat_summary %>%
    dplyr::arrange(desc(Sample), Type) %>%
    
    dplyr::group_by(Type) %>%
    dplyr::mutate(
      stacked_mean = cumsum(mean_count)) %>% dplyr::ungroup() %>%
    dplyr::group_by(Label, Sample, Type) %>%
    dplyr::mutate(stacked_mean = ifelse(grepl("LF\nNo-plex",Type), mean_count,stacked_mean))%>%
    dplyr::mutate(
      er_max = stacked_mean + se,
      er_min = stacked_mean - se
    ) %>%
    dplyr::ungroup()
  
  df_text <- dat_summary %>%
    dplyr::group_by(Type, Sample,Label) %>%
    dplyr::summarise(
      text_value = round(max(er_max) - max(se)),
      text_pos = max(er_max) + 2000
    )

  my_colors <- c("#7ca982","#4e598c","#fed766")
  
  ggplot(dat_summary, aes(x = Type, y = mean_count, fill = Sample)) +
    scale_x_discrete(limits = c("LF\nNo-plex_1col", "LF\nNo-plex_3cols", "3-plexDIA_1col","3-plexDIA_3cols")) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex_1col"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 0.66),
             color = "black", width = 0.66, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex_3cols"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 0.66),
             color = "black", width = 0.66, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA_1col"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.22, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA_3cols"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.22, alpha = 1) +
    scale_y_continuous(labels = comma) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("3-plexDIA_1col","3-plexDIA_3cols")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.18, size = 0.4) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("LF\nNo-plex_1col","LF\nNo-plex_3cols")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.5, size = 0.4,
                  position = position_dodge(0.66)) +
    
    scale_fill_manual(values = my_colors) +
    theme_classic() +
    
    labs(y = "Precursors data points / Run", x = "", fill = "Samples") +
    
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = c(0.16, 0.9),
          legend.direction = "horizontal",
          legend.spacing.x = unit(1.0, 'line')) +
    
    guides(fill = guide_legend(override.aes = list(alpha = 0.7),
                               title = "Samples",
                               label.position = "bottom",
                               title.position = "top", title.vjust = 1, title.hjust = 0.5)) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("3-plexDIA_1col","3-plexDIA_3cols")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("LF\nNo-plex_1col","LF\nNo-plex_3cols")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3,
              width = 0.5, size = 0.4,
              position = position_dodge(width=0.66))
  
  ggsave(paste0(extra,"_","Prec_coverage_","BestChannel_",use_best_channel,".pdf"),width=6,height=5)
  
  
  ##################
  dat <- dat_keep %>% dplyr::mutate(Label = ifelse(grepl("3-plexDIA",Type),"mTRAQ",Label))
  
  dat_counts <- dat %>% dplyr::filter(Protein_Qvalue<filter_level) %>% dplyr::distinct(Sample, Run, Column, Amount, Label, Type, prot,.keep_all=T) %>% 
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>% dplyr::ungroup()
  
  dat_summary <- dat_counts %>%
    dplyr::group_by(Sample, Label, Amount, Type) %>%
    dplyr::summarise(
      mean_count = mean(n, na.rm = TRUE),
      se = sd(n, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% dplyr::ungroup()
  
  dat_summary <- dat_summary %>%
    dplyr::arrange(desc(Sample), Type) %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(
      stacked_mean = cumsum(mean_count)) %>% dplyr::ungroup() %>%
    dplyr::group_by(Label, Sample, Type) %>%
    dplyr::mutate(stacked_mean = ifelse(grepl("LF\nNo-plex",Type), mean_count,stacked_mean))%>%
    dplyr::mutate(
      er_max = stacked_mean + se,
      er_min = stacked_mean - se
    ) %>%
    dplyr::ungroup()
  
  df_text <- dat_summary %>%
    dplyr::group_by(Type, Sample,Label) %>%
    dplyr::summarise(
      text_value = round(max(er_max) - max(se)),
      text_pos = max(er_max) + 310
    )
  
  
  ggplot(dat_summary, aes(x = Type, y = mean_count, fill = Sample)) +
    scale_x_discrete(limits = c("LF\nNo-plex_1col", "LF\nNo-plex_3cols", "3-plexDIA_1col","3-plexDIA_3cols")) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex_1col"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 0.66),
             color = "black", width = 0.66, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex_3cols"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 0.66),
             color = "black", width = 0.66, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA_1col"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.22, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "3-plexDIA_3cols"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.22, alpha = 1) +
    scale_y_continuous(labels = comma) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("3-plexDIA_1col","3-plexDIA_3cols")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.18, size = 0.4) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("LF\nNo-plex_1col","LF\nNo-plex_3cols")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.5, size = 0.4,
                  position = position_dodge(0.66)) +
    
    scale_fill_manual(values = my_colors) +
    theme_classic() +
    
    labs(y = "Precursors data points / Run", x = "", fill = "Samples") +
    
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = c(0.16, 0.9),
          legend.direction = "horizontal",
          legend.spacing.x = unit(1.0, 'line')) +
    
    guides(fill = guide_legend(override.aes = list(alpha = 0.7),
                               title = "Samples",
                               label.position = "bottom",
                               title.position = "top", title.vjust = 1, title.hjust = 0.5)) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("3-plexDIA_1col","3-plexDIA_3cols")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("LF\nNo-plex_1col","LF\nNo-plex_3cols")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3,
              width = 0.5, size = 0.4,
              position = position_dodge(width=0.66))
  
  ggsave(paste0(extra,"_","Prot_coverage_","BestChannel_",use_best_channel,".pdf"),width=6,height=5)
  
  return(dat_counts)
  
  
}

Jmod_FWHMs <- function(uPAC=benchmarking_dat_uPAC$Amount_8, IO=benchmarking_dat_IO$Amount_8){
  
  uPAC$system <- "uPAC"
  IO$system <- "IO"
  joined <- rbind.fill(uPAC, IO)
  joined <- joined %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    )) %>% dplyr::filter(Type =="3-plexDIA\n3-timePlex")
  joined <- joined %>% dplyr::filter(Qvalue<0.01 & decoy==F)
  joined <- joined %>%
    rowwise() %>% 
    dplyr::mutate(FWHM_scans = compute_fwhm(plexfittrace_all))
  
  joined$FWHM <- joined$FWHM_scans*1.3 #seconds
  joined$Column <- joined$Column+1
  
  med_joined <- joined %>% dplyr::group_by(seqcharge,system,Column) %>%
    dplyr::summarize(med_FWHM = median(FWHM,na.rm=T),
                     med_quant = median(MS1_Area,na.rm=T)) %>% dplyr::ungroup() %>% na.omit() %>%
    dplyr::add_count(seqcharge) %>% dplyr::filter(n==6)
  
  y_d0d4 <- med_joined$med_FWHM %>% sort()
  y_d0d4_val <- as.numeric(quantile(y_d0d4,probs=c(0.0,.995)))
  
  ggplot(med_joined, aes(x = med_FWHM, y = as.factor(Column))) +
    facet_grid(rows = vars(system)) +
    geom_density_ridges2(
      aes(), 
      alpha = 0.75, 
      scale = 1.0, 
      color = "black",
      fill="grey70"
    ) +
    stat_summary(
      aes(y = as.factor(Column)), 
      fun = median, 
      geom = "crossbar", 
      width = 0.33, 
      color = "black", 
      fatten = 1.5
    ) + labs(x="FWHM (seconds)")+
    coord_cartesian(xlim = y_d0d4_val) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) 
  ggsave("FWHMs_IO_uPAC.pdf",width=3.5,height=4.5)
  
}

columns_specific_coverage <- function(dat, use_best_channel = F, filter_level=0.01, extra="IO"){
  
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(BestChannel_Qvalue<filter_level)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<filter_level)
  }
  
  dat <- dat %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))
  dat <- dat[dat$Repeatability ==F,]
  
  
  dat <- dat %>% dplyr::mutate(Label = ifelse(Type == "3-plexDIA","mTRAQ",Label))
  dat <- dat[dat$Type=="LF\n3-timePlex",]
  
  
  dat_counts <- dat %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    dplyr::ungroup()
  
  dat_counts$Column <- as.factor(as.numeric(dat_counts$Column)+1)
  colors = c("#edae49", "#d1495b","#00798c")
  ggplot(dat_counts, aes(x=Column, y=n,fill=Sample)) + geom_point(alpha=0.5, size=4,shape=21,height = 0,width=0.05)+
    facet_wrap(~Sample)+ theme_bw()+
    theme(legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_fill_manual(values=colors) + ylim(0,(max(dat_counts$n)+100))
  ggsave(paste0("points_tp_perColumnperSample",extra,".pdf"),width=3,height=3)
}

Missing_species_hist_quant <- function(dat, use_best_channel=F, filter_level=0.01){
  
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(BestChannel_Qvalue<filter_level)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<filter_level)
  }
  
  
  dat_lim <- dat %>% dplyr::add_count(untag_prec, channel)
  dat_lim <- dat_lim %>% dplyr::filter(n==1 & Proteotypic==T) #removes homologous peptides
  dat_lim <- dat_lim %>% dplyr::mutate(species=ifelse(grepl("HUMAN",ProteinName),"Human",
                                                      ifelse(grepl("YEAST",ProteinName),"Yeast","Arabidopsis")))
  
  dat_lim_d <- dcast(dat_lim, untag_prec+species~channel, value.var = "MS1_Int")
  dat_lim_d$t1_t2 <- dat_lim_d$`0`/dat_lim_d$`1`
  
  ggplot(dat_lim_d, aes(x=species,y=log2(t1_t2),color=species)) + geom_violin(alpha=0.1)
  
  
  dat_lim_counts <- dat_lim %>%
    dplyr::filter(MS1_Area > 1) %>%
    dplyr::group_by(species, channel) %>%
    dplyr::summarise(count = n(), .groups = 'drop')
  
  dat_lim_counts <- dat_lim_counts %>%
    dplyr::group_by(channel) %>%
    dplyr::mutate(y_pos = seq(max(log10(dat_lim$MS1_Int), na.rm = TRUE), by = -0.5, length.out = n()))
  
  ggplot(dat_lim, aes(x = log10(MS1_Int), fill = species)) +
    geom_histogram(alpha = 0.75, position = "identity",color="black",size=0.1) +
    facet_grid(~channel) +
    labs(x = "MS1 precursor intensity, log10", title = "H, HY, HY") +
    theme_bw() +
    scale_fill_manual(values=c("#7FDD6C","#68AABC","#E08A6C"))+
    geom_text(data = dat_lim_counts, aes(x = max(log10(dat_lim$MS1_Int), na.rm = TRUE) * 0.9, 
                                         y = y_pos, 
                                         label = paste(species, ": n =", count)),
              hjust = 1, inherit.aes = FALSE) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) #+ylim(0,5800)
  
  ggsave("precursor_counts_H_HY_HY.pdf",width=12,height=2.7)
  
  
  
  
  
  
}

extra_MS1_mzml_tp <- function(mzml_file, target_mz_first = 534.3383, z=2, 
                              mztol=0.0075, iso=14, rt_min=31.6, rt_max=40.6,
                              tp1=31.93026,
                              tp2=36.03937,
                              tp3=40.28072,
                              target_mz_1st = 435.2601,
                              ppm=10,
                              extra=""){
  
  mzdata <- openMSfile(mzml_file)
  hdr <- header(mzdata)
  ms1_scans <- hdr[hdr$msLevel == 1, ]
  
  mztol <- target_mz_1st*ppm/1000000
  
  if (nrow(ms1_scans) > 0) {
    ms1_intensity_data <- data.frame(RT = ms1_scans$retentionTime, TotalIntensity = NA)
    
    for (i in seq_len(nrow(ms1_scans))) {
      peaks_data <- peaks(mzdata, ms1_scans$acquisitionNum[i])  
      ms1_intensity_data$TotalIntensity[i] <- sum(peaks_data[, 2], na.rm = TRUE) # sum intensities per MS1
    }
    
    rect_data <- data.frame(
      ymin = -Inf, ymax = Inf,
      xmin = c(tp1 - 0.13, tp2 - 0.13, tp3 - 0.13),
      xmax = c(tp1 + 0.13, tp2 + 0.13, tp3 + 0.13)
    )
    
    ggplot(ms1_intensity_data, aes(y = (TotalIntensity), x = RT/60)) +  
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) +  
      geom_line(color = "black") +
      labs(
        y = "TIC",
        x = "Retention Time (minutes)") +
      theme_classic() +
      coord_flip()
    ggsave(paste0(extra,"_TIC.pdf"),width=3,height=3.6*1.5)
    
    
    ggplot(ms1_intensity_data, aes(y = (TotalIntensity), x = RT/60)) +  
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) + 
      geom_line(color = "black") +
      labs(
        y = "TIC",
        x = "Retention Time (minutes)") +
      theme_classic()
    ggsave(paste0(extra,"_TIC_xaxis.pdf"),height=2.4,width=3.6*2.4)
    
  }
  
  
  target_mz_values <- c(target_mz_1st)
  mz_tolerance <- target_mz_1st*ppm/1000000
  
  charge <- z
  
  compute_isotopes <- function(mono_mz, charge, num_isotopes = 1) {
    return((mono_mz) + (0:num_isotopes * 1.00335 / charge))
  }
  
  isotope_mz_list <- lapply(target_mz_values, compute_isotopes, charge = charge)
  
  # extract summed intensity
  extract_binned_signal <- function(scan, isotope_mz_set, mzdata, mz_tolerance) {
    peaks <- peaks(mzdata, scan)
    mz_values <- peaks[, 1]
    intensities <- peaks[, 2]
    
    summed_intensity <- sum(sapply(isotope_mz_set, function(mz_center) {
      sum(intensities[abs(mz_values - mz_center) <= mz_tolerance])
    }))
    
    return(summed_intensity)
  }
  
  binned_data <- data.frame(
    scan = ms1_scans$seqNum,
    retention_time = ms1_scans$retentionTime / 60  
  )
  
  for (i in seq_along(target_mz_values)) {
    target_mass <- target_mz_values[i]
    isotopes <- isotope_mz_list[[i]]
    
    summed_intensities <- sapply(ms1_scans$seqNum, function(scan) {
      extract_binned_signal(scan, isotopes, mzdata, mz_tolerance)
    })
    
    binned_data[[paste0("sum_", round(target_mass, 4))]] <- summed_intensities
  }
  
  binned_long <- melt(binned_data, id.vars = c("retention_time", "scan"),
                      variable.name = "mass_category", value.name = "Intensity")

  binned_data <- data.frame(
    scan = ms1_scans$seqNum,
    retention_time = ms1_scans$retentionTime / 60 
  )
  
  for (i in seq_along(target_mz_values)) {
    target_mass <- target_mz_values[i]
    isotopes <- isotope_mz_list[[i]]
    
    summed_intensities <- sapply(ms1_scans$seqNum, function(scan) {
      extract_binned_signal(scan, isotopes, mzdata, mz_tolerance)
    })
    
    binned_data[[paste0("sum_", round(target_mass, 4))]] <- summed_intensities
  }
  

  binned_long <- melt(binned_data, id.vars = c("retention_time", "scan"),
                      variable.name = "mass_category", value.name = "Intensity")
  
  
  
  binned_long <- binned_long %>% dplyr::mutate("timePlex" = ifelse(retention_time<tp1+0.33 & retention_time>tp1-0.33, "1",
                                                                   ifelse(retention_time<tp2+0.33 & retention_time>tp2-0.33, "2",
                                                                          ifelse(retention_time<tp3+0.33 & retention_time>tp3-0.33, "3","none")))) %>%
    dplyr::filter(timePlex!="none")
  
  
  binned_long$timePlex <- factor(binned_long$timePlex, levels = c("3", "2", "1"))
  
  binned_long <- binned_long %>%
    dplyr::mutate(interpolated = FALSE)
  
  
  smooth_intensity <- function(df) {
    df <- df %>% dplyr::arrange(retention_time)  
    if (nrow(df) > 3) {  
      spline_fit <- spline(df$retention_time, df$Intensity, n = 50)  
      df_interp <- data.frame(
        retention_time = spline_fit$x, 
        y_interpolated = spline_fit$y, 
        mass_category = unique(df$mass_category), 
        timePlex = unique(df$timePlex),
        interpolated = TRUE  
      )
      return(df_interp)
    } else {
      return(df %>% dplyr::mutate(y_interpolated = Intensity, interpolated = FALSE))  
    }
  }
  
  #smoothing per mass_category & celltype
  binned_long_smooth <- binned_long %>% 
    dplyr::group_by(mass_category, timePlex) %>%  
    group_split() %>% 
    lapply(smooth_intensity) %>% 
    bind_rows()
  
  binned_long <- binned_long %>% dplyr::rename(Intensity_real = Intensity)
  
  binned_long_smooth <- full_join(binned_long_smooth, binned_long, 
                                  by = c("retention_time", "mass_category", "timePlex", "interpolated")) %>%
    dplyr::mutate(y_interpolated = ifelse(is.na(y_interpolated), Intensity_real, y_interpolated))
  
  binned_long_smooth <- binned_long_smooth %>%
    
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4)  
  
  
  binned_long_smooth <- binned_long_smooth %>%
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4)  
  binned_long_smooth$y_interpolated[binned_long_smooth$y_interpolated<0] <- 0
  binned_long_smooth$timePlex <- factor(binned_long_smooth$timePlex, levels = c("1", "2", "3"))
  
  ggplot(binned_long_smooth, 
         aes(x = retention_time, y = y_interpolated, group = mass_category))+
    geom_hline(yintercept = 0, size = 1) +
    
    facet_grid(cols = vars(timePlex), scales = "free") +
    
    geom_point(data = binned_long_smooth %>% dplyr::filter(interpolated == FALSE), 
               aes(y = y_interpolated), color = "black", size = 2,alpha=0.9) +  
    
    geom_line(size = 1.3) +
    labs(x = "Retention Time (min)", y = "Intensity") +
    scale_alpha(range = c(0.4, 1)) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing.x = unit(1.3, "lines"),
          legend.position = "none")
  ggsave(paste0("Smoothed_Summed_Isotope_Signals_Bulk_TP_",extra,".pdf"), width = 5.5, height = 2.7)
  
}


jusTIC <- function(mzml_file, tp1=31.93026,
                   tp2=36.03937,
                   tp3=40.28072,
                   tp1_2=31.93026,
                   tp2_2=36.03937,
                   tp3_2=40.28072,
                   extra="",log_trans=F){
  
  mzdata <- openMSfile(mzml_file)
  hdr <- header(mzdata)
  ms1_scans <- hdr[hdr$msLevel == 1, ]
  
  
  
  if (nrow(ms1_scans) > 0) {
    ms1_intensity_data <- data.frame(RT = ms1_scans$retentionTime, TotalIntensity = NA)
    
    for (i in seq_len(nrow(ms1_scans))) {
      peaks_data <- peaks(mzdata, ms1_scans$acquisitionNum[i])  
      ms1_intensity_data$TotalIntensity[i] <- sum(peaks_data[, 2], na.rm = TRUE) 
    }
    
    rect_data <- data.frame(
      ymin = -Inf, ymax = Inf,  
      xmin = c(tp1 - 0.13, tp2 - 0.13, tp3 - 0.13, tp1_2 - 0.13, tp2_2 - 0.13, tp3_2 - 0.13),
      xmax = c(tp1 + 0.13, tp2 + 0.13, tp3 + 0.13, tp1_2 + 0.13, tp2_2 + 0.13, tp3_2 + 0.13)
    )
    
    if(log_trans==T){
      
      ggplot(ms1_intensity_data, aes(y = log10(TotalIntensity), x = RT/60)) +  
        geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) +  
        geom_line(color = "black") +
        labs(
          y = "Log10, TIC",
          x = "Retention Time (minutes)") +
        theme_classic() +
        coord_flip()
      ggsave(paste0(extra,"_TIC.pdf"),width=3,height=3.6*1.5)
      
      
      ggplot(ms1_intensity_data, aes(y = log10(TotalIntensity), x = RT/60)) +  
        geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) +  
        geom_line(color = "black") +
        labs(
          y = "Log10, TIC",
          x = "Retention Time (minutes)") +
        theme_classic()
      ggsave(paste0(extra,"_TIC_xaxis.pdf"),height=2.6,width=3.6*1.7)
      
    } else{
      ggplot(ms1_intensity_data, aes(y = (TotalIntensity), x = RT/60)) +  
        geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) +  
        geom_line(color = "black") +
        labs(
          y = "TIC",
          x = "Retention Time (minutes)") +
        theme_classic() +
        coord_flip()
      ggsave(paste0(extra,"_TIC.pdf"),width=3,height=3.6*1.5)
      
      
      ggplot(ms1_intensity_data, aes(y = (TotalIntensity), x = RT/60)) +  
        geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) +  
        geom_line(color = "black") +
        labs(
          y = "TIC",
          x = "Retention Time (minutes)") +
        theme_classic()
      ggsave(paste0(extra,"_TIC_xaxis.pdf"),height=2.4,width=3.6*3)
    }
    
    
    
  }
  
  
}

compare_coverage <- function(mTRAQ, LF){
  
  mTRAQ <- mTRAQ %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))
  
  LF <- LF %>%
    dplyr::mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
      timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
      timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
      TRUE ~ "LF\nNo-plex"
    ))
  
  mTRAQcounts <- mTRAQ %>% dplyr::filter(decoy==F & Qvalue <0.2 & BestChannel_Qvalue<0.01) %>%
    dplyr::add_count(Run) %>% dplyr::select("Run","n","Type") %>% dplyr::distinct(Run,.keep_all=T)
  
  LFcounts <- LF %>% dplyr::filter(decoy==F & Qvalue <0.2 & BestChannel_Qvalue<0.01) %>%
    dplyr::add_count(Run) %>% dplyr::select("Run","n") %>% dplyr::distinct(Run,.keep_all=T) %>% 
    dplyr::rename("LF_n"="n")
  
  joined <- mTRAQcounts %>% inner_join(LFcounts, by = "Run")
  
  joined$change <- 1-(joined$LF_n/joined$n)
  
  # Calculate mean and standard error
  summary_stats <- joined %>%
    dplyr::group_by(Type) %>%
    dplyr::summarise(mean_change = mean(change), 
              se = sd(change) / sqrt(n()))
  
  # Plot
  ggplot(joined, aes(x = Type)) + 
    geom_bar(data = summary_stats, aes(x = Type, y = mean_change), 
             stat = "identity", width = 0.6,fill="gray70",color="black") +
    geom_errorbar(data = summary_stats, aes(x = Type, ymin = mean_change - se, ymax = mean_change + se), 
                  width = 0.2) +
    geom_beeswarm(aes(x = Type,y=change),size = 3, alpha = 0.7, cex = 4) +  
    ylim(-0.01,0.17) + labs(y="Change in # precursors identified")+
    theme_classic()
  
  ggsave("ChangePrecursors_IDed_searchlib.pdf",width=3.4,height=3.4)
  
  
  
  LFcounts <- LF %>% dplyr::filter(decoy==F & Qvalue <0.2 & BestChannel_Qvalue<0.01) %>%
    dplyr::add_count(Run) %>% dplyr::select("Run","n","Type") %>% dplyr::distinct(Run,.keep_all=T) 
  
  LFcounts$lib <- "badlib"
  mTRAQcounts$lib <- "mTRAQ"
  
  both <- rbind(LFcounts, mTRAQcounts)
  
  
  # Calculate mean and standard error
  summary_stats_both <- both %>% dplyr::group_by(Type,lib) %>% dplyr::summarise(mean_n = mean(n), 
                                                                  se = sd(n) / sqrt(n()))
  
  ggplot(both, aes(x = Type, fill=lib)) + 
    geom_col(data = summary_stats_both, aes(x = Type, y = mean_n, fill=lib), 
             position = "dodge", width = 0.6,color="black",alpha=0.7) +
    geom_errorbar(data = summary_stats_both, aes(x = Type, ymin = mean_n - se, ymax = mean_n + se), 
                  position="dodge", width = 0.6) +
    geom_beeswarm(aes(x = Type,y=n,color=lib),size = 2, alpha = 0.5, cex = 2) +  
    ylim(0,145000) + labs(y="# precursors")+
    theme_classic() + 
    scale_fill_manual(values=c("#00BFC4","#F8766D"))+
    scale_color_manual(values=c("#00BFC4","#F8766D"))
  
  
  ggsave("Precursors_IDed_searchlib.pdf",width=4.5,height=3.4)
  
}

pD_weighted_PCA_timePlex <- function(dat){
  
  mat.sc.imp <- dat$imputed.BC
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
  plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")
  pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
  # Re map ...
  pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)
  pca.display$celltype<-no_neg$Celltype[match(pca.display$id, no_neg$run_chan)]
  pca.display$sample<-no_neg$Amount[match(pca.display$id, no_neg$run_chan)]
  
  pca.display$label<-no_neg$Label[match(pca.display$id, no_neg$run_chan)]
  pca.display$Well<-no_neg$V9[match(pca.display$id, no_neg$run_chan)]
  
  pca.display_lim <- pca.display %>% dplyr::select("PC1","Well","label","celltype")
  # PC's to display:
  PCx<-"PC1"
  PCy<-"PC2"
  PCz<-"PC3"
  PCza<-"PC4"
  PCzb<-"PC5"
  
  
  # Display
  pca.display$lab <- paste0("d",pca.display$label)
  pca.display$celltype <- factor(pca.display$celltype)
  pca.display2 <- pca.display %>% left_join(no_neg, by =c("id"="run_chan"))
  pca_lab<-ggplot(pca.display,aes(x =PC1, y = PC2, color =lab),  size = 3, alpha=0.5) +
    labs(color = "Label", shape="",x = paste0(PCx,"  (", round(percent_var[1],0),"%)"), y = paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
    scale_shape_manual(values = c(23,16,15)) + 
    scale_alpha_manual(values = c(0.66,0.66,0.66)) +
    font("ylab",size=20) +
    font("xlab",size=20) +
    font("xy.text", size=15) +
    geom_point(aes(fill=as.character(lab)),color="black",shape=21,size=3,alpha=0.5) +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.background = element_rect(color = "transparent", fill = "transparent")) +
    annotate("text", x=-0.072, y=-0.1, label=paste0(ncol(mat.sc.imp), " cells"), size=3) +
    annotate("text", x=-0.04, y=-0.12, label=paste0(dim(mat.sc.imp)[1], " protein groups"), size=3)+
    guides(size = "none",alpha = "none") 
  
  ggsave("SC_tp_pd_LABEL.pdf",width=4,height=4)
  
  
  pca_LC<-ggplot(pca.display2,aes(x =PC1, y = PC2)) +
    labs(color = "LC_batch", shape="",x = paste0(PCx,"  (", round(percent_var[1],0),"%)"), y = paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
    scale_shape_manual(values = c(23,16,15)) + 
    scale_alpha_manual(values = c(0.66,0.66,0.66)) +
    font("ylab",size=20) +
    font("xlab",size=20) +
    font("xy.text", size=15) +
    geom_point(aes(fill=as.character(LC_batch)),color="black",shape=21,size=3,alpha=0.5)+
    theme_classic() +
    theme(
      legend.position = "top",
      legend.background = element_rect(color = "transparent", fill = "transparent")) +
    annotate("text", x=-0.072, y=-0.1, label=paste0(ncol(mat.sc.imp), " cells"), size=3) +
    annotate("text", x=-0.04, y=-0.12, label=paste0(dim(mat.sc.imp)[1], " protein groups"), size=3)+
    guides(size = "none",alpha = "none") +
    labs(fill = "timePlex column",x="PC1",y = "PC2") 
  pca_LC
  
  ggsave("SC_tp_pd_LC_column.pdf",width=4,height=4)
  
  
  pca.display$sample_celltype <- paste0(pca.display$celltype,"_",pca.display$sample)
  pca<-ggplot(pca.display) +
    labs(shape="",x = paste0(PCx,"  (", round(percent_var[1],0),"%)"), y = paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
    font("ylab",size=20) +
    font("xlab",size=20) +
    font("xy.text", size=15) +
    geom_point(data = pca.display %>% dplyr::filter(!grepl("100x",sample_celltype)),
               aes(x =PC1, y = PC2, fill=celltype),alpha=0.6,colour="black",pch=21, size=3) +
    geom_point(data = pca.display %>% dplyr::filter(grepl("100x",sample_celltype)),
               aes(x =PC1, y = PC2, fill=celltype),alpha=0.6,colour="black",pch=23, size=7)+    
    scale_fill_manual(values=c("#e29578","#006d77"))+
    
    theme_classic() +
    theme(
      legend.position = "top",
      legend.background = element_rect(color = "transparent", fill = "transparent")) +
    annotate("text", x=-0.02, y=-0.21, label=paste0((ncol(mat.sc.imp)-12), " single cells"), size=3) +
    annotate("text", x=-0.02, y=-0.24, label=paste0(dim(mat.sc.imp)[1], " protein groups"), size=3)
  pca
  ggsave("SC_tp_pd.pdf",width=4,height=4)
  
  
  dat_fin <- list(pca = pca, pca_lab = pca_lab, pca_LC=pca_LC)
  return(dat_fin)
}


extra_MS1_mzml_v2 <- function(mzml_file, target_mz_first = 534.3383, z=2, 
                              mztol=0.0075, iso=14, rt_min=31.6, rt_max=40.6,
                              tp1=31.93026,
                              tp2=36.03937,
                              tp3=40.28072,
                              target_mz_1st = 534.3383,
                              target_mz_2nd = 536.3418,
                              target_mz_3rd = 538.3453,
                              target_mz_4th = NA,
                              target_mz_5th = NA,
                              target_mz_6th = NA,
                              target_mz_7th = NA,
                              target_mz_8th = NA,
                              target_mz_9th = NA,
                              pD_plex = 3,
                              ppm=10,
                              extra="",mycolor="#e8a792",SC=F){
  
  mzdata <- openMSfile(mzml_file)
  hdr <- header(mzdata)
  ms1_scans <- hdr[hdr$msLevel == 1, ]
  
  mztol <- target_mz_2nd*ppm/1000000
  
  if (nrow(ms1_scans) > 0) {
    ms1_intensity_data <- data.frame(RT = ms1_scans$retentionTime, TotalIntensity = NA)
    
    for (i in seq_len(nrow(ms1_scans))) {
      peaks_data <- peaks(mzdata, ms1_scans$acquisitionNum[i])  
      ms1_intensity_data$TotalIntensity[i] <- sum(peaks_data[, 2], na.rm = TRUE) 
    }
    
    rect_data <- data.frame(
      ymin = -Inf, ymax = Inf,
      xmin = c(tp1 - 0.13, tp2 - 0.13, tp3 - 0.13),
      xmax = c(tp1 + 0.13, tp2 + 0.13, tp3 + 0.13)
    )
    
    ggplot(ms1_intensity_data, aes(y = log10(TotalIntensity), x = RT/60)) +  
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) +  
      geom_line(color = "black") +
      labs(
        y = "Log10, TIC",
        x = "Retention Time (minutes)") +
      theme_classic() +
      coord_flip()
    ggsave(paste0(extra,"_TIC.pdf"),width=3,height=3.6*1.5)
    
    ggplot(ms1_intensity_data, aes(y = log10(TotalIntensity), x = RT/60)) +  
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) + 
      geom_line(color = "black") +
      labs(
        y = "Log10, TIC",
        x = "Retention Time (minutes)") +
      theme_classic()
    ggsave(paste0(extra,"_TIC_xaxis.pdf"),height=2.6,width=3.6*1.7)
    
    
    ggplot(ms1_intensity_data, aes(y = (TotalIntensity), x = RT/60)) +  
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) +  
      geom_line(color = "black") +
      labs(
        y = "TIC",
        x = "Retention Time (minutes)") +
      theme_classic() +
      coord_flip()
    ggsave(paste0(extra,"_TIC_noLog.pdf"),width=3,height=3.6*1.5)
    
    ggplot(ms1_intensity_data, aes(y = (TotalIntensity), x = RT/60)) +  
      geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = "mediumpurple1", alpha = 1, inherit.aes = FALSE) + 
      geom_line(color = "black") +
      labs(
        y = "TIC",
        x = "Retention Time (minutes)") +
      theme_classic()
    ggsave(paste0(extra,"_TIC_xaxis_noLog.pdf"),height=2.6,width=3.6*1.7)
    
    
  } 
  
  if(pD_plex==3){
    target_mz_values <- c(target_mz_1st, target_mz_2nd, target_mz_3rd)
  } else{
    target_mz_values <- c(target_mz_1st, target_mz_2nd, target_mz_3rd, target_mz_4th,
                          target_mz_5th,
                          target_mz_6th,
                          target_mz_7th,
                          target_mz_8th,
                          target_mz_9th)
    
  }
  
  
  mz_tolerance <- target_mz_2nd*ppm/1000000
  
  charge <- z  #Charge state
  
  compute_isotopes <- function(mono_mz, charge, num_isotopes = 1) {
    return((mono_mz) + (0:num_isotopes * 1.00335 / charge))
  }
  
  isotope_mz_list <- lapply(target_mz_values, compute_isotopes, charge = charge,num_isotopes = 1)
  
  extract_binned_signal <- function(scan, isotope_mz_set, mzdata, mz_tolerance) {
    peaks <- peaks(mzdata, scan)
    mz_values <- peaks[, 1]
    intensities <- peaks[, 2]
    
    summed_intensity <- sum(sapply(isotope_mz_set, function(mz_center) {
      sum(intensities[abs(mz_values - mz_center) <= mz_tolerance])
    }))
    
    return(summed_intensity)
  }
  
  binned_data <- data.frame(
    scan = ms1_scans$seqNum,
    retention_time = ms1_scans$retentionTime / 60  # Convert to minutes
  )
  
  for (i in seq_along(target_mz_values)) {
    target_mass <- target_mz_values[i]
    isotopes <- isotope_mz_list[[i]]
    
    summed_intensities <- sapply(ms1_scans$seqNum, function(scan) {
      extract_binned_signal(scan, isotopes, mzdata, mz_tolerance)
    })
    
    binned_data[[paste0("sum_", round(target_mass, 4))]] <- summed_intensities
  }

  binned_long <- melt(binned_data, id.vars = c("retention_time", "scan"),
                      variable.name = "mass_category", value.name = "Intensity")

  binned_long <- binned_long %>% dplyr::mutate("timePlex" = ifelse(retention_time<tp1+0.13 & retention_time>tp1-0.13, "1",
                                                                   ifelse(retention_time<tp2+0.13 & retention_time>tp2-0.13, "2",
                                                                          ifelse(retention_time<tp3+0.13 & retention_time>tp3-0.13, "3","none")))) %>%
    dplyr::filter(timePlex!="none")
  
  
  binned_long$timePlex <- factor(binned_long$timePlex, levels = c("3", "2", "1"))
  #binned_long
  
  if(pD_plex==3){
    df_map <- data.frame("mapper" = c(paste0("sum_",target_mz_1st,"_1"), paste0("sum_",target_mz_2nd,"_1"), paste0("sum_",target_mz_3rd,"_1"),
                                      paste0("sum_",target_mz_1st,"_2"), paste0("sum_",target_mz_2nd,"_2"), paste0("sum_",target_mz_3rd,"_2"),
                                      paste0("sum_",target_mz_1st,"_3"), paste0("sum_",target_mz_2nd,"_3"), paste0("sum_",target_mz_3rd,"_3")),
                         "celltype" = c("U937","K562","K562","U937","K562","U937","K562","U937","U937"))
   } else{
    df_map <- data.frame("mapper" = c(paste0("sum_",target_mz_1st,"_1"), paste0("sum_",target_mz_2nd,"_1"), paste0("sum_",target_mz_3rd,"_1"),
                                      paste0("sum_",target_mz_4th,"_1"), paste0("sum_",target_mz_5th,"_1"), paste0("sum_",target_mz_6th,"_1"),
                                      paste0("sum_",target_mz_7th,"_1"), paste0("sum_",target_mz_8th,"_1"), paste0("sum_",target_mz_9th,"_1"),
                                      
                                      paste0("sum_",target_mz_1st,"_2"), paste0("sum_",target_mz_2nd,"_2"), paste0("sum_",target_mz_3rd,"_2"),
                                      paste0("sum_",target_mz_4th,"_2"), paste0("sum_",target_mz_5th,"_2"), paste0("sum_",target_mz_6th,"_2"),
                                      paste0("sum_",target_mz_7th,"_2"), paste0("sum_",target_mz_8th,"_2"), paste0("sum_",target_mz_9th,"_2"),
                                      
                                      paste0("sum_",target_mz_1st,"_3"), paste0("sum_",target_mz_2nd,"_3"), paste0("sum_",target_mz_3rd,"_3"),
                                      paste0("sum_",target_mz_4th,"_3"), paste0("sum_",target_mz_5th,"_3"), paste0("sum_",target_mz_6th,"_3"),
                                      paste0("sum_",target_mz_7th,"_3"), paste0("sum_",target_mz_8th,"_3"), paste0("sum_",target_mz_9th,"_3")),
                         "celltype" = c("4x","1x","2x","3x","4x","3x","1x","2x","3x",
                                        "4x","1x","2x","3x","4x","3x","1x","2x","3x",
                                        "4x","1x","2x","3x","4x","3x","1x","2x","3x"))
   }
  
  binned_long$mapper <- paste0(binned_long$mass_category,"_",binned_long$timePlex)
  binned_long <- binned_long %>% left_join(df_map, by =c("mapper"="mapper"))
  
  binned_long <- binned_long %>%
    dplyr::mutate(interpolated = FALSE)  
  
  smooth_intensity <- function(df) {
    df <- df %>% dplyr::arrange(retention_time)  
    if (nrow(df) > 3) {  
      spline_fit <- spline(df$retention_time, df$Intensity, n = 50)  
      df_interp <- data.frame(
        retention_time = spline_fit$x, 
        y_interpolated = spline_fit$y, 
        mass_category = unique(df$mass_category), 
        celltype = unique(df$celltype),
        timePlex = unique(df$timePlex),
        interpolated = TRUE  
      )
      return(df_interp)
    } else {
      return(df %>% dplyr::mutate(y_interpolated = Intensity, interpolated = FALSE))  
    }
  }
  
  binned_long_smooth <- binned_long %>% 
    dplyr::group_by(mass_category, celltype, timePlex) %>%  
    group_split() %>% 
    lapply(smooth_intensity) %>% 
    bind_rows()
  
  binned_long <- binned_long %>% dplyr::rename(Intensity_real = Intensity)
  binned_long_smooth <- full_join(binned_long_smooth, binned_long, 
                                  by = c("retention_time", "mass_category", "celltype", "timePlex", "interpolated")) %>%
    dplyr::mutate(y_interpolated = ifelse(is.na(y_interpolated), Intensity_real, y_interpolated))
  
  binned_long_smooth <- binned_long_smooth %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4)  

  binned_long_smooth <- binned_long_smooth %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(norm_int = (y_interpolated - min(y_interpolated, na.rm = TRUE)) / 
             (max(y_interpolated, na.rm = TRUE) - min(y_interpolated, na.rm = TRUE)),
           alpha1 = norm_int * 0.6 + 0.4,  
           color_group = case_when(  
             celltype == "K562" ~ norm_int,   
             celltype == "U937" ~ -norm_int   
           ))  
  binned_long_smooth$y_interpolated[binned_long_smooth$y_interpolated<0] <- 0
  binned_long_smooth$timePlex <- factor(binned_long_smooth$timePlex, levels = c("1", "2", "3"))
  
  if(SC==T){
    ggplot(binned_long_smooth, 
           aes(x = retention_time, y = y_interpolated, group = mass_category, 
               color = celltype))+
      geom_hline(yintercept = 0, size = 1) +
      facet_grid(cols = vars(timePlex),rows=vars(mass_category), scales = "free") +
      geom_point(data = binned_long_smooth %>% dplyr::filter(interpolated == FALSE), 
                 aes(y = y_interpolated), color = "black", size = 2,alpha=0.9) +  
      geom_line(size = 1.3) +
      scale_color_manual(values=c("#e29578","#006d77"))+
      labs(x = "Retention Time (min)", y = "Intensity") +
      scale_alpha(range = c(0.4, 1)) +  
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 14),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text.y = element_text(size = 14),
            strip.text.x = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing.x = unit(1.3, "lines"),
            legend.position = "none") +
      ggsave(paste0("Smoothed_Summed_Isotope_Signals_",extra,".pdf"), width = 4, height = 5)
  }
  
  
  # 
  ggplot(binned_long_smooth,
         aes(x = retention_time, y = y_interpolated, group = mass_category))+#color = celltype))+
    geom_hline(yintercept = 0, size = 1) +
    facet_grid(cols = vars(timePlex),rows=vars(mass_category), scales = "free_x") +
    geom_point(data = binned_long_smooth %>% dplyr::filter(interpolated == FALSE),
               aes(y = y_interpolated), color = "black", size = 2,alpha=0.9) +
    geom_line(size = 1.3,color=mycolor) +
    #scale_color_manual(values=c("#ec008c", "#2bb673","#00aeef", "#fbb040"))+
    labs(x = "Retention Time (min)", y = "Intensity") +
    scale_alpha(range = c(0.4, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing.x = unit(1.3, "lines"),
          legend.position = "none")#,
  ggsave(paste0("Smoothed_Summed_Isotope_Signals_raw",extra,".pdf"), width = 4, height = 12)
  
  
  binned_long_smooth_tp_norm <- binned_long_smooth %>% dplyr::group_by(timePlex) %>%
    dplyr::mutate("tp_norm" = y_interpolated/max(y_interpolated,na.rm=T)) %>% 
    dplyr::group_by(timePlex,mass_category) %>%
    dplyr::mutate("AUC" = sum(tp_norm,na.rm=T)) %>% ungroup() %>%
    dplyr::group_by(timePlex) %>%
    dplyr::mutate(AUC_rel = AUC/min(AUC)) %>%
    ungroup()
  
  binned_long_smooth_tp_norm_dis <- binned_long_smooth_tp_norm %>% dplyr::distinct(timePlex, mass_category,.keep_all = T)
  
  
  ggplot(binned_long_smooth_tp_norm,
         aes(x = retention_time, y = tp_norm, group = mass_category,color = celltype))+
    geom_hline(yintercept = 0, size = 1) +
    facet_grid(cols = vars(timePlex),rows=vars(mass_category), scales = "free_x") +
    geom_point(data = binned_long_smooth_tp_norm %>% dplyr::filter(interpolated == FALSE),
               aes(y = tp_norm), color = "black", size = 2,alpha=0.9) +
    geom_line(size = 1.3, color=mycolor) +
    #scale_color_manual(values=c("#ec008c", "#2bb673","#00aeef", "#fbb040"))+
    labs(x = "Retention Time (min)", y = "Intensity") +
    scale_alpha(range = c(0.4, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing.x = unit(1.3, "lines"),
          legend.position = "none")#,
  
  ggsave(paste0("Smoothed_Summed_Isotope_Signals_normalized",extra,".pdf"), width = 4, height = 12)
    
}



cor_bulk_SC <- function(BC_MLFQ, meta_path){
  
  dat <- BC_MLFQ$unimputed.BC
  dat$norm_prot <- 2^dat$norm_prot
  m <- fread(meta_path)
  m$run_chan <- paste0(m$Run,"_",(m$Column-1),"_",m$Label)
  dat <- dat %>% left_join(m,by =c("run_chan" = "run_chan"))
  dat_consensus <- dat %>% 
    dplyr::group_by(Sample, Amount, Genes) %>% dplyr::add_count(name="total_possible") %>% dplyr::select(-c("Rep")) %>%
    na.omit() %>% dplyr::add_count(name="total_quantified") %>% dplyr::mutate("frac_quant" = total_quantified/total_possible) %>%
    dplyr::filter(frac_quant>0.2) %>%
    dplyr::summarize("mean_prot" = median(norm_prot,na.rm=T))
  dat_consensus_d <- dcast(dat_consensus, Genes~Sample+Amount, value.var = "mean_prot")
  dat_consensus_d$bulk_KU <- log2(dat_consensus_d$K562_100x)-log2(dat_consensus_d$U937_100x)
  dat_consensus_d$SC_KU <- log2(dat_consensus_d$K562_SC)-log2(dat_consensus_d$U937_SC)
  dat_consensus_d <- dat_consensus_d %>% na.omit()
  fit <- tls_fit(dat_consensus_d$bulk_KU, dat_consensus_d$SC_KU)
  
  highlight_genes <- subset(dat_consensus_d, Genes %in% c("HMGA1_HUMAN", "FABP5_HUMAN"))
  
  ggplot(dat_consensus_d, aes(x = bulk_KU, y = SC_KU)) + 
    geom_point(alpha = 1, shape = 21, size = 2, color = "black") + 
    geom_pointdensity(alpha = 0.5, size = 2) + 
    geom_text_repel(data = highlight_genes, 
                    aes(x = bulk_KU, y = SC_KU, label = Genes),
                    size = 3.5, fontface = "bold", color = "red",
                    box.padding = 0.5, point.padding = 0.25,
                    segment.color = "red", segment.size = 0.4, segment.alpha = 0.8,
                    max.overlaps = Inf) +
    
    labs(x = "100x bulk, LF 3-timePlex\nLog2, K562 / U937",
         y = "Single cell, 3-plexDIA & 3-timePlex\nLog2, K562 / U937",
         title = paste0(nrow(dat_consensus_d), " proteins"),
         subtitle = paste0("Pearson cor = ", round(cor(dat_consensus_d$bulk_KU, dat_consensus_d$SC_KU, method = "pearson"), 2))) +
    
    scale_color_viridis() +
    theme_classic() +
    theme(legend.position = "none")
  
  
  ggsave("SC_bulk_cor_prots.pdf",width=3.2,height=3.6)
  
}

tls_fit <- function(x, y) {
  fit <- prcomp(~ x + y)
  intercept <- fit$center[2] - fit$rotation[2,1] / fit$rotation[1,1] * fit$center[1]
  slope <- fit$rotation[2,1] / fit$rotation[1,1]
  list(intercept = intercept, slope = slope)
}

pD_batchCorrect_labs_LC <- function(dat, meta_fpath){
  sc.imp <- as.matrix(dat$sc.imp)
  sc.imp<-cr_norm_log((sc.imp))
  m <- fread(meta_fpath)
  m$BioRep <- "1"
  m$run_chan <- paste0(m$Run,"_",(m$Column-1),"_",m$Label)
  m$Celltype <- m$Sample
  m$LC_batch <- m$Column
  no_neg <- m[!grepl("neg|Carrier",m$Celltype),]
  no_neg$lab_set <- paste0(no_neg$Label,no_neg$LC_batch)
  uni_BR <- unique(no_neg$BioRep)
  temp_list <- list()
  
  for(i in 1:length(uni_BR)){
    no_neg_temp <- no_neg[no_neg$BioRep==uni_BR[i],]
    temp <- sc.imp[,colnames(sc.imp)%in%no_neg_temp$run_chan]
    batch.covs <-no_neg_temp$lab_set[match(colnames(temp), no_neg_temp$run_chan)]
    lab_set_counts <- table(batch.covs)  # Count the number of cells in each lab_set
    valid_lab_sets <- names(lab_set_counts[lab_set_counts >= 3]) # Filter lab_set groups that have 3 or more cells
    batch.covs <- batch.covs[batch.covs%in%valid_lab_sets]
    no_neg_temp <- no_neg_temp[no_neg_temp$lab_set%in%valid_lab_sets]
    temp <- temp[,colnames(temp)%in%no_neg_temp$run_chan] 
    temp<-cr_norm_log((temp))
    batch.covs <-no_neg_temp$lab_set[match(colnames(temp), no_neg_temp$run_chan)]
    matrix.sc.batch_temp <- ComBat(temp, batch=batch.covs)#)
    matrix.sc.batch_temp<-cr_norm_log((matrix.sc.batch_temp))
    temp_list[[i]] <- matrix.sc.batch_temp
  }
  matrix.sc.batch <- do.call("cbind",temp_list)
  imputed.BC<-cr_norm_log((matrix.sc.batch))
  imputed.BC_updated<-imputed.BC
  
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

join_post_imp <- function(SC, bulk){
  SC_pre_imp <- SC$pre_impute
  Bulk_pre_imp <- bulk$pre_impute
  pre_impute <- merge(SC_pre_imp, Bulk_pre_imp, by = "row.names", all = FALSE)
  rownames(pre_impute) <- pre_impute$Row.names
  pre_impute$Row.names <- NULL
  
  SC_imp <- SC$sc.imp
  Bulk_imp <- bulk$sc.imp
  sc.imp <- merge(SC_imp, Bulk_imp, by = "row.names", all = FALSE)
  rownames(sc.imp) <- sc.imp$Row.names
  sc.imp$Row.names <- NULL
  
  joined <- list(pre_impute = as.matrix(pre_impute), sc.imp=as.matrix(sc.imp))
  return(joined)
}

pD_impute_log2 <- function(df, missing.prot.frac=0.7, missing.cell.frac=0.95, k=10){
  dat_fin1 <- df %>% dplyr::select("Genes", "run_chan", "norm_prot")
  dat_fin1$Genes <- as.factor(dat_fin1$Genes)
  dat_fin1$norm_prot <- as.numeric(dat_fin1$norm_prot)
  dat_fin1$run_chan <- gsub('\\.', '-', dat_fin1$run_chan)
  dat_fin1$run_chan <- as.factor(dat_fin1$run_chan)
  dat_fin2<-dcast(dat_fin1, Genes~run_chan, value.var = "norm_prot", fill=NA)
  dat3<-as.matrix(dat_fin2[,-1]); row.names(dat3)<-dat_fin2[,1]
  
  ## Impute single celldata
  pre_impute <-(filt.mat.rc(dat3,missing.prot.frac,missing.cell.frac))
  sc.imp<-cr_norm_log(pre_impute)
  sc.imp <- hknn_weighted_EucNorm(sc.imp, k)
  sc.imp<-cr_norm_log((sc.imp))
  
  dat <- list(pre_impute = pre_impute, sc.imp = sc.imp)
  
  return(dat)
}

pD_prot_norm <- function(df, Quant = "prot_norm_adj"){
  
  df$Quantnew <- df[[Quant]]
  df$Genes <- df$prot
  dat_fin <- df %>% dplyr::group_by(run_chan) %>%
    dplyr::mutate("norm_prot" = Quantnew - median(Quantnew, na.rm=T)) %>%
    ungroup() %>% dplyr::group_by(Genes) %>%
    dplyr::mutate("norm_prot" = norm_prot - mean(norm_prot, na.rm=T)) %>% ungroup()
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


pD_MsEmpire_run_specific_jmod <- function(dat_fpaths_msEmp){
  dat <- dat_fpaths_msEmp$dat_fpaths_fin
  mappings <- dat_fpaths_msEmp$mappings_fpaths_fin

  BRs <- dat_fpaths_msEmp$BRs
  empty_list <- list()
  for(i in 1:length(mappings)){
    data <- msEmpiRe::read.standard(dat[i], mappings[i],
                                    prot.id.generator=function(pep) unlist(strsplit(pep, "\\."))[1],
                                    signal_pattern="c.*rep.*")
    
    #extract the first two conditions
    conditions <- extract_conditions(data)
    conditions_temp <- conditions[, c(1,2)]
    
    #removing peptides that are detected in less than 2 samples per condition
    data_temp <- msEmpiRe::filter_detection_rate(data, condition=conditions_temp, rate=2)
    
    system.time(data_temp <- msEmpiRe::normalize(data_temp))
    set.seed(1234)
    #main analysis, system.time is just to measure the time of processing
    system.time(result <- de.ana(data_temp))
    result$BioRep <- paste0(BRs[i])
    result$Cond <- paste0(dat[i])
    empty_list[[i]] <- data.frame(result)
  }
  df <- do.call("rbind", empty_list)
  df$Gene <- df$prot.id
  write.table(df, file = "KU_msEmpiRe.tsv", row.names = FALSE, sep = "\t")
  
  return(df)
}

pD_MsEmpire_prep_jmod <- function(Nuc_dat, BR=c(1,2,3,6),conds=c("NT","min_10","min_30","min_60")){
  Nuc_dat$BioRep <- "1"
  dat_fpaths_fin <- as.vector(NULL)
  mappings_fpaths_fin <- as.vector(NULL)
  BRs <- as.vector(NULL)
  Nuc_dat$id <- paste0(Nuc_dat$prot,".",Nuc_dat$seqcharge)
  Nuc_dat$Celltype <- Nuc_dat$Sample
  for(i in 1:length(BR)){
    Nuc_temp <- Nuc_dat[Nuc_dat$BioRep==i,]
    for(k in 2:length(conds)){
      conds_temp <- conds[c(1,k)]
      Nuc <- Nuc_temp[Nuc_temp$Celltype%in%conds_temp,]
      #### join with meta data:
      Nuc <- Nuc %>% dplyr::mutate("cond" = ifelse(grepl("K562", Celltype), "c1","c2"))
      Nuc$sample <- paste0(Nuc$cond, ".rep", Nuc$Rep, Nuc$BioRep, Nuc$time_channel)
      print(unique(Nuc$sample))
      Nuc <- Nuc %>% dplyr::select("id","MS1_Int","sample")
      Nuc_d <- reshape2::dcast(Nuc,id~sample, value.var = "MS1_Int")
      Nuc_d[is.na(Nuc_d)] <- 0
      mappings <- Nuc  %>% dplyr::distinct(sample) %>% dplyr::mutate("condition" = ifelse(grepl("c1", sample), "0", "1"))
      time <- Sys.time()
      dat_fpaths <- paste0("bulk_BR",BR[i],"_cond",conds[k],time,".tsv")
      mappings_fpath <- paste0("bulk_BR_mapping",BR[i],"_cond",conds[k],time,".tsv")
      write.table(Nuc_d,file=dat_fpaths, row.names = FALSE,sep="\t")
      write.table(mappings,file=mappings_fpath, row.names = FALSE,sep="\t")
      dat_fpaths_fin <- c(dat_fpaths,dat_fpaths_fin)
      mappings_fpaths_fin <- c(mappings_fpath,mappings_fpaths_fin)
      BRs <- c(paste0(BR[i]),BRs)
    }
  }
  fin <- list(dat_fpaths_fin=dat_fpaths_fin, mappings_fpaths_fin=mappings_fpaths_fin,BRs=BRs)
  return(fin)
}

plot_counts <- function(sc, sc_ms1_filt, filt=F){
  uniruns <- sc_ms1_filt %>% dplyr::distinct(run_chan, .keep_all=T)
  
  if(filt==T){
    prec <- sc$prec_counts %>% dplyr::filter(!is.na(Sample) & run_chan%in%uniruns$run_chan)
    prec$type <- "precursors"
    prot <- sc$prot_counts %>% dplyr::filter(!is.na(Sample) & run_chan%in%uniruns$run_chan)
    prot$type <- "proteins"
  } else{
    prec <- sc$prec_counts %>% dplyr::filter(!is.na(Sample))
    prec$type <- "precursors"
    prot <- sc$prot_counts %>% dplyr::filter(!is.na(Sample))
    prot$type <- "proteins"
  }
  
  all <- rbind.fill(prot,prec)
  
  sum <- all %>% dplyr::group_by(type,Sample) %>% dplyr::summarize(med_n=median(n,na.rm=T))
  sum <- all %>% dplyr::group_by(type) %>% dplyr::summarize(max_n=max(n,na.rm=T))
  
  ggplot(all %>% dplyr::filter(!is.na(Sample)), aes(x=Sample,y=n,fill=Sample)) + geom_boxplot(alpha=0.75,outliers = F) + 
    geom_beeswarm(alpha=0.5) +
    facet_grid(~type) + theme_bw() + 
    theme(          panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    legend.position="none")+
    scale_fill_manual(values=c("#e29578","#006d77")) + ylim(0,max(all$n))
  ggsave("counts_prec_prot.pdf",width=3,height=4)
  
}

sc_quality_control <- function(dat, min_counts = 250){
  
  dat <- dat %>% dplyr::group_by(run_chan) %>%
    dplyr::add_count() %>% dplyr::filter(n>250) %>%
    dplyr::ungroup() %>% dplyr::select(-c("n"))
  
  return(dat)
  
}

pD_ms1_quant <- function(df, adjacent_scans = 1,dontapplyLF=T) {
  df_lim <- df %>% dplyr::filter(grepl("mTRAQ", seq)) %>% dplyr::filter(BestChannel_Qvalue < 0.01 & Qvalue<0.2)
  df_lim2 <- df_lim %>%
    dplyr::mutate(
      MS1_fitted = pmap_chr(list(plexfittrace_ps_all, plexfittrace_all), function(ps, tr) {
        extract_highest_pearson_vec(ps, tr, adjacent_scans)
      }),
      MS1_fitted_sum = map_dbl(MS1_fitted, sum_fitted_values)
    )
  
  # for LF get the apex and adj scans.
  df_complement <- df %>% dplyr::filter(!(grepl("mTRAQ", seq) & BestChannel_Qvalue < 0.01  & Qvalue<0.2))
  if(dontapplyLF==F){
    df_complement <- df_complement %>%
      dplyr::mutate(
        MS1_fitted = pmap_chr(list(all_ms1_iso0vals, all_ms1_iso0vals), function(ps, tr) {
          extract_highest_pearson_vec(ps, tr, adjacent_scans)
        }),
        MS1_fitted_sum = map_dbl(MS1_fitted, sum_fitted_values)
      )
  }
  all <- bind_rows(df_complement, df_lim2) %>%
    dplyr::mutate(MS1_bestfit = ifelse(is.na(MS1_fitted_sum), plexfitMS1, MS1_fitted_sum))
  
  return(all)
}

compare_carryover <- function(dat_all_carryover){
  dat <- dat_all_carryover %>% dplyr::filter(BestChannel_Qvalue<0.01 & Qvalue<0.2 & decoy==F)
  dat <- dat[grepl("JD0407",dat$Run),]
  dat_d <- dcast(dat, untag_prec~Run+time_channel, value.var = "MS1_Area")
  dat_d$rat <- log10(dat_d$JD0407_1/dat_d$JD0407_0)
  10^median(dat_d$rat,na.rm=T)
  
  dat_d <- dat_d %>% dplyr::mutate(rat=ifelse(is.na(JD0407_0), log10(1000000),
                                              ifelse(is.na(JD0407_1), log10(1/1000000), rat)) )
  
  dat_d <- dat_d %>% dplyr::mutate(presence=ifelse(is.na(JD0407_0), "Unintended time channel only",
                                                   ifelse(is.na(JD0407_1), "Intended time channel only", "Both time channels")) )
  
  ggplot(dat_d, aes(x=rat, fill=presence)) + geom_histogram(alpha=0.5, color="black",size=0.2) + 
    theme_classic() + labs(x="Log10, Unintended time channel /  Intended time channel")
  ggsave("timechannel_carryover.pdf",width=7,height=6)
}

pD_ms1_quantINT <- function(df, adjacent_scans = 1) {
  df_lim <- df %>% dplyr::filter(grepl("mTRAQ", seq)) %>% dplyr::filter(BestChannel_Qvalue < 0.01)
  df_lim2 <- df_lim %>%
    dplyr::mutate(
      MS1_fitted = pmap_chr(list(all_ms1_iso0vals, all_ms1_iso0vals), function(ps, tr) {
        extract_highest_pearson_vec(ps, tr, adjacent_scans)
      }),
      MS1_fitted_sum = map_dbl(MS1_fitted, sum_fitted_values)
    )
  
  df_complement <- df %>% dplyr::filter(!(grepl("mTRAQ", seq) & BestChannel_Qvalue < 0.01))
  
  all <- bind_rows(df_complement, df_lim2) %>%
    dplyr::mutate(MS1_bestfit = ifelse(is.na(MS1_fitted_sum), plexfitMS1, MS1_fitted_sum))
  
  return(all)
}

throughput_plot <- function(){
  
  pD_tp <- data.frame(pd = seq(1:7),
                      tp = seq(1:7))
  
  pD_tp$tp_pd <-pD_tp$pd*pD_tp$tp
  pD_tp$plex <- as.factor(as.character(pD_tp$pd))
  pD_tp_m <- melt(pD_tp)
  pD_tp_m$plex <- as.numeric(as.character(pD_tp_m$plex))
  ggplot(pD_tp_m, aes(x=plex, y=value, color=variable,group=variable))+geom_line(size=2)  +
    geom_point(color="black")+ theme_classic() +
    labs(x="Plex", y="Throughput (samples / run)") + 
    theme(legend.position="none")
  ggsave("Throughput_samplesPerRun.pdf",width=2.5,height=3.9)
  
  
}

compute_prec_prot_qvalues_v2 <- function(df, toPlot=F){
  
  df_filt <- df %>% dplyr::filter(BestChannel_Qvalue<0.01 & decoy==F & Qvalue<0.2)
  df_filt_count <- df_filt %>%
    dplyr::group_by(run_chan,channel) %>% dplyr::add_count() %>% dplyr::ungroup() %>% dplyr::distinct(run_chan,channel,.keep_all=T)
  
  df_filt <- df_filt %>% dplyr::group_by(prot) %>% dplyr::mutate(min_prot_qval = min(Protein_Qvalue)) %>% dplyr::ungroup()
  df_filt_prot_count <- df_filt %>% dplyr::filter(min_prot_qval <0.01) %>% 
    dplyr::distinct(run_chan,channel,prot, .keep_all=T) %>%
    dplyr::group_by(run_chan,channel) %>% dplyr::add_count() %>% dplyr::ungroup() %>% dplyr::distinct(run_chan, channel, .keep_all=T)
  
  med_prec <- median(df_filt_count$n)
  med_prot <- median(df_filt_prot_count$n)
  
  if(toPlot==T){
    
    ggplot(df_filt_prot_count, aes(x=n)) + geom_histogram(color="black",fill="grey") + theme_classic() +
      labs(x="Proteins", y="Single cells",
           subtitle = paste0("Median: ",med_prot," proteins/cell"))
    ggsave("Proteins_per_SC.pdf",width=3.8,height=2.7)
    
    ggplot(df_filt_count, aes(x=n)) + geom_histogram(color="black",fill="grey") + theme_classic() +
      labs(x="Precursors", y="Single cells",
           subtitle = paste0("Median: ",med_prec," precursors/cell"))
    ggsave("Precursors_per_SC.pdf",width=3.8,height=2.7)
    
  }
  
  
  df_filt <- df_filt %>% dplyr::filter(BestChannel_Qvalue <0.01 & min_prot_qval <0.01 & Qvalue<0.2)
  
  fin_list <- list(df = df_filt, prec_counts = df_filt_count, prot_counts = df_filt_prot_count)
  return(fin_list)
  
}

'%!in%' <- function(x,y)!('%in%'(x,y))


gaussian_model <- function(x, A, mu, sigma) {
  A * exp(-((x - mu)^2) / (2 * sigma^2))
}

compute_fwhm <- function(seq_str) {
  intensity <- as.numeric(unlist(strsplit(seq_str, ";")))
  scans <- seq_along(intensity)
  
  if (sum(intensity > 0) < 3) return(NA_real_)
  
  df <- data.frame(x = scans, y = intensity)
  
  tryCatch({
    start_vals <- list(
      A = max(df$y),
      mu = df$x[which.max(df$y)],
      sigma = length(df$y) / 5
    )
    
    fit <- nls(y ~ gaussian_model(x, A, mu, sigma), data = df, start = start_vals, control = list(warnOnly = TRUE))
    
    coefs <- coef(fit)
    A <- coefs["A"]
    mu <- coefs["mu"]
    sigma <- coefs["sigma"]
    
    x_fine <- seq(min(df$x), max(df$x), by = 0.01)
    y_fine <- gaussian_model(x_fine, A, mu, sigma)
    #plot(x_fine,y_fine)
    #plot(scans,intensity)
    half_max <- A / 2
    
    above_half <- which(y_fine >= half_max)
    if (length(above_half) < 2) return(NA_real_)
    
    left_idx <- min(above_half)
    right_idx <- max(above_half)
    
    fwhm_precise <- x_fine[right_idx] - x_fine[left_idx]
    return(fwhm_precise)
    
  }, error = function(e) {
    return(NA_real_)
  })
}


#from HS SCoPE2
filt.mat.rc<-function(mat, pct.r,pct.c){
  
  kr<-c()
  for(k in 1:nrow(mat)){
    
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)
    
  }
  
  mat<-mat[kr,]
  
  kc<-c()
  for(k in 1:ncol(mat)){
    
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)
    
  }
  
  mat<-mat[,kc]
  
  return(mat)
  
}

hknn_weighted_EucNorm <-function(dat, k){
  
  #create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  #similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  counts <- pairwiseCount(dat)
  dist.mat <- dist.mat/counts
  
  cnames<-colnames(dist.mat)
  
  for(X in cnames){
    
    distances<-dist.mat[, X]
    
    # reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # reorder the data matrix columns, smallest distance to largest from the column of interest
    dat.reordered<-dat[ , names(distances.ordered ) ]
    vec<-dat[, X]
    na.index<-which( is.na(vec) )
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      distances.ordered_mat <- as.matrix(distances.ordered)
      
      #from the samples that are ordered in terms of their similarity, select the ones which have a non-NA value
      distances.ordered_mat <- distances.ordered_mat[rownames(distances.ordered_mat)%in%closest.columns,] 
      closest_vec <- distances.ordered_mat[1:(k)] #select the first 'k' samples with their euclidian distances
      closest_vec_inv <- 1/closest_vec #inverted euclidian distances (small distances = more similar = should have greater weight)
      total <- sum(closest_vec_inv) #sum the inverted euclidian distances
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the weighted average of the same row in k of the most similar columns
        vec[i]<-sum(closest_vec_inv*as.matrix(dat[ i, closest.columns[1:k] ] ))/total
        
      }
      
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        j <- length(closest.columns)
        closest_vec <- distances.ordered_mat[1:(j)]
        closest_vec_inv <- 1/closest_vec
        total <- sum(closest_vec_inv)
        # Replace NA in column X with the weighted avg of the same row in all of the most similar columns
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

quant <- function(dat, numerator = "A", denominator = "B", exp_FC = 2,
                  title_actual = "3-plexDIA & 3-timePlex: 8 ng HY (uPAC)",alphaScaler=1){
  quant_col <- paste0(numerator,denominator)
  dat$FCs <- dat[[quant_col]]
  dat$denominator <- dat[[denominator]]
  dat <- dat %>% filter(!is.na(FCs) & is.finite(FCs) ) %>% filter(FCs!=0)
  H_med <- dat %>% filter(grepl("Human",species))
  H_med <- median(H_med$FCs)
  dat$FCs_corrected <- dat$FCs-H_med
  
  box <- ggplot(dat, aes(x=species,y=FCs_corrected,fill=species)) + geom_boxplot(alpha=0.85,outlier.alpha=0.01*alphaScaler) + 
    geom_hline(yintercept = exp_FC,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_classic()+
    scale_fill_manual(values=c("#68aabc","#e08a6c")) +
    labs(y=paste0("Log2, ",numerator," / ",denominator), x ="Species", 
         subtitle = ""
    ) +
    theme(legend.position="none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) + coord_cartesian(ylim = c(-8,8))
  
  
  
  
  dense <- ggplot(dat, aes(fill=species,y=FCs_corrected)) + geom_density(alpha=0.75) + 
    geom_hline(yintercept = exp_FC,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_classic()+
    scale_fill_manual(values=c("#68aabc","#e08a6c")) +
    labs(
      subtitle = ""
    ) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank()) + coord_cartesian(ylim = c(-8,8))
  
  dat <- dat %>% add_count(species)
  dat <- dat %>% dplyr::mutate("Alpha"=0.5*(1-n/nrow(.)))
  dat <- dat %>% dplyr::mutate("Alpha" = ifelse(grepl("Human",species),0.05,0.2))
  
  scatter <- ggplot(dat, aes(x=denominator, y=FCs_corrected, color=species, alpha=Alpha*alphaScaler)) + geom_point() + 
    geom_hline(yintercept = exp_FC,linetype="dashed",color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + theme_classic()+
    scale_color_manual(values=c("#68aabc","#e08a6c")) +
    scale_alpha_identity() +
    labs(#title = title_actual, 
      subtitle = paste0("n = ",format(nrow(dat),big.mark=",",scientific=FALSE)," precursor ratios"),
      y=paste0("Log2, ",numerator," / ",denominator),
      x=paste0("Log2, ",denominator))+
    theme(legend.position="none") + coord_cartesian(ylim = c(-8,8))
  
  quantplot <- ggarrange(scatter, box, dense,
                         ncol = 3, nrow = 1,
                         widths = c(0.8, 0.2, 0.33))
  
  quantplot <- annotate_figure(quantplot, top = text_grob(title_actual, 
                                                          color = "black", face = "bold", size = 14))
  return(quantplot)
}

tp_coverage_3x9 <- function(dat,extra="", use_best_channel = F, filter_level=0.01, plex3x9=F){
  
  if(use_best_channel==T){
    dat <- dat %>% dplyr::filter(BestChannel_Qvalue<filter_level) %>% dplyr::filter(Qvalue < 0.2)
  } else{
    dat <- dat %>% dplyr::filter(Qvalue<filter_level)
  }
  
  if(plex3x9==T){
    dat <- dat %>%
      dplyr::mutate(Type = case_when(
        timePlex == TRUE & Label != "LF" ~ "9-plexDIA\n3-timePlex",
        timePlex == FALSE & Label != "LF" ~ "9-plexDIA",
        TRUE ~ "LF\nNo-plex"
      ))
  } else{
    dat <- dat %>%
      dplyr::mutate(Type = case_when(
        timePlex == TRUE & Label == "LF" ~ "LF\n3-timePlex",
        timePlex == TRUE & Label != "LF" ~ "3-plexDIA\n3-timePlex",
        timePlex == FALSE & Label != "LF" ~ "3-plexDIA",
        TRUE ~ "LF\nNo-plex"
      ))
  }
  
  dat <- dat[dat$Repeatability ==F,]
  dat_keep <- dat
  
  
  
  dat$Label <- paste0(dat$Label,"_",dat$time_channel)
  
  dat_counts <- dat %>%
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>%
    dplyr::ungroup()
  

  if(plex3x9==T){
    dat_counts$Sample <- factor(dat_counts$Sample, levels = c("J", "K", "L","M"))
  } else{
    dat_counts$Sample <- factor(dat_counts$Sample, levels = c("A", "B", "C"))
  }
  
  
  # Compute mean and standard error
  dat_summary <- dat_counts %>%
    dplyr::group_by(Sample, Label, Amount, Type) %>%
    dplyr::summarise(
      mean_count = mean(n, na.rm = TRUE),
      se = sd(n, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% dplyr::ungroup()
  
  dat_summary <- dat_summary %>%
    dplyr::arrange(desc(Sample), Type) %>%
    
    dplyr::group_by(Type) %>%
    dplyr::mutate(
      stacked_mean = cumsum(mean_count)) %>% dplyr::ungroup() %>%
    dplyr::group_by(Label, Sample, Type) %>%
    dplyr::mutate(stacked_mean = ifelse(Type=="LF\nNo-plex", mean_count,stacked_mean))%>%
    dplyr::mutate(
      er_max = stacked_mean + se,
      er_min = stacked_mean - se
    ) %>%
    dplyr::ungroup()
  
  df_text <- dat_summary %>%
    dplyr::group_by(Type, Sample,Label) %>%
    dplyr::summarise(
      text_value = round(max(er_max) - max(se)),
      text_pos = max(er_max) + 2000
    )
  
  my_colors <- c("#ec008c", "#2bb673","#00aeef", "#fbb040")
  
  ggplot(dat_summary, aes(x = Type, y = mean_count, fill = Sample)) +
    scale_x_discrete(limits = c("LF\nNo-plex","9-plexDIA", "9-plexDIA\n3-timePlex")) +
    
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "9-plexDIA\n3-timePlex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +

    geom_bar(data = dat_summary %>% dplyr::filter(Type == "9-plexDIA"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 1.17),
             color = "black", width = 1, alpha = 1) +
 
    scale_y_continuous(labels = comma) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("9-plexDIA\n3-timePlex", "9-plexDIA")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.18, size = 0.4) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.5, size = 0.4,
                  position = position_dodge(1.17)) +
    
    scale_fill_manual(values = my_colors) +
    theme_classic() +
    
    labs(y = "Precursors data points / Run", x = "", fill = "Samples") +
    
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = c(0.16, 0.9),
          legend.direction = "horizontal",
          legend.spacing.x = unit(1.0, 'line')) +
    
    guides(fill = guide_legend(override.aes = list(alpha = 0.7),
                               title = "Samples",
                               label.position = "bottom",
                               title.position = "top", title.vjust = 1, title.hjust = 0.5)) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("9-plexDIA\n3-timePlex", "9-plexDIA")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3) +
    
    geom_text(data = df_text %>% dplyr::filter(Type == "LF\nNo-plex"), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3,
              width = 0.5, size = 0.4,
              position = position_dodge(1.17))
    
  
  ggsave(paste0(extra,"_","Prec_coverage_","BestChannel_",use_best_channel,".pdf"),width=5.5,height=5.5)
  
  
  
  ##################

  dat_counts <- dat %>% dplyr::filter(Protein_Qvalue<filter_level) %>% dplyr::distinct(Sample, Run, Column, Amount, Label, Type, prot_original,.keep_all=T) %>% 
    dplyr::group_by(Sample, Run, Column, Amount, Label, Type) %>%
    dplyr::count() %>% dplyr::ungroup()
  
  dat_summary <- dat_counts %>%
    dplyr::group_by(Sample, Label, Amount, Type) %>%
    dplyr::summarise(
      mean_count = mean(n, na.rm = TRUE),
      se = sd(n, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>% dplyr::ungroup()
  
  dat_summary <- dat_summary %>%
    dplyr::arrange(desc(Sample), Type) %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(
      stacked_mean = cumsum(mean_count)) %>% ungroup() %>%
    dplyr::group_by(Label, Sample, Type) %>%
    dplyr::mutate(stacked_mean = ifelse(Type=="LF\nNo-plex", mean_count,stacked_mean))%>%
    dplyr::mutate(
      er_max = stacked_mean + se,
      er_min = stacked_mean - se
    ) %>%
    dplyr::ungroup()
  
  df_text <- dat_summary %>%
    dplyr::group_by(Type, Sample,Label) %>%
    dplyr::summarise(
      text_value = round(max(er_max) - max(se)),
      text_pos = max(er_max) + 310
    )
  
  
  ggplot(dat_summary, aes(x = Type, y = mean_count, fill = Sample)) +
    scale_x_discrete(limits = c("LF\nNo-plex","9-plexDIA", "9-plexDIA\n3-timePlex")) +
    
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "9-plexDIA\n3-timePlex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "9-plexDIA"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", color = "black", width = 0.35, alpha = 1) +
    
    geom_bar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
             aes(x = Type, y = mean_count, fill = Sample), alpha=0.8,
             stat = "identity", position = position_dodge(width = 1.17),
             color = "black", width = 1, alpha = 1) +
    
    scale_y_continuous(labels = comma) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type %in% c("9-plexDIA\n3-timePlex", "9-plexDIA")),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.18, size = 0.4) +
    
    geom_errorbar(data = dat_summary %>% dplyr::filter(Type == "LF\nNo-plex"),
                  aes(ymin = er_min, ymax = er_max),
                  width = 0.5, size = 0.4,
                  position = position_dodge(1.17)) +
    
    scale_fill_manual(values = my_colors) +
    theme_classic() +
    
    labs(y = "Protein data points / Run", x = "", fill = "Samples") +
    
    theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position = c(0.16, 0.9),
          legend.direction = "horizontal",
          legend.spacing.x = unit(1.0, 'line')) +
    
    guides(fill = guide_legend(override.aes = list(alpha = 0.7),
                               title = "Samples",
                               label.position = "bottom",
                               title.position = "top", title.vjust = 1, title.hjust = 0.5)) +
    
    geom_text(data = df_text %>% dplyr::filter(Type %in% c("9-plexDIA\n3-timePlex", "9-plexDIA")), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3) +
    
    geom_text(data = df_text %>% dplyr::filter(Type == "LF\nNo-plex"), 
              aes(x = Type, y = text_pos, label = comma(text_value)),
              col = 'black', size = 3,
              width = 0.5, size = 0.4,
              position = position_dodge(1.17))
  
  ggsave(paste0(extra,"_","Prot_coverage_","BestChannel_",use_best_channel,".pdf"),width=5.5,height=5.5)
  
  return(dat_counts)
  
}

pD_jmod_maxlfq <- function(df, group.header="prot", id.header = "seqcharge", quantity.header = "quant_val", min_num_reps=2,minimum_b_count=2){
  df <- df[df$Protein_Qvalue <= 0.01,]

  df$last_aa <- str_sub(df$stripped_seq, -1)
  
  df <- df %>%
    filter(
      (last_aa == "K" & Qvalue < 0.2 & BestChannel_Qvalue < 0.01) |
        (last_aa != "K" & Qvalue < 0.05 & BestChannel_Qvalue < 0.01)
    )
  
  df <- df %>%
    mutate(Type = case_when(
      timePlex == TRUE & Label == "LF" ~ "tp",
      timePlex == TRUE & Label != "LF" ~ "tp_pd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))
  
  df$bcount <- str_count(df$frag_names, "b")
  df <- df %>%
    filter(!(grepl("tag6", seq)) | (grepl("tag6", seq) & (!grepl("K", seq)) & bcount > (minimum_b_count - 1)))

  protein.groups <- diann_maxlfq(df, group.header=paste0(group.header), id.header = paste0(id.header), quantity.header = paste0(quantity.header))
  
  protein.groups <- as.data.frame(protein.groups)  
  protein.groups$prot <- row.names(protein.groups)
  return(protein.groups)
}


pD_melt_jmod_MaxLFQ <- function(dat,meta_bulk_path,infer_sample=T){
  pD_1 <- dat
  pD_3_m <- melt(as.matrix(pD_1))
  m <- fread(meta_bulk_path)
  m$Column <- as.numeric(m$Column)-1
  m$Run <- sub("_.*", "", m$Run)
  
  m <- m %>% dplyr::mutate("TP" = ifelse(timePlex==T,Column,
                                         ifelse(Label=="LF","NoPlex","")))
  m$run_chan <- paste0(m$Run, "_", m$TP,"_",m$Label)
  m$batch <- paste0(m$Label,"_",m$Column)
  m$run_timechannel <- paste0(m$Run,"_",m$Column)
  
  pD_3_m <- pD_3_m %>% left_join(m, by=c("Var2"="run_chan"))
  pD_3_m <- pD_3_m %>% dplyr::rename("prot"="Var1", "run_chan"="Var2")
  pD_3_m$norm_prot <- log2(as.numeric(pD_3_m$value))
  
  return(pD_3_m)
}

prot_normalizeQuant <- function(dat, quant_col="MS1_Int"){
  
  print(time)
  #normalization using only HUMAN proteins
  dat$value <- as.numeric(dat$value)
  datH <- dat %>%
    filter(!grepl("YEAST", prot)) %>%
    group_by(run_chan) %>%
    summarise(med_H = median(value, na.rm = TRUE), .groups = "drop")
  
  dat <- dat %>%
    left_join(datH, by = "run_chan") %>%
    mutate(med_norm = log2(.data[[quant_col]]) - log2(med_H)) %>%
    group_by(prot) %>%
    mutate(prot_mean = mean(.data[[quant_col]], na.rm = TRUE)) %>%
    mutate(norm = med_norm - mean(med_norm, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(norm) & norm != 0 & is.finite(norm))
  
  dat$quant_val <- (2^dat$norm)*dat$prot_mean
  return(dat)
}


Protein_level_quant <- function(normed_quant, grouper="prot", datatype = "protein",min_reps=2){

  normed_quant <- normed_quant %>%
    mutate(Type = case_when(
      timePlex == TRUE & Label != "LF" ~ "tp_pd",
      timePlex == FALSE & Label != "LF" ~ "pd",
      TRUE ~ "np"
    ))
  
  normed_quant$ProteinName <- normed_quant$prot
  normed_quant <- normed_quant %>% dplyr::mutate("species" = ifelse(grepl("HUMAN",prot),"Human",
                                                                    ifelse(grepl("YEAST",prot),"Yeast","remove"))) %>% filter(!grepl("remove",species))
  
  np <- normed_quant[normed_quant$Type=="np",]
  pD <- normed_quant[normed_quant$Type=="pd",] %>% filter(Rep == 2)
  tp_pd <- normed_quant[normed_quant$Type=="tp_pd",] %>% filter(Rep == 2)

  joined <- rbind(np,tp_pd, pD)
  
  normed_quant_avg_ints <- joined %>% dplyr::group_by(Type, .data[[grouper]], Sample) %>%
    dplyr::summarize("quant_avg" = mean(quant_val,na.rm=T)) %>% dplyr::rename("denom_sample"="Sample")
  
  human_medians <- joined %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(run_chan) %>%
    dplyr::summarise(median_Human = median(norm, na.rm = TRUE), .groups = "drop") %>% ungroup()
  
  joined <- joined %>%
    dplyr::left_join(human_medians, by = "run_chan") %>%  
    dplyr::mutate(norm = norm-median_Human) %>%  
    dplyr::select(-median_Human) 
  
  
  joined_quant <- joined %>% group_by(Type, Sample,.data[[grouper]],species) %>%
    dplyr::add_count() %>% dplyr::filter(n>(min_reps-1)) %>% #require atleast certain number of reps
    dplyr::summarize(med_quant = median(norm,na.rm=T))
  
  
  dynamic_column <- grouper  
  
  formula_str <- paste(dynamic_column, "species", sep = "+")
  joined_quant_d <- dcast(
    joined_quant, 
    formula = as.formula(paste(formula_str, "~ Type+Sample")), 
    value.var = "med_quant"
  )

  joined_quant_d$pd_JM <- joined_quant_d$pd_J-joined_quant_d$pd_M

  joined_quant_d$tp_pd_JM <- joined_quant_d$tp_pd_J-joined_quant_d$tp_pd_M

  joined_quant_d$np_JM <- joined_quant_d$np_J-joined_quant_d$np_M

  joined_quant_m <- joined_quant_d %>% dplyr::select(.data[[grouper]],species, 
                                                     pd_JM, tp_pd_JM, np_JM) %>% melt()
  
  human_medians <- joined_quant_m %>%
    dplyr::filter(species == "Human") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(median_Human = median(value, na.rm = TRUE), .groups = "drop") %>% dplyr::ungroup()
  
  joined_quant_m <- joined_quant_m %>%
    dplyr::left_join(human_medians, by = "variable") %>%  
    dplyr::mutate(value = value-median_Human) %>%  
    dplyr::select(-median_Human) %>% na.omit()
  
  joined_quant_m <- joined_quant_m %>% 
    dplyr::mutate(denom_sample = ifelse(grepl("J",variable),"J",
                                        ifelse(grepl("M",variable),"M", "K")))
  joined_quant_m$Type <- sub("_[^_]*$", "", joined_quant_m$variable)
  
  joined_quant_m <- joined_quant_m %>% dplyr::left_join(normed_quant_avg_ints, by =c("Type","denom_sample",paste0(grouper)))
  
  allcounts  <- joined_quant_m %>% dplyr::count(variable)
  
  
  ##### AB
  joined_AB <- joined_quant_m[grepl("JM", joined_quant_m$variable),]
  
  
  counts_AB <- allcounts[grepl("JM",allcounts$variable),]
  
  int_AB <- joined_AB %>% na.omit() %>% dplyr::count(.data[[grouper]]) %>% dplyr::filter(n==max(n))
  
  df_int <- data.frame(variable = "Intersected", n=nrow(int_AB))
  counts_AB <- rbind(counts_AB, df_int)
  
  
  counts_AB$variable <- factor(counts_AB$variable, 
                               levels=c("np_JM",
                                        "pd_JM",
                                        "tp_pd_JM",
                                        "Intersected"))
  
  counts_AB_plot <- ggplot(counts_AB, aes(x=as.factor(variable), y=n, fill=variable)) + 
    geom_bar(stat="identity", colour="black", alpha =0.7) +
    geom_text(data=counts_AB, aes(x=as.factor(variable), y=n+30, label=n), col='black', size=4.5) +
    theme_classic() + 
    theme(
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14),
      axis.ticks.x = element_blank(),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.title = element_text(size = 13.5, hjust = 0.5, face="bold"),
      legend.position = "none",
      strip.text.x = element_text(size=12, face = "italic"),
      legend.text = element_text(size=12))+
    labs(x = "", y = "Protein Ratios (A/B)", fill = "") +
    scale_fill_manual(values = c("#413c58", "#b6db8a", "#d87c7c","grey40"))
  counts_AB_plot
  
  
  joined_AB_int <- joined_AB %>% dplyr::filter(.data[[grouper]]%in%int_AB[[grouper]]) #%>% na.omit()
  joined_AB_int$variable <- factor(joined_AB_int$variable, 
                                   levels=c("np_JM",
                                            "pd_JM",
                                            "tp_pd_JM",
                                            "Intersected"))
  
  medians <- joined_AB_int %>%
    dplyr::group_by(species, variable) %>%
    dplyr::summarize(median_value = median(value), .groups = "drop")
  
  denseplot <- ggplot(joined_AB_int, aes(fill=species, y=value)) +
    facet_grid(~variable) +
    geom_density(alpha=0.75) + 
    geom_hline(data=medians, aes(yintercept=median_value), linetype="solid", color="black") + 
    geom_hline(yintercept = -2, linetype="dashed", color="#e56730") +
    geom_hline(yintercept = 0,linetype="dashed",color="#109abf") + 
    labs(subtitle="")+
    theme_bw() +
    scale_fill_manual(values=c("#68aabc","#e08a6c"))+
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  
  quantplot_AB <- ggarrange(counts_AB_plot, denseplot,
                            ncol = 2, nrow = 1,
                            widths = c(0.4, 0.7))
  quantplot_AB
  ggsave(paste0("prot_level_3x9",datatype,".pdf"),width=6,height=4)
  
  df1 <- data.frame(type = c("1np","9pd","9pd_3tp"),
                   plex = c(1, 9, 27))
  
  ggplot(df1, aes(x=type,y=plex,fill=type)) + geom_point(size=4,color="black",shape=21) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position="none"
    ) +
  scale_fill_manual(values = c("#413c58", "#b6db8a", "#d87c7c"))
  ggsave(paste0("plex",datatype,".pdf"),width=4.1,height=1.7)
  
}





