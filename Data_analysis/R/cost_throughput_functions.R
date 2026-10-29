################################################################################
# Fig. 1: samples per run and cost per sample with plexDIA x timePlex
################################################################################

throughput_plot <- function(){
  
  pD_tp <- expand.grid(
    pd = 1:10,
    tp = 1:10
  )
  
  pD_tp$samples_per_run <- pD_tp$pd * pD_tp$tp
  
  ggplot(pD_tp, aes(x = pd, y = tp, fill = samples_per_run)) +
    geom_tile() +
    scale_fill_viridis_c(option = "mako") +
    theme_classic() +
    labs(
      x = "plexDIA",
      y = "timePlex",
      fill = "Samples / run"
    ) +
    coord_fixed()
  
  ggsave("Throughput_samplesPerRun.pdf", width = 3.5, height = 3.5)
}

pD_cost <- function(d0=778.64, d4=2348.04, d8=2343.62, numcolumns=3,pack10columns=10160,runspercolumn=2500,
                    out_file = "Figs/Fig1/Fig1b.pdf"){
  #Figure 1b:
  paper_figure(out_file)
  
  total_prot <- 50*3*0.1 #50 aliquots * 3 labels * 0.1 mg protein/aliquot
  
  total_cost <- d0+d4+d8
  
  cost_per_ug <- total_cost/(total_prot*1000)
  cost_extra_columns <- ((pack10columns/10)*(numcolumns-1))/runspercolumn
  # costs per run divided by samples per run: numcolumns for timePlex, 3 * numcolumns combined with 3-plexDIA
  plexDIA_timePlex <- (100*1.5)/(3*numcolumns)
  timePlex <- (100*1.5)/numcolumns
  plexDIA <- (100*1.5)/3
  LFDIA <- (100*1.5)
  cost <- c(LFDIA, plexDIA, (cost_per_ug), timePlex, cost_extra_columns/numcolumns, plexDIA_timePlex, (cost_per_ug), cost_extra_columns/(3*numcolumns))
  Application_Reagent <- c("LC-MS/MS cost","LC-MS/MS cost","Mass tag cost" ,"LC-MS/MS cost", "Extra columns cost","LC-MS/MS cost","Mass tag cost","Extra columns cost")
  type <- c("LF-DIA", "3-plexDIA","3-plexDIA","3-timePlex","3-timePlex", "3-timePlex & 3-plexDIA","3-timePlex & 3-plexDIA","3-timePlex & 3-plexDIA")
  
  df <- data.frame(cost = cost, Application_Reagent = Application_Reagent,type = type)
  
  # first level is stacked at the end of each bar
  df$Application_Reagent <- factor(df$Application_Reagent, levels=c("Extra columns cost", "Mass tag cost", "LC-MS/MS cost"))
  df$type <- factor(df$type, levels=c("LF-DIA", "3-plexDIA", "3-timePlex", "3-timePlex & 3-plexDIA"))
  
  p <- ggplot(df, aes(x=type, y=cost, fill=Application_Reagent)) + 
    geom_col(width = 0.45) +
    scale_y_continuous(labels=scales::dollar_format(), expand = expansion(mult = c(0.02, 0.03))) + 
    scale_fill_manual(values = c("LC-MS/MS cost" = "#A7A9AC", "Mass tag cost" = "red",
                                 "Extra columns cost" = "blue")) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = NULL, y = "Estimated cost (USD) per sample") +
    coord_flip() +
    paper_theme() +
    theme(legend.position = "inside",          # as in the manuscript: right of the short bars
          legend.position.inside = c(0.58, 1),
          legend.justification = c(0, 1),
          legend.title = element_blank(),
          legend.key.size = unit(14, "pt"),
          legend.key.spacing.y = unit(0, "pt"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 13, face = "bold"),
          axis.title.x = element_text(size = paper_size("title")))
  
  invisible(list(fig1b = paper_panel(p, out_file, width = 6.1, height = 1.6)))
}

pD_cost_plexDIA <- function(d0=778.64, allothers=2348.04,midplex=10,highplex=100){
  #Figure 1d:
  
  ## 3-plexDIA
  total_prot <- 50*3*0.1 #50 aliquots * 3 labels * 0.1 mg protein/aliquot
  total_cost <- d0+(3-1)*allothers
  
  plex_cost_per_ug_total <- (total_cost/(total_prot*1000))
  
  # 10-plexDIA
  midplex_total_prot <- 50*midplex*0.1
  midplex_total_cost <- d0+(midplex-1)*allothers
  
  midplex_cost_per_ug_total <- (midplex_total_cost/(midplex_total_prot*1000))
  
  # 100-plexDIA
  highplex_total_prot <- 50*highplex*0.1
  highplex_total_cost <- d0+(highplex-1)*allothers
  
  highplex_cost_per_ug_total <- (highplex_total_cost/(highplex_total_prot*1000))
  
  LFDIA <- (100*1.5)
  cost <- c(LFDIA, LFDIA/3, plex_cost_per_ug_total, LFDIA/midplex, midplex_cost_per_ug_total, LFDIA/highplex,highplex_cost_per_ug_total)
  Application_Reagent <- c("LC-MS/MS cost","LC-MS/MS cost","Label cost" ,"LC-MS/MS cost","Label cost","LC-MS/MS cost","Label cost")
  type <- c("LF-DIA", "3-plexDIA","3-plexDIA",paste0(midplex,"-plexDIA"),paste0(midplex,"-plexDIA"), 
            paste0(highplex,"-plexDIA"),paste0(highplex,"-plexDIA"))
  
  df <- data.frame(cost = cost, Application_Reagent = Application_Reagent,type = type)
  
  df$Application_Reagent <- factor(df$Application_Reagent, levels=c("Extra columns cost", "Label cost", "LC-MS/MS cost"))
  df$type <- factor(df$type, levels=c("LF-DIA", "3-plexDIA", paste0(midplex,"-plexDIA"), paste0(highplex,"-plexDIA")))
  
  total_costs_vec <- c(LFDIA, LFDIA/3 + plex_cost_per_ug_total, LFDIA/midplex + midplex_cost_per_ug_total, LFDIA/highplex + highplex_cost_per_ug_total)
  type_unique <- c("LF-DIA", "3-plexDIA", paste0(midplex,"-plexDIA"), paste0(highplex,"-plexDIA"))
  df_unique <- data.frame(totalcosts = total_costs_vec, type=type_unique)
  
  ggplot(df, aes(x=type, y=cost, fill=type)) + 
    theme_classic() + labs(x="", y="Estimated cost (USD) per sample",subtitle = "Assuming 1 ug protein per sample") +
    geom_bar(data = df %>% dplyr::filter(type == "LF-DIA"),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == "3-plexDIA"),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(midplex,"-plexDIA")),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(highplex,"-plexDIA")),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    scale_y_continuous(labels=scales::dollar_format()) + 
    scale_fill_manual(values = c("grey39","red","green")) +
    theme(legend.position = c(0.7,0.8),
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          axis.text.x = element_text(size=13, face = "bold", color = "black"),
          axis.text.y = element_text(size=14, face = "bold", color = "black"),
          axis.title.x = element_text(size=12, color = "black")) + coord_flip()

  
  cost <- c(LFDIA, LFDIA/3, plex_cost_per_ug_total/10, LFDIA/midplex, midplex_cost_per_ug_total/10, LFDIA/highplex,highplex_cost_per_ug_total/10)
  Application_Reagent <- c("LC-MS/MS cost","LC-MS/MS cost","Label cost" ,"LC-MS/MS cost","Label cost","LC-MS/MS cost","Label cost")
  type <- c("LF-DIA", "3-plexDIA","3-plexDIA",paste0(midplex,"-plexDIA"),paste0(midplex,"-plexDIA"), 
            paste0(highplex,"-plexDIA"),paste0(highplex,"-plexDIA"))
  
  df <- data.frame(cost = cost, Application_Reagent = Application_Reagent,type = type)
  
  df$Application_Reagent <- factor(df$Application_Reagent, levels=c("Extra columns cost", "Label cost", "LC-MS/MS cost"))
  df$type <- factor(df$type, levels=c("LF-DIA", "3-plexDIA", paste0(midplex,"-plexDIA"), paste0(highplex,"-plexDIA")))
  
  
  ggplot(df, aes(x=type, y=cost, fill=type)) + 
    theme_classic() + labs(x="", y="Estimated cost (USD) per sample",subtitle = "Assuming 100 ng protein per sample") +
    geom_bar(data = df %>% dplyr::filter(type == "LF-DIA"),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == "3-plexDIA"),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(midplex,"-plexDIA")),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    geom_bar(data = df %>% dplyr::filter(type == paste0(highplex,"-plexDIA")),
             stat="identity", aes(x=type, y =cost, fill=Application_Reagent), width = 0.45, alpha =1) + 
    scale_y_continuous(labels=scales::dollar_format()) + 
    scale_fill_manual(values = c("grey39","red","blue")) +
    theme(legend.position = c(0.7,0.8),
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          axis.text.x = element_text(size=13, face = "bold", color = "black"),
          axis.text.y = element_text(size=14, face = "bold", color = "black"),
          axis.title.x = element_text(size=12, color = "black")) + coord_flip()
  
  ggsave("Figs/Figure1_cost_plexDIAOnly.pdf", width=4.1, height=1.9, dpi=500)
}

