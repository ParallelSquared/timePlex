# Fig. 1: multiplexing in the time and mass domains is combinatorial

list(

  ## ---- samples per run for plexDIA x timePlex ---------------------------------
  tar_target(throughputplot, throughput_plot(), format="qs"),

  ## ---- Fig. 1b: cost per sample -----------------------------------------------
  tar_target(cost, pD_cost(d0=778.64, d4=2348.04, d8=2343.62, numcolumns=3, pack10columns=10160,runspercolumn=1000), format="qs"),
  tar_target(cost_pd_only, pD_cost_plexDIA(d0=778.64, allothers=2348.04,midplex=10,highplex=100), format="qs"),

  ## ---- assembled figure (Figs/Fig1/Fig1_full.pdf; a is a schematic, left blank) ----
  tar_target(fig1_full, paper_fig1(cost$fig1b), format = "file")
)

