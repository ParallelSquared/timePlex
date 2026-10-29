################################################################################
# Manuscript figure style: Helvetica text, print-size panels, shared colors
################################################################################


## ---- Colors ------------------------------------------------------------------

paper_cols <- list(
  species  = c(Human = "#68AABC", Yeast = "#E08A6C", Arabidopsis = "#7FDD6C"),
  expected = c(Human = "#109ABF", Yeast = "#E56730"),   # dashed expected log2 ratios
  samples3 = c(A = "#7CA982", B = "#4E598C", C = "#FED766"),
  samples4 = c(A = "#EC008C", B = "#2BB673", C = "#00AEEF", D = "#FBB040")
)


## ---- Text sizes (pt) -----------------------------------------------------------
# Sizes each text has when its figure is 180 mm wide, measured from the
# manuscript (Derks_timePlex.pdf, Oct 2026) relative to the figure width.
# They differ by figure: Figs. 3-4 are the smallest, Extended Data figures are
# 1.3-2x larger.
paper_fig_sizes <- list(
  Fig1 = c(text = 10, title = 11, tag = 19),
  Fig2 = c(text = 6,  title = 8,  tag = 14),
  Fig3 = c(text = 6,  title = 7,  tag = 15),
  Fig4 = c(text = 7,  title = 8,  tag = 14),   
  Fig5 = c(text = 7,  title = 8,  tag = 14),
  Fig6 = c(text = 7,  title = 9,  tag = 18),
  ED1  = c(text = 8,  title = 10, tag = 18),   
  ED2  = c(text = 8,  title = 9,  tag = 19),
  ED3  = c(text = 7,  title = 8,  tag = 17),
  ED4  = c(text = 9,  title = 11, tag = 20),
  ED5  = c(text = 8,  title = 9,  tag = 20),
  ED6  = c(text = 6,  title = 7,  tag = 13),
  ED7  = c(text = 8,  title = 9,  tag = 17),   
  ED8  = c(text = 11, title = 14, tag = 21),   
  ED10 = c(text = 8,  title = 9,  tag = 14)    
)

# Set the figure whose sizes paper_size() returns, until the calling function
# returns: a name ("Fig5", "ED3") or an output path ("Figs/ExtDataFig3/x.pdf").
# Figures not in the table (e.g. new Extended Data figures) get the Fig. 4 sizes.
paper_figure <- function(fig, env = parent.frame()){
  parts <- sub("^ExtDataFig", "ED", strsplit(fig[[1]], "/", fixed = TRUE)[[1]])
  key <- utils::tail(c("Fig4", parts[parts %in% names(paper_fig_sizes)]), 1)
  withr::local_options(list(paper.figure = key), .local_envir = env)
  invisible(key)
}
paper_size <- function(what = "text") paper_fig_sizes[[getOption("paper.figure", "Fig4")]][[what]]

paper_pt <- paper_fig_sizes$Fig4             # Fig. 4 sizes; use paper_size() in plotting code
paper_geom_text <- function(pt = paper_size("text")) pt / ggplot2::.pt   # geom_text() sizes are in mm


## ---- Themes ------------------------------------------------------------------

paper_theme_common <- function(base_size, frame = 0.25) {
  text <- paper_size("text")
  theme(text              = element_text(colour = "black", size = base_size),
        axis.text         = element_text(colour = "black", size = text),
        axis.title        = element_text(size = base_size),
        axis.ticks        = element_line(colour = "black", linewidth = frame),
        axis.ticks.length = unit(2, "pt"),
        legend.title      = element_text(size = base_size),
        legend.text       = element_text(size = text),
        legend.key.size   = unit(1.1 * text, "pt"),
        legend.key        = element_blank(),
        legend.background = element_blank(),
        legend.margin     = margin(0, 0, 0, 0),
        strip.text        = element_text(size = base_size, colour = "black",
                                         margin = margin(1.5, 1.5, 1.5, 1.5)),
        plot.title        = element_text(size = base_size, hjust = 0.5,
                                         margin = margin(0, 0, 2, 0)),
        plot.subtitle     = element_text(size = text, hjust = 0,
                                         margin = margin(0, 0, 2, 0)),
        plot.tag          = element_text(size = paper_size("tag"), face = "bold"),
        plot.background   = element_blank(),
        panel.background  = element_blank(),
        plot.margin       = margin(2, 2, 2, 2))
}

# Axis lines, no panel border (bar charts, scatter plots, TICs)
paper_theme <- function(base_size = paper_size("title"), frame = 0.25) {
  theme_classic(base_size = base_size, base_line_size = frame, base_rect_size = frame) +
    paper_theme_common(base_size, frame) +
    theme(axis.line = element_line(colour = "black", linewidth = frame))
}

# Panel border and grey strips, no grid (faceted panels, ridgelines, XICs)
paper_theme_box <- function(base_size = paper_size("title"), grid = FALSE, frame = 0.25) {
  t <- theme_bw(base_size = base_size, base_line_size = frame, base_rect_size = frame) +
    paper_theme_common(base_size, frame) +
    theme(panel.border     = element_rect(colour = "grey20", fill = NA, linewidth = frame),
          strip.background = element_rect(fill = "grey85", colour = "grey20", linewidth = frame))
  if (!grid) t <- t + theme(panel.grid = element_blank())
  else t <- t + theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
                      panel.grid.minor = element_blank())
  t
}


## ---- Axis helpers ------------------------------------------------------------

# Whole-number breaks every `by` units, labelled without decimals
paper_int_breaks <- function(by = 2) {
  function(lim) seq(ceiling(lim[1] / by) * by, floor(lim[2] / by) * by, by = by)
}
paper_int_labels <- function(x) scales::comma(x, accuracy = 1)


## ---- Saving ------------------------------------------------------------------
paper_pdf_device <- function(filename, width, height, bg = "transparent", ...) {
  grDevices::pdf(filename, width = width, height = height, bg = bg,
                 family = paper_helvetica(), encoding = "CP1253", useDingbats = FALSE)
}

paper_helvetica <- function(){
  dir <- file.path(tempdir(), "paper_helvetica")
  faces <- c("Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique")
  files <- file.path(dir, paste0(faces, ".afm"))
  if (!all(file.exists(files))) {
    dir.create(dir, showWarnings = FALSE)
    after <- c("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "R")
    for (k in seq_along(faces)) {
      x <- readLines(gzfile(file.path(R.home("library"), "grDevices", "afm", paste0(faces[k], ".afm.gz"))))
      metric <- function(g) x[grep(paste0(" N ", g, " ;"), x)[1]]
      left <- as.numeric(sub(".* B (-?[0-9]+) .*", "\\1", sapply(after, metric)))   # left side bearings
      kern <- round((if (grepl("Bold", faces[k])) 482 else 475) - left)
      x <- append(x, c(paste("KPX Delta", after, kern), "KPX Delta Delta 0"), after = grep("^EndKernPairs", x) - 1)
      n <- grep("^StartKernPairs", x)
      x[n] <- paste("StartKernPairs", length(grep("^KPX", x)))
      writeLines(x, files[k])
    }
  }
  grDevices::Type1Font("Helvetica", files)
}

paper_pdf <- function(file, width, height, draw) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  old <- Sys.getlocale("LC_CTYPE")
  suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
  on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE", old)), add = TRUE)
  paper_pdf_device(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  draw()
  invisible(file)
}

# Save a plot at print size (inches)
paper_save <- function(plot, file, width, height) {
  g <- paper_grob(plot)
  paper_pdf(file, width, height, function() grid::grid.draw(g))
}

# Drop x tick labels, ticks and title (rows that share the x axis below them)
paper_drop_x <- function(){
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
}


## ---- Panels for the assembled figures ----------------------------------------

paper_grob <- function(plot){
  old <- Sys.getlocale("LC_CTYPE")
  suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
  on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE", old)), add = TRUE)
  grDevices::pdf(NULL, family = paper_helvetica(), encoding = "CP1253")
  on.exit(grDevices::dev.off(), add = TRUE)
  g <- if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot)
       else if (inherits(plot, "ggplot")) ggplot2::ggplotGrob(plot)
       else plot
  paper_dashes(g)
}

paper_dashes <- function(g){
  if (inherits(g, "text") && is.character(g$label)) {
    g$label <- gsub("\u2212", "-", g$label, fixed = TRUE)
    g$label <- gsub("(?<=[[:alnum:]])-(?=[[:alpha:]])", "\u00ad", g$label, perl = TRUE)
  }
  if (inherits(g, "gtable")) {
    g$grobs <- lapply(g$grobs, paper_dashes)
  } else if (inherits(g, "gTree") && length(g$children)) {
    g$children <- do.call(grid::gList, lapply(g$children, paper_dashes))
  }
  g
}

# Save a manuscript panel at print size and keep its grob and size, so the
# assembled figure (figure_assembly.R) draws exactly the same panel
paper_panel <- function(plot, file, width, height){
  g <- paper_grob(plot)
  paper_pdf(file, width, height, function() grid::grid.draw(g))
  list(grob = g, width = width, height = height, file = file)
}


## ---- Dense layers as images --------------------------------------------------


paper_raster <- function(layer){
  if (!isTRUE(getOption("paper.raster", TRUE))) return(layer)
  geom <- layer$geom
  ggplot2::ggproto(NULL, layer, geom = ggplot2::ggproto(NULL, geom,
    draw_panel = function(self, ...) {
      grid::gTree(children = grid::gList(ggplot2::ggproto_parent(geom, self)$draw_panel(...)),
                  cl = "paper_raster")
    }))
}


makeContext.paper_raster <- function(x){
  w <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE)
  h <- grid::convertHeight(grid::unit(1, "npc"), "in", valueOnly = TRUE)
  dev <- grDevices::dev.cur()
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = w, height = h, units = "in", res = getOption("paper.raster.dpi", 600),
                background = "transparent", bitsize = 16)
  grid::pushViewport(grid::viewport())     # points are in native units: 0-1 across the panel
  grid::grid.draw(grid::gTree(children = x$children))
  grid::popViewport()
  grDevices::dev.off()
  grDevices::dev.set(dev)
  img <- png::readPNG(tmp)
  unlink(tmp)
  grid::rasterGrob(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = FALSE)
}

paper_box_outliers <- function(alpha, size = 0.5, colour = "grey20"){
  paper_raster(geom_boxplot(alpha = 0, colour = NA, outlier.colour = colour, outlier.alpha = alpha,
                            outlier.size = size, outlier.stroke = 0, show.legend = FALSE))
}
