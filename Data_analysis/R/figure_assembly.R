################################################################################
# Assembled manuscript figures: code-made panels placed at print size, with
# blank space (and the panel letter) left for schematics, tables and photos
################################################################################
# Figures are laid out at the journal's double-column width (180 mm = 7.09 in).
# Positions are in inches from the top-left corner of the figure. Each panel
# comes from paper_panel(), so it is the same object as its own PDF.


## ---- Generic placement --------------------------------------------------------

# One placed panel: `panel` from paper_panel() (NULL leaves the space empty),
# its top-left corner (x, y) and the panel letter's position
paper_place <- function(panel, x, y, letter = NULL, lx = x, ly = y){
  list(panel = panel, x = x, y = y, letter = letter, lx = lx, ly = ly)
}

# Items can be named (list(a = paper_place(...))) so that decorate(), run after
# the panels are drawn, can find their plot panels with paper_panel_boxes(name)
paper_assemble <- function(file, width, height, items, decorate = NULL){
  paper_figure(file)                                  # panel letters at the figure's size
  # panels come from other targets; ggtext labels (Figs. 4c-e v2, 5g) are sized
  # and drawn by gridtext's grob methods, which a fresh session has not loaded
  requireNamespace("gridtext", quietly = TRUE)
  ids <- if (is.null(names(items))) seq_along(items) else ifelse(nzchar(names(items)), names(items), seq_along(items))
  paper_pdf(file, width, height, function(){
    for (i in seq_along(items)) {
      it <- items[[i]]
      if (is.null(it$panel)) next
      grid::pushViewport(grid::viewport(x = grid::unit(it$x, "in"), y = grid::unit(height - it$y, "in"),
                                        width = grid::unit(it$panel$width, "in"),
                                        height = grid::unit(it$panel$height, "in"),
                                        just = c("left", "top"), name = paste0("paper-item-", ids[i])))
      grid::grid.draw(it$panel$grob)
      grid::upViewport()      # kept in the viewport tree for decorate()
    }
    if (is.function(decorate)) decorate()
    for (it in items) {
      if (is.null(it$letter)) next
      grid::grid.text(it$letter, x = grid::unit(it$lx, "in"), y = grid::unit(height - it$ly, "in"),
                      just = c("left", "top"),
                      gp = grid::gpar(fontsize = paper_size("tag"), fontface = "bold"))
    }
  })
}


# Page positions (inches from the bottom-left) of a placed item and of the
# plot panels drawn inside it; for use inside a decorate() function
paper_item_box <- function(item){
  grid::upViewport(0)
  grid::downViewport(paste0("paper-item-", item))
  ll <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
  ur <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
  grid::upViewport(0)
  data.frame(x0 = ll$x, x1 = ur$x, y0 = ll$y, y1 = ur$y)
}

paper_panel_boxes <- function(item){
  vp <- paste0("paper-item-", item)
  grid::upViewport(0)
  grid::downViewport(vp)
  tree <- paste(utils::capture.output(print(grid::current.vpTree(all = FALSE))), collapse = "")
  nm <- regmatches(tree, gregexpr("viewport\\[panel[^]]*\\]", tree))[[1]]
  nm <- unique(sub("^viewport\\[", "", sub("\\]$", "", nm)))
  boxes <- do.call(rbind, lapply(nm, function(v){
    grid::upViewport(0); grid::downViewport(vp); grid::downViewport(v)
    ll <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    ur <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    data.frame(name = v, x0 = ll$x, x1 = ur$x, y0 = ll$y, y1 = ur$y)
  }))
  grid::upViewport(0)
  boxes[order(boxes$x0, -boxes$y1), ]
}

# Grey zoom wedge from a span (x_left..x_right at y_top) down to a point, in page inches
paper_wedge <- function(x_left, x_right, y_top, x_apex, y_apex){
  grid::grid.polygon(x = grid::unit(c(x_left, x_right, x_apex), "in"),
                     y = grid::unit(c(y_top, y_top, y_apex), "in"),
                     gp = grid::gpar(fill = grid::linearGradient(c("white", "grey72"),
                                                                 x1 = 0.5, y1 = 1, x2 = 0.5, y2 = 0),
                                     col = "grey60", lwd = 0.4))
}

# Zoom wedges from each XIC panel column (items left to right, columns left to
# right) to its band in the TIC (tic$bands, in the same order)
paper_zoom_wedges <- function(xic_items, tic_item, tic, gap = 0.03){
  tb <- paper_panel_boxes(tic_item)[1, ]
  band_x <- tb$x0 + (tic$bands - tic$x_range[1]) / diff(tic$x_range) * (tb$x1 - tb$x0)
  cols <- do.call(rbind, lapply(xic_items, function(it){
    b <- paper_panel_boxes(it)
    b$col <- cumsum(c(TRUE, diff(round(b$x0, 2)) > 0))        # panels sharing an x position
    top <- paper_item_box(it)$y0 - gap
    do.call(rbind, lapply(split(b, b$col), function(cb) data.frame(x0 = min(cb$x0), x1 = max(cb$x1), y = top)))
  }))
  n <- min(nrow(cols), length(band_x))
  for (i in seq_len(n)) paper_wedge(cols$x0[i], cols$x1[i], cols$y[i], band_x[i], tb$y1)
}


## ---- Fig. 1: a is a schematic -------------------------------------------------

paper_fig1 <- function(b, file = "Figs/Fig1/Fig1_full.pdf"){
  paper_assemble(file, width = 7.09, height = 6.2 + b$height, items = list(
    paper_place(NULL, 0, 0, "a", 0, 0),
    paper_place(b, 0.25, 6.2, "b", 0, 6.2)))        # b ends at ~90% of the width, as in the manuscript
}


## ---- Fig. 2: a-c are schematics, d-e come from the RT model (Python) -------------

paper_fig2 <- function(f, file = "Figs/Fig2/Fig2_full.pdf"){
  paper_assemble(file, width = 7.09, height = 4.85 + f$height, items = list(
    paper_place(NULL, 0, 0, "a", 0, 0),
    paper_place(NULL, 0, 0, "b", 4.85, 0),
    paper_place(NULL, 0, 0, "c", 0, 2.6),
    paper_place(NULL, 0, 0, "d", 5.05, 2.6),
    paper_place(NULL, 0, 0, "e", 0, 4.85),
    paper_place(f, 7.09 - f$width, 4.85, "f", 7.09 - f$width - 0.1, 4.85)))
}


## ---- Fig. 3: c is a table; the zoom wedges in b are drawn by hand -----------------

paper_fig3 <- function(a, b_human, b_yeast, b_tic, d, e, file = "Figs/Fig3/Fig3_full.pdf"){
  # b as in the manuscript: the XIC groups centred, the TIC below them and grey
  # zoom wedges from each XIC panel to its band in the TIC
  W <- 7.09
  y_b   <- a$height + 0.05                       # XIC titles start right under a
  y_tic <- y_b + b_human$height + 0.3            # room for the wedges
  y_d   <- y_tic + b_tic$height + 0.08
  x_h   <- (W - b_human$width - b_yeast$width - 0.1) / 2
  paper_assemble(file, width = W, height = y_d + d$height, items = list(
    a       = paper_place(a,       0,                       0,     "a", 0, 0),
    b_human = paper_place(b_human, x_h,                     y_b,   "b", 0, y_b),
    b_yeast = paper_place(b_yeast, x_h + b_human$width + 0.1, y_b),
    b_tic   = paper_place(b_tic,   W - b_tic$width - 0.08,  y_tic),
    c       = paper_place(NULL,    0, 0,                           "c", 0, y_d),
    d       = paper_place(d,       W - e$width - d$width,   y_d,   "d", W - e$width - d$width, y_d),
    e       = paper_place(e,       W - e$width,             y_d,   "e", W - e$width, y_d)),
    decorate = function() paper_zoom_wedges(c("b_human", "b_yeast"), "b_tic", b_tic))
}


## ---- Fig. 4: quantitative accuracy (all panels from code) --------------------

# Left edge (inches) of the first plot in a panel, right of its letter or label
# column: the panel is drawn on a null device and its plot panels located as in
# paper_panel_boxes() (gtable widths alone miss markdown titles sized at draw time)
paper_plot_left <- function(panel){
  requireNamespace("gridtext", quietly = TRUE)
  grDevices::pdf(NULL, width = panel$width, height = panel$height, family = paper_helvetica(), encoding = "CP1253")
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(name = "paper-measure"))
  grid::grid.draw(panel$grob)
  tree <- paste(utils::capture.output(print(grid::current.vpTree(all = FALSE))), collapse = "")
  nm <- unique(sub("^viewport\\[", "", sub("\\]$", "", regmatches(tree, gregexpr("viewport\\[panel-[0-9][^]]*\\]", tree))[[1]])))
  x0 <- sapply(nm, function(v){
    grid::upViewport(0); grid::downViewport("paper-measure"); grid::downViewport(v)
    grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)$x
  })
  sort(unique(round(x0, 4)))[2]
}

# Resize a panel's first (letter or label) column so its first plot starts `left` inches in
paper_align_left <- function(panel, left){
  col <- min(panel$grob$layout$l[grepl("^panel-[0-9]", panel$grob$layout$name)])
  panel$grob$widths[col] <- panel$grob$widths[col] + grid::unit(left - paper_plot_left(panel), "in")
  panel
}

paper_fig4 <- function(a, b, ce, f, file = "Figs/Fig4/Fig4_full.pdf"){
  # a and c-e share the left edge, b and f the right edge; letters sit in the
  # panels' empty top-left corners (c-e carries its own letters); a's plots start where c-e's do
  W <- 7.09
  a <- paper_align_left(a, paper_plot_left(ce))
  paper_assemble(file, width = W, height = a$height + ce$height, items = list(
    paper_place(a,  0,             0,        "a", 0, 0),
    paper_place(b,  W - b$width,   0,        "b", W - b$width, 0),
    paper_place(ce, 0,             a$height),
    paper_place(f,  W - f$width,   a$height, "f", W - f$width, a$height)))
}


## ---- Fig. 5: a is a table; the zoom wedges in b are drawn by hand -----------------

paper_fig5 <- function(b_yeast, b_human, b_tic, c_prec, c_prot, d, e, f, g, h,
                       file = "Figs/Fig5/Fig5_full.pdf", gap = 0.09){
  # small gaps between neighbouring panels; zoom wedges from the XIC columns
  # down to their bands in the TIC, as in the manuscript
  W <- 7.09
  x_b   <- W - b_yeast$width - b_human$width     # b fills the right ~55%
  y_c   <- 1.55
  y_d   <- y_c + c_prec$height + gap
  y_tic <- b_yeast$height + 0.32                 # room for the wedges
  y_e   <- max(y_tic + b_tic$height, y_d + d$height) - 0.05   # just below the axis titles of d and the TIC
  sp    <- (W - e$width - f$width - g$width - h$width) / 3
  x_f   <- e$width + sp; x_g <- x_f + f$width + sp; x_h <- x_g + g$width + sp
  paper_assemble(file, width = W, height = y_e + max(e$height, f$height, g$height, h$height), items = list(
    a       = paper_place(NULL,    0, 0,                  "a", 0, 0),
    b_yeast = paper_place(b_yeast, x_b, 0,                "b", x_b, 0),
    b_human = paper_place(b_human, x_b + b_yeast$width, 0),
    b_tic   = paper_place(b_tic,   W - b_tic$width, y_tic),
    c_prec  = paper_place(c_prec,  0, y_c,                "c", 0, y_c),
    c_prot  = paper_place(c_prot,  c_prec$width + gap, y_c),
    d       = paper_place(d,       0, y_d,                "d", 0, y_d),
    e       = paper_place(e,       0,   y_e,              "e", 0,   y_e),
    f       = paper_place(f,       x_f, y_e,              "f", x_f, y_e),
    g       = paper_place(g,       x_g, y_e,              "g", x_g, y_e),
    h       = paper_place(h,       x_h, y_e,              "h", x_h, y_e)),
    decorate = function() paper_zoom_wedges(c("b_yeast", "b_human"), "b_tic", b_tic))
}


## ---- Fig. 6: a is a schematic, b a photo -------------------------------------------

paper_fig6 <- function(cf, file = "Figs/Fig6/Fig6_full.pdf"){
  # letters c-f sit at the tops of the stacked panels (legends + c, d, e, f)
  y0 <- 1.6
  ys <- y0 + cumsum(c(0, 0.7 + 1.75, 1.3, 1.25))     # legends (0.7 in) + c, d, e, f
  paper_assemble(file, width = 7.09, height = y0 + cf$height, items = list(
    paper_place(NULL, 0, 0, "a", 0, 0),
    paper_place(NULL, 0, 0, "b", 4.9, 0),
    paper_place(cf,   0, y0, "c", 0, ys[1] + 0.2),     # level with the legend rows, as in the manuscript
    paper_place(NULL, 0, 0,  "d", 0, ys[2]),
    paper_place(NULL, 0, 0,  "e", 0, ys[3]),
    paper_place(NULL, 0, 0,  "f", 0, ys[4])))
}


## ---- Extended Data figures (all panels from code) --------------------------------

paper_ed1 <- function(rows, file = "Figs/ExtDataFig1/ExtDataFig1_full.pdf"){
  y <- c(0, cumsum(sapply(rows, `[[`, "height")))
  items <- lapply(seq_along(rows), function(i) paper_place(rows[[i]], 0, y[i], letters[i], 0, y[i]))
  paper_assemble(file, width = 7.09, height = y[length(y)], items = items)
}

paper_ed2 <- function(a, b, c, file = "Figs/ExtDataFig2/ExtDataFig2_full.pdf"){
  h <- a$height + b$height
  paper_assemble(file, width = 7.09, height = h, items = list(
    paper_place(a, 0, 0,        "a", 0, 0),
    paper_place(b, 0, a$height, "b", 0, a$height),
    paper_place(c, 7.09 - c$width, (h - c$height) / 2, "c", 7.09 - c$width, (h - c$height) / 2)))
}

paper_ed3 <- function(a, b, file = "Figs/ExtDataFig3/ExtDataFig3_full.pdf"){
  paper_assemble(file, width = 7.09, height = max(a$height, b$height), items = list(
    paper_place(a, 0, 0,       "a", 0, 0),
    paper_place(b, a$width, 0, "b", a$width, 0)))
}

paper_ed4 <- function(a, b, file = "Figs/ExtDataFig4/ExtDataFig4_full.pdf"){
  # b is aligned with the bottom of a (a carries the legend on top); both
  # letters level with the legend, above b's y title, as in the manuscript
  paper_assemble(file, width = 7.09, height = a$height, items = list(
    paper_place(a, 0, 0,                         "a", 0, 0.25),
    paper_place(b, a$width, a$height - b$height, "b", a$width, 0.25)))
}

paper_ed5 <- function(a, b, c, file = "Figs/ExtDataFig5/ExtDataFig5_full.pdf"){
  paper_assemble(file, width = 7.09, height = max(a$height, b$height, c$height), items = list(
    paper_place(a, 0, 0,                 "a", 0, 0),
    paper_place(b, a$width, 0,           "b", a$width, 0),
    paper_place(c, a$width + b$width, 0, "c", a$width + b$width, 0)))
}

paper_ed10 <- function(a, b, c, file = "Figs/ExtDataFig10/ExtDataFig10_full.pdf"){
  paper_assemble(file, width = 7.09, height = max(a$height, b$height, c$height), items = list(
    paper_place(a, 0, 0,                 "a", 0, 0),
    paper_place(b, a$width, 0,           "b", a$width, 0),
    paper_place(c, a$width + b$width, 0, "c", a$width + b$width, 0)))
}
