# timePlex (Derks et al. 2026) analysis pipeline

Run with `targets::tar_make()` from this folder.

- `R/` functions
- `pipelines/` targets, one file per manuscript figure
- `Figs/` output: one folder per figure, with each panel and the assembled `*_full.pdf`

## Figures

Each `R/*_functions.R` builds panels; `R/figure_assembly.R` places them on a 180 mm page.
Text sizes per figure are in `paper_fig_sizes` (`R/paper_theme.R`), measured from the manuscript.
PDFs use R's `pdf()` device with Helvetica and Windows-Greek encoding.
Dense point clouds are embedded as 600-dpi images (`paper_raster()`).
