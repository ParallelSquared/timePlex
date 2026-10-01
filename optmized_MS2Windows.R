
install.packages("mzR")
install.packages("dplyr")
library(mzR)
library(dplyr)
library(ggplot2)


# this code takes a mzML file as input, reads the MS1 intensities across m/z range, and defines MS2 windows that evenly split intentity across the user specified number of windows. 


############## user input ##############
file <- "/path/to/data.mzML" #path to mzML
bin_width <- 0.5 # Th bin width for intensities
x <- 12  # number of bins
starting_mz <- 478 #starting m/z


mzml <- openMSfile(file)

header <- header(mzml)
ms1_indices <- which(header$msLevel == 1)

# function to bin by m/z and sum intensities
bin_and_sum <- function(scan) {
  data <- peaks(mzml, scan)
  df <- data.frame(mz = data[, 1], intensity = data[, 2])
  
  df <- df %>%
    dplyr::mutate(bin = floor(mz / bin_width) * bin_width) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(total_intensity = sum(intensity)) %>%
    ungroup()
  
  return(df)
}

result <- do.call(rbind, lapply(ms1_indices, bin_and_sum))

result <- result %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(total_intensity = sum(total_intensity)) %>%
  ungroup()


ggplot(result, aes(x=bin,y=log10(total_intensity))) + geom_point(alpha=0.5) + geom_smooth() + 
  theme_classic() + labs(x="m/z",y="Log10, Intensity", title="mTRAQ: summed intensities across gradient, binned by m/z")
ggsave("mTRAQ_380_TIC.png",width=5,height=3.7)



result <- result[result$bin>starting_mz,] 

# running sum of total intensity across bins
binned_data <- result %>%
  arrange(bin) %>%
  dplyr::mutate(cumulative_intensity = cumsum(total_intensity))

total_intensity <- sum(binned_data$total_intensity)
target_intensity <- total_intensity / x # amount of intensity each bin should get

bin_edges <- c(1)  #

for (i in 1:x) {
  target <- target_intensity * i
  edge <- which(binned_data$cumulative_intensity >= target)[1]
  bin_edges <- c(bin_edges, edge)
}

bin_edges <- unique(bin_edges)

# find the start and end of each bin
binned_data <- binned_data %>%
  dplyr::mutate(variable_bin = cut(cumulative_intensity, 
                                   breaks = c(binned_data$cumulative_intensity[bin_edges], Inf),
                                   labels = FALSE, 
                                   include.lowest = TRUE))

result_mzs <- binned_data %>%
  group_by(variable_bin) %>%
  summarise(
    start_mz = min(bin) - 0.5,
    end_mz = max(bin) + 0.5,
    center_mz = mean(bin),
    mz_window = end_mz-start_mz,
    total_intensity = sum(total_intensity)
  ) %>%
  ungroup()

result_mzs <- data.frame(result_mzs)

ggplot(result, aes(x=bin, y=log10(total_intensity))) +
  geom_point(alpha=0.5) +
  geom_smooth() +
  geom_rect(
    data = result_mzs,
    aes(xmin = start_mz, xmax = end_mz, ymin = -Inf, ymax = Inf),
    fill = "orange", alpha = 0.1,  
    color = "black",
    size = 0.5,       
    inherit.aes = FALSE  
    
  ) +
  theme_classic() +
  labs(
    x = "m/z",
    y = "Log10, Intensity",
    title = "mTRAQ: summed intensities across gradient, binned by m/z"
  )
ggsave("mTRAQ_",starting_mz,"_TIC_",x,"windows.png",width=5,height=3.7)

write.table(result_mzs, paste0("mTRAQ_",starting_mz,"_1356_",x,"windows_binned_intensities.tsv",sep="\t"))


