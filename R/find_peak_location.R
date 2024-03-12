#' Find the approximate location of peaks
#' 
#' This function helps to visualize the approximate location of peaks in CFSE data.
#' 
#' @param sce SingleCellExperiment object containing recombined data
#' @param donor_list Vector containing donor identifiers
#' @param division_peaks_graphs_list List of ggplot objects
#' @param bd Beginning dip, as in the dip between division 0 and division 1
#' @param dn Donor number
#' @param num_peaks Number of peaks to find (default is 5)
#' @param peak_boundary_uncertainty The uncertainty for peak boundary approximation (default is 0.2)
#' @param approximate_peak_distance The approximate distance between peaks (default is 1)
#' @param soi Subset of interest
#' @param condition Condition of interest
#' @param moi Marker of interest (default is "CFSE")
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom viridis viridis
#' @export
find_peak_location <- function(sce, donor_list, division_peaks_graphs_list, bd, dn, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, soi = "CD8_CM", condition = "stim_0x_Treg", moi = "CFSE") {
  
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(viridis)
  
  # subset to donor
  donor <- donor_list[dn]
  sce_temp1 <- filterSCE(sce, patient_id == donor)
  
  plot.data <- ggplot_build(division_peaks_graphs_list[[dn]])
  plot.data <- as.data.frame(plot.data$data)
  
  # Find peaks
  plot.data2 <- plot.data[plot.data$x > bd - peak_boundary_uncertainty & plot.data$x < bd + peak_boundary_uncertainty,]
  intersect1 <- plot.data2$x[plot.data2$density == min(plot.data2$density)]
  
  plot.data2 <- plot.data[plot.data$x > bd - (approximate_peak_distance + peak_boundary_uncertainty) & plot.data$x < bd - (approximate_peak_distance - peak_boundary_uncertainty),]
  intersect2 <- plot.data2$x[plot.data2$density == min(plot.data2$density)]
  
  plot.data2 <- plot.data[plot.data$x > bd - (2 * approximate_peak_distance + peak_boundary_uncertainty) & plot.data$x < bd - (2 * approximate_peak_distance - peak_boundary_uncertainty),]
  intersect3 <- plot.data2$x[plot.data2$density == min(plot.data2$density)]

  first_graph <- division_peaks_graphs_list[[dn]] + geom_vline(xintercept = c(intersect1, intersect2, intersect3))
  
  ave_ints <- ((intersect1 - intersect2) + (intersect2 - intersect3)) / 2
  
  second_graph <- division_peaks_graphs_list[[dn]] + geom_vline(xintercept = c(intersect1 + ave_ints, intersect1, intersect1 - 1 * ave_ints, intersect1 - 2 * ave_ints, intersect1 - 3 * ave_ints, intersect1 - 4 * ave_ints))
  
  # Calculate the differences for all conditions
  differences <- intersect1 + ave_ints - (0:(num_peaks - 1)) * ave_ints

  # extend the limit of the highest CFSE peak
  differences[1] <- differences[1] + 2
  
  # Sort the differences in non-decreasing order
  sorted_differences <- sort(differences)
  
  # Find the interval indices
  div_levels <- findInterval(sce_temp1@assays@data$exprs[moi, ], sorted_differences, rightmost.closed = TRUE)
  
  # Calculate the maximum division level
  max_div_level <- length(sorted_differences)
  
  # Convert division levels to start at "div0"
  converted_div_levels <- max_div_level - div_levels
  
  # Adjust the division levels to match the desired range ("div0" to num_peaks)
  adjusted_div_levels <- converted_div_levels - 1
  
  # Assign the adjusted division levels based on the comparisons
  sce_temp1$div <- paste0("div", adjusted_div_levels)
  
  color_intervals <- differences
  line_intervals <- color_intervals
  color_intervals[num_peaks + 1] <- -1
  
  sce_use <- filterSCE(sce_temp1, condition == condition & subsets == soi)
  
  new_tab <- data.frame(exprs = sce_use@assays@data$exprs[moi, ], div = sce_use$div, sample_id = sce_use$sample_id, condition = sce_use$condition, patient_id = sce_use$patient_id, subsets = sce_use$subsets)
  new_tab <- melt(new_tab)
  
  p <- ggplot(new_tab, aes(x = value, y = subsets)) + 
    geom_density_ridges(aes(x = value, y = subsets)) +
    theme_bw(base_size = 14) +
    theme_ridges()
  
  plot.data <- ggplot_build(p)
  
  new_tab2 <- data.frame(x = plot.data$data[[1]]$x, density = plot.data$data[[1]]$density)
  
  third_graph <- ggplot(new_tab2, aes(x = x, y = density)) + 
    geom_area(aes(fill = cut(x, breaks = color_intervals)), color = "black") +
    scale_fill_manual(values = rev(viridis(num_peaks))) +
    geom_vline(xintercept = line_intervals, color = "black", linetype = "dashed") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    ggtitle(paste0(donor, "_", soi))
  
  return(list(first_graph = first_graph, second_graph = second_graph, third_graph = third_graph))
}

