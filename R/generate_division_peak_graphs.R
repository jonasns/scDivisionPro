#' Generate Division Peaks Recap Graphs
#'
#' This function generates graphs that recapitulate division peaks per donor.
#'
#' @param sce A SingleCellExperiment object containing the data.
#' @param moi Marker of interest.
#' @param patient_id_col Name of the column containing patient IDs.
#' @param subsets_col Name of the column containing subset information.
#' @param soi Vector of subsets of interest names to filter data.
#' @param condition_col Name of the column containing condition information.
#' @param conditions Vector of condition names to filter data.
#' @return A list of graphs recapitulating division peaks per donor.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges
#' @export
generate_division_peak_graphs <- function(sce, 
                                          moi = "CFSE",
                                          soi = c("CD8_CM", "CD8_EM", "Th1", "eTreg"),
                                          conditions = c("stim_0x_Treg", "stim_10x_eTreg", "stim_10x_nTreg", "stim_10x_Tfr", "unstim_0x_Treg"),
                                          patient_id_col = "patient_id",
                                          subsets_col = "subsets",
                                          condition_col = "condition") {
  
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(CATALYST)
  
  # Selecting conditions to use
  donor_list <- rownames(table(sce[[patient_id_col]]))
  
  # Assign donor_list to global environment
  assign("donor_list", donor_list, envir = .GlobalEnv)
  
  # Start of loops
  pp_list <- list()
  for (i in seq_along(donor_list)) {
    donor <- donor_list[i]
    nam <- paste0("p", i) 
    
    sce_temp1 <- filterSCE(sce, sce[[subsets_col]] %in% soi & sce[[patient_id_col]] == donor & sce[[condition_col]] %in% conditions)

    test1 <- as.data.frame(t(sce_temp1@assays@data$exprs))
    
    test1 <- test1[[moi]]
    test1 <- melt(test1)
    test1$culture_condition <- 1
    test1$culture_condition <- as.factor(test1$culture_condition)
    
    # Plotting
    p <- ggplot(test1, aes(x = value, y = culture_condition, group = culture_condition)) + 
      geom_density_ridges(fill = "#00AFBB") +
      ylab("% of max") +
      xlab(paste(moi, " expression")) +
      theme_bw(base_size = 14) +
      ggtitle(paste(moi, donor, sep = "_"))
    
    assign(nam, p)
    pp_list[[i]] <- p
  }
  
  return(pp_list)
}

