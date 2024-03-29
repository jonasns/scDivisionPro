% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_peak_location.R
\name{find_peak_location}
\alias{find_peak_location}
\title{Find the approximate location of peaks}
\usage{
find_peak_location(
  sce,
  donor_list,
  division_peaks_graphs_list,
  bd,
  dn,
  num_peaks = 5,
  peak_boundary_uncertainty = 0.2,
  approximate_peak_distance = 1,
  soi = "CD8_CM",
  condition = "stim_0x_Treg",
  moi = "CFSE"
)
}
\arguments{
\item{sce}{SingleCellExperiment object containing recombined data}

\item{donor_list}{Vector containing donor identifiers}

\item{division_peaks_graphs_list}{List of ggplot objects}

\item{bd}{Beginning dip, as in the dip between division 0 and division 1}

\item{dn}{Donor number}

\item{num_peaks}{Number of peaks to find (default is 5)}

\item{peak_boundary_uncertainty}{The uncertainty for peak boundary approximation (default is 0.2)}

\item{approximate_peak_distance}{The approximate distance between peaks (default is 1)}

\item{soi}{Subset of interest}

\item{condition}{Condition of interest}

\item{moi}{Marker of interest (default is "CFSE")}
}
\description{
This function helps to visualize the approximate location of peaks in CFSE data.
}
