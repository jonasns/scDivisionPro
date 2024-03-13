# scDivisionPro

scDivisionPro is an R package designed for single-cell division profiling analysis. It offers functions to generate division peaks recap graphs, find the approximate location of peaks, and perform single-cell division profiling.

## Installation

You can install the scDivisionPro package from GitHub using the `devtools` package:

```r
devtools::install_github("jonasns/scDivisionPro")
```

## Dependencies

please install the following dependencies separatenly according to their instructions:

ggplot2

reshape2

ggridges

viridis

CATALYST

## Functions

### 1. generate_division_peak_graphs

This function generates graphs that recapitulate division peaks per donor.

```r
generate_division_peak_graphs(sce, moi = "CFSE", soi = c("CD8_CM", "CD8_EM", "Th1", "eTreg"), conditions = c("test", "ctrl"), patient_id_col = "patient_id", subsets_col = "subsets", condition_col = "condition")
```

#### parameters

**sce** A SingleCellExperiment object containing the data.

**moi** Marker of interest.

**patient_id_col** Name of the column containing patient IDs.

**subsets_col** Name of the column containing subset information.

**soi** Vector of subsets of interest names to filter data.

**condition_col** Name of the column containing condition information.

**conditions** Vector of condition names to filter data.

#### return

A list of graphs recapitulating division peaks per donor.

### 2. find_peak_location

This function helps to visualize the approximate location of peaks in CFSE data.

```r
find_peak_location(sce, donor_list, division_peaks_graphs_list, bd, dn, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, soi = "CD8_CM", condition = "stim_0x_Treg", moi = "CFSE")
```

#### parameters

**sce** SingleCellExperiment object containing recombined data

**donor_list** Vector containing donor identifiers

**division_peaks_graphs_list** List of ggplot objects

**bd** Beginning dip, as in the dip between division 0 and division 1

**dn** Donor number

**num_peaks** Number of peaks to find (default is 5)

**peak_boundary_uncertainty** The uncertainty for peak boundary approximation (default is 0.2)

**approximate_peak_distance** The approximate distance between peaks (default is 1)

**soi** Subset of interest

**condition** Condition of interest

**moi** Marker of interest (default is "CFSE")

#### return

Two graphs showing the division_peak_graph (from function 1) with:

1. Lines indicating where the local dips inbetween CFSE peaks have been found

2. Lines that are equally distanced, calculated as the average of the first 2 peaks

A graph like 2, but for a specific subset of interest (soi) and colourcoded per division.

### 3. scDivPro

This function annotates each cell subset by division number and generates accompanying graphs.

```r
scDivPro(sce, donor_list, division_peaks_graphs_list, bd, dn, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, condition = "test", moi = "CFSE", moi2 = "Ki67", soi_list = c("CD8_CM", "CD8_EM", "Th1", "eTreg"))
```

#### parameters

**sce** SingleCellExperiment object

**donor_list** Vector containing donor identifiers

**division_peaks_graphs_list** List of division peaks graphs

**bd** Beginning dip, as in the dip between division 0 and division 1

**dn** Donor number

**num_peaks** Number of peaks

**peak_boundary_uncertainty** The uncertainty for peak boundary approximation (default is 0.2)

**approximate_peak_distance** The approximate distance between peaks (default is 1)

**condition** Condition of interest

**moi** Marker of interest (default is CFSE)

**moi2** 2nd marker of interest to use for scatter plot (default is Ki67)

**soi_list** A vector containing all subsets of interest.

#### return

A list containing sce annotated for each cell subset by division data and accompanying graphs related to the analysis.

There are three types of graphs: 

- density_plot
  
- violin_plot
  
- scatter_plot


## License

This project is licensed under the MIT License - see the LICENSE.md file for details.
```
