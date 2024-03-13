# scDivisionPro

scDivisionPro is an R package designed for single-cell division profiling analysis. It offers functions to generate division peaks recap graphs, find the approximate location of peaks, and perform single-cell division profiling.

## Installation

You can install the scDivisionPro package from GitHub using the `devtools` package:

```r
devtools::install_github("jonasns/scDivisionPro")
```

## Functions
### 1. generate_division_peak_graphs
This function generates graphs that recapitulate division peaks per donor.

```r
generate_division_peak_graphs(sce, moi = "CFSE", soi = c("CD8_CM", "CD8_EM", "Th1", "eTreg"), conditions = c("stim_0x_Treg", "stim_10x_eTreg", "stim_10x_nTreg", "stim_10x_Tfr", "unstim_0x_Treg"), patient_id_col = "patient_id", subsets_col = "subsets", condition_col = "condition")
```

### 2. find_peak_location
This function helps visualize the approximate location of peaks in CFSE data.

```r
find_peak_location(sce, donor_list, division_peaks_graphs_list, bd, dn, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, soi = "CD8_CM", condition = "stim_0x_Treg", moi = "CFSE")
```

### 3. scDivPro
```r
scDivPro(sce, donor_list, division_peaks_graphs_list, bd, dn, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, condition = "stim_0x_Treg", moi = "CFSE", moi2 = "Ki67", soi_list = c("CD8_CM", "CD8_EM", "Th1", "eTreg"))
```

## Dependencies
ggplot2
reshape2
ggridges
viridis
CATALYST
License

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
```
