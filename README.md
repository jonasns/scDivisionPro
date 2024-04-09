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

## Quick start

Load the package

```r
library(scDivisionPro)
```

Load you SingleCellExperiment (sce) file

```r
sce <- readRDS("~/Destop/240321_sce_scDivisionPro_testing.rds")
```

If you do not have one, the sce file can be generated from fcs files using the workflow from CATALYST (https://bioconductor.org/packages/devel/bioc/vignettes/CATALYST/inst/doc/preprocessing.html). In short:

Define working directory

```{r setup}
###This need to be set to where the FCS files are
knitr::opts_knit$set(root.dir = '~/Desktop/all_fcs/renamed/')
```

Import fcs files as a flowset

```{r}
fcs_files <- list.files(pattern = ".fcs$")
fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)
```

Import metadata. Metadata file should be in the in working directory. Needs to have sample_id, condition, and patient_id. It is possible to add even more if desired

```{r}
md <- read_excel("Totalcell_metadata.xlsx")                                  
md
```

Import panel data

```{r}
panel <- "Totalcell_panel.xlsx"                           
panel <- read_excel(panel)
panel 
```

Build the SingleCellExperiment from the fcs files, metadata, and panel information

```{r}
sce <- prepData(fs, panel, md)
```

### 1. generate_division_peak_graphs

This function generates graphs that recapitulate division peaks per donor.

```{r}
division_peaks_graphs_list <- generate_division_peak_graphs(sce, moi = "CFSE", soi = c("CD8_CM", "CD8_EM", "Th1", "eTreg"), conditions = c("stim_0x_Treg", "stim_10x_eTreg"), patient_id_col = "patient_id", subsets_col = "subsets", condition_col = "condition")
division_peaks_graphs_list
```

### 2. find_peak_location

This function helps to visualize the approximate location of peaks in CFSE data.

From the 'division_peaks_graphs_list' you can get a first estimate of the location of the beginning dip (bd). It may be different from each donor, so please specify the donor number (dn). Keep iterating this function untill you have values that you are satisfied with.

Additional tweaks can be done to the peak_boundary_uncertainty and the approximate_peak_distance. bd defines the center of the area where the local minimum is calculated, and peak_boundary_uncertainty extends this area. E.g. 4.8 +/- 0.2. 

```r
find_peak_location(sce, donor_list, division_peaks_graphs_list, bd = 4.8, dn = 2, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, soi = "CD8_CM", condition = "stim_0x_Treg", moi = "CFSE")
```

### 3. scDivPro

This function annotates each cell subset by division number and generates accompanying graphs.

Using the 'bd' for each 'dn' in section 2, we can now run the full division annotation on the entire dataset:

```r
sce_d1_results = scDivPro(sce, donor_list, division_peaks_graphs_list, bd = 4.8, dn = 1, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, condition = "test", moi = "CFSE", moi2 = "Ki67", soi_list = c("CD8_CM", "CD8_EM", "Th1", "eTreg"))
```

```r
sce_d2_results = scDivPro(sce, donor_list, division_peaks_graphs_list, bd = 4.5, dn = 2, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, condition = "test", moi = "CFSE", moi2 = "Ki67", soi_list = c("CD8_CM", "CD8_EM", "Th1", "eTreg"))
```

```r
sce_d3_results = scDivPro(sce, donor_list, division_peaks_graphs_list, bd = 4.5, dn = 3, num_peaks = 5, peak_boundary_uncertainty = 0.2, approximate_peak_distance = 1, condition = "test", moi = "CFSE", moi2 = "Ki67", soi_list = c("CD8_CM", "CD8_EM", "Th1", "eTreg"))
```

Make a pdf of the division profiling from each donor

```{r}
pdf(paste0("~/Desktop/CFSE_peaks_CD8_CM_all_donors.pdf"),width=8,height=5,paper='special') 
plot_grid(sce_d1_results[[1]], sce_d2_results[[1]], sce_d3_results[[1]])
dev.off()
```

Remerge the sce's

```{r}
sce_recombined = cbind(sce_d1,sce_d2,sce_d3, deparse.level=1)
sce_recombined$cluster_id = 1
```


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
