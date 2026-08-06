# README

Source code for "Hierarchical drivers of ribbed mussel (*Geukensia demissa*) 
population dynamics and spatial distribution across marsh landscapes" submitted
to Diversity and Distributions. Here, we measured ribbed mussel individual 
growth rates, predation, and recruitment across the marsh landscape and examined 
the hierarchical drivers of these vital rates. Additionally, we estimated the 
spatial distribution of population vital rates to compare with spatial 
distribution of mussel densities.

## Contact information and citation

```bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```
Cite as:

> CITE

## Installation

Download and unzip all required data from the 
[Zenodo repository](https://doi.org/10.5281/zenodo.21828495) into the 
working directory. <ins>IMPORTANT:</ins> Retain existing folder structures

<br>
Install all required software and packages and ensure they are the proper 
version:

```text
JAGS Version: 4.3.2
R Version: 4.3.1

R packages:
- 'corrplot' Version: 0.95
- 'dplyr' Version: 1.1.4
- 'ggplot2' Version: 4.0.2
- 'jagsUI' Version: 1.6.2
- 'raster' Version: 3.6.32
- 'readxl' Version: 1.4.5
- 'rjags' Version: 4.17
- 'rstatix' Version: 0.7.2
- 'scales' Version: 1.4.0
- 'spdep' Version: 1.3.11
- 'stringr' Version: 1.5.1

```



## Workflow
R scripts are listed in order of workflow.

### Structural Equation Model Analysis 
**Scripts:** `01_mussel_growth_SEM`, `02_mussel_predation_SEM`,`03_mussel_recruitment_SEM`

<ins>Purpose:</ins> Conduct structural equation modelling (SEM) of ribbed 
mussel population vital rates using the SEM methodology of 
[Grace et al., 2012.]( https://doi.org/10.1890/ES12-00048.1) 
This script exports JAGs model code to be ran using Rjags and jagsUI. 
Standardized effect queries are calculated for comparison of both direct and 
indirect effect sizes.

<ins>Input:</ins> `mussel_population_vital_data.xlsx`

<ins>Output:</ins> `model_out`


### Spatial Analysis
**Script:** `04_pop_vital_spatial_analysis.R`

<ins>Purpose:</ins> Creates predicted spatial distributions of ribbed mussel 
growth, predation, and recruitment using SEM outputs and rasters representing 
marsh-level explanatory variables. The predictor rasters are then used in 
correlation analyses with predicted mussel densities from 
[Annis et al., 2022](https://doi.org/10.1007/s12237-022-01090-w), 
where Spearman’s correlation is calculated between all vitals and densities. 
Finally, SEM residuals are tested for spatial autocorretion using Moran's I.

<ins>Input:</ins> `mussel_population_vital_data.xlsx`, `model_out`, 
`Mussel SDM rasters`, `plot_parameters/site_coors.csv`

<ins>Output:</ins> `Prediction rasters`, `plots/cor`, `plots/moran_tables`


### Plotting
**NOTE:** Raw R plots and tables used to create manuscript figures can be found 
in the `plots` directory.

**Script:** `05_summary_plots.R`

<ins>Purpose:</ins> Plots the summary statistics from the three population vital 
rate experiments. 

<ins>Input:</ins> `mussel_population_vital_data.xlsx`

<ins>Output:</ins> `plots/summary`

<br>

**Script:** `06_direct_effect_plots.R`

<ins>Purpose:</ins> Creates direct effect plots for every pathway in each vital 
rate SEM. 

<ins>Input:</ins> `mussel_population_vital_data.xlsx`, `model_out`, 
`plot_parameters`

<ins>Output:</ins> `plots/direct`

<br>

**Script:** `07_total_effect_plots.R`

<ins>Purpose:</ins> Creates total effect predictions of temperature and 
marsh-level predictors on each vital rates, holding all other variables 
constant

<ins>Input:</ins> `mussel_population_vital_data.xlsx`, `model_out`

<ins>Output:</ins> `plots/total_effect`


