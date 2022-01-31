# Maths Honours Code

## Repository Structure
All .m files within the main repository are the main MATLAB code files. The following is a list of all folders within the repository

**ArchivedCode:** Folder with old MATLAB code files using old versions of models.

**cbrewer:** [MATLAB colorbrewer schemes](https://au.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) - categorical colour schemes.

**Colormaps:** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps) - colormaps of viridis, magma etc. 

**DataSources:** Data used in the model, and data created from model for supporting information.

**DataSources/archive:** Old data sources that were obtained but not used.

**IdentifyKeySources:** Connectivity matrices data from [Bode et al. (2012)](https://www.int-res.com/abstracts/meps/v466/p155-166) and available [here](https://github.com/MikeBode/IdentifyKeySources).

**Plots:** All plots created using model.

**Plots/00_archive:** Any old plots no longer used.

**Plots/01_progress_report:** Plots created for Honours Sem 1 Progress Report.

**Plots/02_thesis:** Plots created for Honours Thesis.

**Plots/03_paper:** Plots created for writing of paper.

## MATLAB Code Files

### Functions
#### match_reefs_to_areas.m
This function takes the intersection of reefs from the connectivity matrices in `IdentifyKeySources` to the reefs with areas in `DataSources/ReefGazette`, and returns information about the new reef dataset: the number of reefs, the longitude and latitude coordinates, the reef area, and the connectivity matrix. This function is called in both scripts at the beginning to setup the reef structure.

#### calculate_cots_age_1_reproduction.m
This function calculates and returns the proportion of age 1 COTS that can reproduce, or the parameter mu_s in the model, using equations from two papers: [Lucas (1984)]() and [Babcock et al. (2016)](). This function is called in both scripts to calculate this parameter.

#### simulate_reefs_v3.m
This function 

#### calculate_population_box.m

### Scripts

#### compare_scenarios_v3.m

#### sensitivity_analysis_v2.m

