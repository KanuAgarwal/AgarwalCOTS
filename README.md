# Optimising crown-of-thorns starfish control effort on the Great Barrier Reef

## Repository Structure
All .m files within the main repository are the main [MATLAB code files](https://github.com/KanuAgarwal/MathsHonoursCode#matlab-code-files). The following is a list of all folders within the repository

- **ArchivedCode:** Folder with old MATLAB code files using old versions of models.
- **cbrewer:** [MATLAB colorbrewer schemes](https://au.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) - categorical colour schemes.
- **Colormaps:** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps) - colormaps of viridis, magma etc. 
- **DataSources:** Data used in the model, and data created from model for supporting information.
	- **DataSources/archive:** Old data sources that were obtained but not used.
- **IdentifyKeySources:** Connectivity matrices data from [Bode et al. (2012)](https://www.int-res.com/abstracts/meps/v466/p155-166) and available [here](https://github.com/MikeBode/IdentifyKeySources).
- **Plots:** All plots created using model.

## MATLAB Code Files

### Functions
#### match_reefs_to_areas.m
This function takes the intersection of reefs from the connectivity matrices in `IdentifyKeySources` to the reefs with areas in `DataSources/ReefGazette.csv`, and returns information about the new reef dataset: the number of reefs, the longitude and latitude coordinates, the reef area, and the connectivity matrix. This function is called in both [scripts](https://github.com/KanuAgarwal/MathsHonoursCode#scripts) at the beginning to setup the reef structure.

#### calculate_cots_age_1_reproduction.m
This function calculates and returns the proportion of age 1 COTS that can reproduce, or the parameter mu_s in the model, using equations from two papers: [Lucas (1984)](https://www.sciencedirect.com/science/article/abs/pii/0022098184902144) and [Babcock et al. (2016)](https://link.springer.com/article/10.1007/s00227-016-3009-5). This function is called in both [scripts](https://github.com/KanuAgarwal/MathsHonoursCode#scripts) to calculate this parameter.

#### simulate_reefs_v3.m
This function runs a simulation of the GBR model given all of the required parameters, initial conditions, etc. and returns the coral and starfish populations over time. See [file](https://github.com/KanuAgarwal/MathsHonoursCode/blob/main/simulate_reefs_v3.m) for more information on function inputs and outputs. 

#### calculate_population_box.m
This function calculates the population size of coral, and age-structured COTS within the outbreak initiation box over time, given the coral and age-structured COTS population size on the whole GBR over time. This function is called in both [scripts](https://github.com/KanuAgarwal/MathsHonoursCode#scripts) after each simulation.

### Scripts

#### compare_scenarios_v3.m
This script runs all the simulations required to compare all of the control scenarios for the paper, and it creates all plots for the paper (plus some extra plots that aren't used). There are also lots of others plots that were created but have been commented out.

#### sensitivity_analysis_v2.m
This script runs the sensitivity analysis of the initial adult starfish population at each reef at the initiation box, and creates the plot used in the paper. 

## Usage
### Before Running in MATLAB
To ensure the plots are created with the correct colour schemes and the code doesn't throw an error, you may need to right click on the folders `cbrewer` and `Colormaps` in the 'Current Folder' window within MATLAB and select the option 'Add to Path' and then select the option 'Selected Folders and Subfolders'.

### Color Schemes
Most of the colour schemes used in the plots have been taken from [cbrewer](https://au.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) or [Colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps). However some viridis categorical palettes were taken from [this website](https://waldyrious.net/viridis-palette-generator/).
