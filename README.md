## Simulation and Analysis Demo

This repository contains MATLAB scripts that can be used as templates for:
1) simulation of post-synaptic currents (PSC) rising from detection of co-packaging and independent release of glutamate and GABA neurotransmitters
2) analysis of simulated PSC dataset
3) analysis of electrophysiology experimental PSC dataset 

Please see the publication for detailed description of analysis workflow (https://www.biorxiv.org/content/10.1101/2021.03.23.436594v1)

## Citation info:
Kim SA, Wallace ML, El-Rifai M, Knudsen AR, Sabatini BL (2021). Biophysical demonstration of co-packaging of glutamate and GABA in individual synaptic vesicles in the central nervous system. bioRxiv (preprint)

## MATLAB package requirements:
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox

Add ons: 
- chenxinfeng4/scalebar
- hline and vline

## How to use the scripts:
### 1. Download all MATLAB functions under ```functions/```

### 2. Simulation
   - Open **Fig2_SimulationAnalysis.m**.
   - Run the first section, which shows how to set parameters and perform simulation using **runSimCorelease** function. 

### 3. Analysis
   ![schematic_lt](../master/Images/AnalysisPipelineIm.png)
#### Reproduce figures 2 and 5 from Kim et al. 2021
  #### Figure 2: 
  - Run the full script **Fig2_SimulationAnalysis.m** 
  #### Figure 5:
  - Download raw data under ```exampleData/``` folder in the repo. There are two example datasets, one for a co-packaging and the other for an independent site.
  - Go to the folder containing raw files with AD0_***.mat generic file name.
  - Open and run the full script **Fig5_ExperimentalDataAnalysis.m**, directory path needs to match the location of the files to automatically load individual files

## Where to get more help:
Contact Seul Ah Kim - Graduate student in Sabatini Lab at Harvard Medical School (seulah_kim@g.harvard.edu)

## Lincense
This project is under the MIT License
