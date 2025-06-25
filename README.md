# nct-mci

This repository contains code used in ***"Disrupted Energy Landscape in Individuals with Mild Cognitive Impairment: Insights from Network Control Theory"*** (https://doi.org/10.1101/2025.06.19.660545).

Analysis was conducted using a combination of R and MATLAB. The R script is designed to be executed sequentially. Within the script, there are instructions on when to run portions of the MATLAB code. MATLAB code block are numbered 1-17, and are referenced by these block numbers in the R script.

Please organize your project folder as follows. In the ```Configuration``` blocks at the top of both scripts, specify the file paths corresponding to each directory. These directories are referenced throughout both scripts. MATLAB packages can be found at https://github.com/singlesp/energy_landscape and https://github.com/kjamison/atlasblobs. 

<pre>
<code>
Project Folder
├── R                              
│   └── NCT.MCI.Github.code.R       # Main R script
├── MATLAB                         
│   └── NCT.MCI.Github.code.m       # Main MATLAB script
│   ├── Packages/                   # MATLAB packages
├── Data/                         
│   ├── Time Series/                # Folder containing time series data (.txt)
│   └── Demographics.RData          
│   └── Tau.csv                     # File containing tau data (.csv)
│   └── Plaque.csv                  # File containing amyloid beta data (.csv)
│   └── SC.mat                      # File containing structural connectome data
├── Results/                      
│   ├── Clusters/                  
│   ├── Figures/                    # Plots/images generated from analysis
</code>
</pre>
