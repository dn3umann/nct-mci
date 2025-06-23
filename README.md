# nct-mci
<pre>
<code>
Project Folder
├── R                               # All R code and files
│   └── NCT.MCI.code.R              # Main R script
├── MATLAB                          # All MATLAB code and files
│   └── NCT.MCI.code.m              # Main MATLAB script
│   ├── Packages/                   # MATLAB packages
├── Data/                           # Input data files
│   └── Demographics.RData          # Demographic data
│   ├── Time Series/                # Folder containing time series data (.txt)
│   └── Tau.csv                     # File containing tau data (.csv)
│   └── Plaque.csv                  # File containing amyloid beta data (.csv)
│   └── SC.mat                      # File containing structural connectome data
├── Results/                        # Output results
│   ├── Clusters/                   # Clustering results
│   ├── Figures/                    # Plots/images generated from analysis
├── README.md                       # Project overview and instructions
├── .gitignore                      # Excludes large/private files
└── NCT_Project.Rproj               # RStudio project file
</code>
</pre>


#Demographic.Data.RData should have the following columns: ID, Age, Sex (female = 1, male = 0), Status (HC or MCI), T1, one column per region containing
