# nct-mci
<pre>
<code>
Project Folder
├── R                              
│   └── NCT.MCI.Github.code.R       # Main R script
├── MATLAB                         
│   └── NCT.MCI.Github.code.m       # Main MATLAB script
│   ├── Packages/                   # MATLAB packages
├── Data/                         
│   └── Demographics.RData        
│   ├── Time Series/                # Folder containing time series data (.txt)
│   └── Tau.csv                     # File containing tau data (.csv)
│   └── Plaque.csv                  # File containing amyloid beta data (.csv)
│   └── SC.mat                      # File containing structural connectome data
├── Results/                      
│   ├── Clusters/                  
│   ├── Figures/                    # Plots/images generated from analysis
├── README.md                       # Project overview and instructions
</code>
</pre>


#Demographic.Data.RData should have the following columns: ID, Age, Sex (female = 1, male = 0), Status (HC or MCI), T1, one column per region containing
