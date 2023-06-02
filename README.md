# Multiscale Embedded Gene Network Analysis (MEGENA)

This repositiry contains code to run a Multiscale Embedded Gene Co-Expression Network using the 'MEGENA' package (Song et al., 2015).

MEGENA revelas differnetually expressed modules of genes.

There are *three* analyses in this repository
1) MEGENA Analaysis
2) Plotting Sunbursts of Modules
3) Plotting Module Networks

There are *three* files in this repository
1) MEGENA_plotting (plotting_megena.Rmd)
2) MEGENA (megena.R)
3) VMEGENA.sh (megena.sh)

Note: 
- Due to MEGENA being a computationally heavy analysis this script is a R file, prepared for workload manager submission. 
- The shell file (megena.sh) will run the MEGENA on a slurm workspace manager using 24 cores.
