# Copy back defective viral genome (cbDVG) analysis 
Here we provide codes and parameters used in the analysis of measles virus copy back defective viral genome in SSPE brain in the paper: 

# A measles virus collective infectious unit that caused lethal human brain disease includes many locally restricted and few widespread copy-back defective genomes 

di_tector.sh : slurm scheduler setup and parameters used during the DI-tector analysis

combine_results.sh : This script organaizes and sumerize DI-tector output data from each specimen 

DItector_downstream_analysisV3.Rmd :  This script perform the downstream analysis and filtering of cbDVGs with rule of 6 and more than 100 copies, and scatter plot fig 4B-H

plots.py : Scripts for volcano plot, GO biological process heatmap (fig 6 A and B), and a script that generate sequnces for WeLogo analysis (Fig 4I)

cbDVG_forONT_pipeline.sh: a new method for the detection of cbDVGs in Oxford Nanopore Technologies (ONT) sequencing Fastq data.





 
