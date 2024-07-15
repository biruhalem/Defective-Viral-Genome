#---------------------------------------------------------------------------------------
# Spesifications for slurm scheduler on High performance computer cluster (HPC cluster) 
#---------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --mail-user=email-adress
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=DI-tector                              #- Also
#SBATCH --partition=cpu-med                           #- Queue also
#SBATCH --nodes=1                                       #- Number of nodes (optional)
#SBATCH --tasks=1                                       #- --ntasks number of cores
#SBATCH --time=1-00:00:00                                 #- time requested
#SBATCH --mem=64G                                        #- Memory/node (or mem-per-cpu)
#SBATCH --export==ALL
#SBATCH --chdir SSPE/sgDI-tector                 #- Also -D
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error  logs/%x.%N.%j.stderr
#SBATCH --signal=USR1@60                                #- Notify <signal>@<time>

#----------------------------------------
# Defining the directories and file paths
#----------------------------------------

DI=SSPE/sgDI-tector  # This direcory contains the DI-tector_06.py program from Beauclair et al., 2018 (DOI: 10.1261/rna.066910.118).
Input_Data=SSPE/SSPE_fastp/fastp_UBS_R1.fastq  # RNA-seq data raw reads from UBS SSPE specimen. 
Host_Ref=SSPE/SSPE_fastp/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  # Human reference genome
Virus_Ref=SSPE/SSPE_fastp/virus_ref/MVSSPE_ancestor.fa.fasta  #Measles virus reference genome 
OUTDIR=SSPE/sgDI-tector/DI_results #Output directory for the analysis 

#-----------------------------------
# loading modules using slurm script
#-----------------------------------

module load python/3.10.7
module load bwa/0.7.17
module load samtools/1.16

#---------------------
# DI-tector parameters 
#---------------------

python $DI/DI-tector_06.py -g $Host_Ref -o $OUTDIR -t UBS_2_DI -p 1 -x 16 -k -d $Virus_Ref $Input_Data

#-g- host reference genome
#-o - output directory 
#-t - prefix for thr output result
#-p - (+/-) strand genome 
#-x number of threads 
#-k get all the intermediate files 
#-d get all the DVG sequnces 

