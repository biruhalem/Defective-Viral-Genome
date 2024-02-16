#!/bin/bash
#SBATCH --mail-user=beyene.biruhalem@mayo.edu
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=di_test                              #- Also
#SBATCH --partition=cpu-med                           #- Queue also
#SBATCH --nodes=1                                       #- Number of nodes (optional)
#SBATCH --tasks=1                                       #- --ntasks number of cores
#SBATCH --time=1-00:00:00                                 #- time requested
#SBATCH --mem=64G                                        #- Memory/node (or mem-per-cpu)
#SBATCH --export==ALL
#SBATCH --chdir /mforge/research/labs/molecmed/cattaneor/m243773/biruh/Fastq/SSPE/sgDI-tector                 #- Also -D
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error  logs/%x.%N.%j.stderr
#SBATCH --signal=USR1@60                                #- Notify <signal>@<time>

DI=/mforge/research/labs/molecmed/cattaneor/m243773/biruh/Fastq/SSPE/sgDI-tector
Input_Data=/mforge/research/labs/molecmed/cattaneor/m243773/biruh/Fastq/SSPE/SSPE_fastp/fastp_UBS_R1.fastq
Host_Ref=/mforge/research/labs/molecmed/cattaneor/m243773/biruh/Fastq/SSPE/SSPE_fastp/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
Virus_Ref=/mforge/research/labs/molecmed/cattaneor/m243773/biruh/Fastq/SSPE/SSPE_fastp/virus_ref/MVSSPE_ancestor.fa.fasta
OUTDIR=/mforge/research/labs/molecmed/cattaneor/m243773/biruh/Fastq/SSPE/sgDI-tector/DI_results

module load python/3.10.7
module load bwa/0.7.17
module load samtools/1.16

python $DI/DI-tector_06.py -g $Host_Ref -o $OUTDIR -t UBS_2_DI -p 1 -x 16 -k -d $Virus_Ref $Input_Data

