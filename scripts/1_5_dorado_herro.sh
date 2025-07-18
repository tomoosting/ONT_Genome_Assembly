#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=200
#SBATCH --mem=500G
#SBATCH --time=1-0:00:00
#SBATCH --job-name=herro
#SBATCH -o /nfs/scratch/oostinto/stdout/herro.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/herro.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

# Code to set up GPU
module use /home/software/tools/eb_modulefiles/all/Core
module unuse /home/software/tools/modulefiles
module load GCC/10.2.0
module load CUDA/11.1.1
module load OpenMPI/4.0.5

#program
export PATH=~/bin/dorado/dorado-0.9.0-linux-x64/bin:$PATH
export PATH=~/bin/htslib-1.18:$PATH

#variables
PROJECT=$1

#fasta file
FASTQ=$SCRATCH/projects/$PROJECT/raw_data/dorado/$PROJECT
NANO=$SCRATCH/projects/$PROJECT/raw_data/dorado/nanoplot

# download model - need to do ones
# dorado download --model herro-v1

# perform correction
dorado correct -t 200 -m herro-v1 $FASTQ'_simplex_filtered_10kbp_q10.fastq.gz' > $FASTQ'_herro_hifi.fasta'

# run NanoPlot on all fastq files 
source ~/python/NanoPlot_env/bin/activate
    NanoPlot -t 10 --fasta $FASTQ'_herro_hifi.fasta' -o $NANO/herro_hifi
deactivate
