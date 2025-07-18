#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --ntasks=10
#SBATCH --mem=300G
#SBATCH --time=2-0:00:00
#SBATCH --job-name=merge&qc
#SBATCH -o /nfs/scratch/oostinto/stdout/merge&qc.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/merge&qc.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#load modules 
module load Python/3.8.2

#binaries
export PATH=/nfs/home/oostinto/bin/htslib-1.18:$PATH
export PATH=/nfs/home/oostinto/bin/samtools-1.18:$PATH

#run parameters
PROJECT=$1

#PATHS
READ_DIR=$SCRATCH/projects/$PROJECT/raw_data/dorado

#merge fastq.gz files
cat $READ_DIR/*/*_simplex_porechop_q8_hc50_tc20_l100_clean.fastq.gz > $READ_DIR/$PROJECT'_simplex_filtered.fastq.gz'

#run NanoPlot on all fastq files
mkdir $READ_DIR/nanoplot
source ~/python/NanoPlot_env/bin/activate
    NanoPlot -t 10 --fastq $READ_DIR/$PROJECT'_simplex_filtered.fastq.gz'    -o $READ_DIR/nanoplot/simplex   
deactivate


