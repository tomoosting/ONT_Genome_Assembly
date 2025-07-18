#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --partition=parallel
#SBATCH --time=3-0:00
#SBATCH --job-name=filter_longreads
#SBATCH -o /nfs/scratch/oostinto/stdout/filter_longreads.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/filter_longreads.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#binaries
export PATH=/nfs/home/oostinto/bin/htslib-1.18:$PATH

#variables
PROJECT=$1

#paths
FASTQ=$SCRATCH/projects/$PROJECT/raw_data/dorado/$PROJECT
NANO=$SCRATCH/projects/$PROJECT/raw_data/dorado/nanoplot

# run chopper to filter reads over 10,20 and 50 kbp
bgzip -cd $FASTQ'_simplex_filtered.fastq.gz'           | chopper -t 10 -l 1000   -q 10 | bgzip -c > $FASTQ'_simplex_filtered_1kbp_q10.fastq.gz'
bgzip -cd $FASTQ'_simplex_filtered_1kbp_q10.fastq.gz'  | chopper -t 10 -l 10000        | bgzip -c > $FASTQ'_simplex_filtered_10kbp_q10.fastq.gz'
bgzip -cd $FASTQ'_simplex_filtered_10kbp_q10.fastq.gz' | chopper -t 10 -l 20000        | bgzip -c > $FASTQ'_simplex_filtered_20kbp_q10.fastq.gz'
bgzip -cd $FASTQ'_simplex_filtered_20kbp_q10.fastq.gz' | chopper -t 10 -l 50000        | bgzip -c > $FASTQ'_simplex_filtered_50kbp_q10.fastq.gz'

# run NanoPlot on all fastq files 
source ~/python/NanoPlot_env/bin/activate
    NanoPlot -t 10 --fastq $FASTQ'_simplex_filtered_1kbp_q10.fastq.gz'  -o $NANO/simplex_1kbp
    NanoPlot -t 10 --fastq $FASTQ'_simplex_filtered_10kbp_q10.fastq.gz' -o $NANO/simplex_10kbp
    NanoPlot -t 10 --fastq $FASTQ'_simplex_filtered_20kbp_q10.fastq.gz' -o $NANO/simplex_20kbp
    NanoPlot -t 10 --fastq $FASTQ'_simplex_filtered_50kbp_q10.fastq.gz' -o $NANO/simplex_50kbp
deactivate


