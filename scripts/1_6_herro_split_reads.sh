#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --mem=30G
#SBATCH --partition=parallel
#SBATCH --time=0-5:00
#SBATCH --job-name=herro_split_reads
#SBATCH -o /nfs/scratch/oostinto/stdout/herro_split_reads.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/herro_split_reads.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#variables
PROJECT=$1
GSIZE=$2 #in Gb e.g. 0.7

# hapuka  = 0.76Gb
# kahawai = 0.66Gb
# moki    = 0.60Gb

#fasta file
FASTA=$SCRATCH/projects/$PROJECT/raw_data/dorado/$PROJECT
NANO=$SCRATCH/projects/$PROJECT/raw_data/dorado/nanoplot

#get stats on herro hifi reads
seqkit stats $FASTA'_herro_hifi.fasta' > $FASTA'_herro_hifi_stats.tsv'

#split herro hifi reads into 30kbp segments
seqkit sliding -j 10 -g -s 30000 -W 30000 $FASTA'_herro_hifi.fasta' | seqkit seq -m 10000 -w 0 > $FASTA'_herro_hifi_split.fasta'
seqkit stats $FASTA'_herro_hifi_split.fasta' > $FASTA'_herro_hifi_split_stats.tsv'

#Nanoplot
source ~/python/NanoPlot_env/bin/activate
    NanoPlot -t 10 --fasta $FASTA'_herro_hifi_split.fasta' -o $NANO/herro_hifi_split
deactivate

#get total number of bp in hifi fasta
echo "estimated genome size is $GSIZE" > $FASTA'_subsample_stats.txt'
BP_TOT=$(tail -n 1 "${FASTA}_herro_hifi_split_stats.tsv" | sed -E 's/[[:space:]]+/ /g' | cut -d' ' -f5 | sed 's/,//g')
echo "total number of bp is $BP_TOT" >> $FASTA'_subsample_stats.txt'

#calculate coverage based on estimated genome size
COV=$(echo "scale=0; $BP_TOT / ($GSIZE * 1000000000)" | bc -l)
echo "estmated coverage is $COV" >> $FASTA'_subsample_stats.txt'

#based on herro paper - subset down to 35x of coverage is above 38x
if [ "$COV" -ge 38 ]; then
    echo "covrage is >= 38" >> $FASTA'_subsample_stats.txt'
    #calculate fraction to subsample for getting cov ~35x
    FRAC=$(echo "scale=2; 35 / $COV" | bc -l)
    echo "fraction used: $FRAC" >> $FASTA'_subsample_stats.txt'
    
    #use output seqkit stats and genome size estimation to decide
    seqkit sample -j 10 -p $FRAC $FASTA'_herro_hifi_split.fasta' > $FASTA'_herro_hifi_split_35x.fasta'
    seqkit stats $FASTA'_herro_hifi_split_35x.fasta' > $FASTA'_herro_hifi_split_35x_stats.tsv'

    #only subsample and NOT split reads for assembly using verkko
    seqkit sample -j 10 -p $FRAC $FASTA'_herro_hifi.fasta' > $FASTA'_herro_hifi_35x.fasta'
    seqkit stats $FASTA'_herro_hifi_35x.fasta' > $FASTA'_herro_hifi_35x_stats.tsv'

    #Nanoplot
    source ~/python/NanoPlot_env/bin/activate
        NanoPlot -t 10 --fasta $FASTA'_herro_hifi_split_35x.fasta' -o $NANO/herro_hifi_split_35x
        NanoPlot -t 10 --fasta $FASTA'_herro_hifi_35x.fasta'       -o $NANO/herro_hifi_35x
    deactivate
fi

