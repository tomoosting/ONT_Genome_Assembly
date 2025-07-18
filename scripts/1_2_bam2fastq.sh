#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G
#SBATCH --partition=parallel
#SBATCH --time=3-0:00
#SBATCH --job-name=bam2fastq
#SBATCH -o /nfs/scratch/oostinto/stdout/bam2fastq.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/bam2fastq.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#for mem request about 2x the size of your sam folder

#load modules
module load Python/3.8.2
module load Java/17.0.4

#binaries
export PATH=/nfs/home/oostinto/bin/samtools-1.18:$PATH
export PATH=/nfs/home/oostinto/bin/htslib-1.18:$PATH
export PATH=/nfs/home/oostinto/.local/bin:$PATH

#scripts
#porechop=/nfs/home/oostinto/bin/Porechop/porechop-runner.py

#variables
PROJECT=$1
LIB=$PROJECT'_'$2

#paths
READ_DIR=$SCRATCH/projects/$PROJECT/raw_data/dorado/$LIB
EXT=$READ_DIR/$LIB

#merge bam files
samtools merge --threads 20 -o $EXT'_5mCG_5hmCG.bam' $READ_DIR/bam/*.bam

#bam2fastq
samtools fastq --threads 20 $EXT'_5mCG_5hmCG.bam' | bgzip -c > $EXT.fastq.gz

##filter adapters
echo "removing adapters with porechop"
porechop -i $EXT.fastq.gz -o $EXT'_simplex_porechop.fastq.gz'

#quality filtering
echo "starting filtering with chopper"
bgzip -cd $EXT'_simplex_porechop.fastq.gz' | chopper -q 8  -t 20 -l 100 --headcrop 50 --tailcrop 20 | bgzip -c > $EXT'_simplex_porechop_q8_hc50_tc20_l100.fastq.gz'

#clean
CLEAN_DIR=$READ_DIR/clean
nextflow run rki-mf1/clean -r v1.1.2 --input_type nano                                                  \
                                     --input $EXT'_simplex_porechop_q8_hc50_tc20_l100.fastq.gz' 	    \
							         --control dcs 								                        \
							         --output $CLEAN_DIR             		                            \
                                     --cores 20                                                         \
                                     --memory '200.GB'                                                  \
							         -profile conda
mv $CLEAN_DIR/clean/*.fastq.gz  $EXT'_simplex_porechop_q8_hc50_tc20_l100_clean.fastq.gz'
mv $CLEAN_DIR/qc/multiqc_report.html $EXT'_clean_multiqc_report.html'
#rm -r $CLEAN_DIR

#NanoPlot
source ~/python/NanoPlot_env/bin/activate
    mkdir $READ_DIR/nanoplot
    echo "starting NanoPlot" 
    #NanoPlot -t 10 --fastq $EXT.fastq.gz                                            -o $READ_DIR/nanoplot/$LIB'_raw'
    NanoPlot -t 10 --fastq $EXT'_simplex_porechop_q8_hc50_tc20_l100_clean.fastq.gz' -o $READ_DIR/nanoplot/$LIB'_final'
deactivate