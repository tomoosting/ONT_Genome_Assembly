#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --partition=parallel
#SBATCH --time=2-0:00:00
#SBATCH --job-name=purge_haplotigs
#SBATCH -o /nfs/scratch/oostinto/stdout/purge_haplotigs.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/purge_haplotigs.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#binaries
export PATH=/home/software/apps/minimap2/2.16:$PATH
export PATH=/nfs/home/oostinto/bin/samtools-1.18:$PATH
export PATH=/nfs/home/oostinto/bin:$PATH

#parameters
PROJECT=$1
VERSION=$2
ASSEMBLER=$3
ASSEMBLY=$PROJECT'_'$VERSION

#PATHS
READ_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY
MEDAKA_DIR=$READ_DIR/$ASSEMBLER/2_medaka
PURGE_DIR=$READ_DIR/$ASSEMBLER/3_purge_haplotigs
ORDER_DIR=$READ_DIR/$ASSEMBLER/4_reordered
STATS_DIR=$SCRATCH/scripts/genome_assembly/4_assembly_stats
mkdir $PURGE_DIR

#map reads to assembly
#minimap2 -ax map-ont $MEDAKA_DIR/$ASSEMBLY'_polished_'$ASSEMBLER'_assembly.fasta' $READ_DIR/$ASSEMBLY.fastq.gz | samtools view -b | samtools sort > $PURGE_DIR/$ASSEMBLY.bam

#activate conda environment
source activate /nfs/scratch/oostinto/conda/purge_haplotypes
cd $PURGE_DIR

#create histogram
#module load R/4.0.0 
#purge_haplotigs hist -b $PURGE_DIR/$ASSEMBLY.bam -g $MEDAKA_DIR/$ASSEMBLY'_polished_'$ASSEMBLER'_assembly.fasta' -t 10 -d 300

#!!! manually check the histogram to set low (-l), mid (-m), and high (-h) !!!#

#contig coverage
purge_haplotigs contigcov -i $PURGE_DIR/$ASSEMBLY.bam.gencov  -l 48  -m 96  -h 200  [-o $PURGE_DIR/$ASSEMBLY'_coverage_stats.csv' -j 80  -s 80 ]
mv coverage_stats.csv $PURGE_DIR/$ASSEMBLY'_coverage_stats.csv'

#purge
purge_haplotigs purge 	-g $MEDAKA_DIR/$ASSEMBLY'_polished_'$ASSEMBLER'_assembly.fasta'	\
						-c $PURGE_DIR/$ASSEMBLY'_coverage_stats.csv' 					\
						-b $PURGE_DIR/$ASSEMBLY.bam										\
						-o $PURGE_DIR/$ASSEMBLY'_purged_'$ASSEMBLER'_assembly'			\
						-t 10

#reorder and rename contigs
mkdir $ORDER_DIR
seqkit sort --by-length --reverse $PURGE_DIR/$ASSEMBLY'_purged_'$ASSEMBLER'_assembly.fasta' | seqkit replace --pattern '.+' --replacement 'Contig{nr}' > $ORDER_DIR/$ASSEMBLY'_reordered_'$ASSEMBLER'_assembly.fasta'
samtools faidx $ORDER_DIR/$ASSEMBLY'_reordered_'$ASSEMBLER'_assembly.fasta'

#clean up
rm -r $PURGE_DIR/tmp_purge_haplotigs $PURGE_DIR/Rplot001.png $PURGE_DIR/$ASSEMBLY.bam*

#run circleplot and BUSCO
cd $STATS_DIR
sbatch circle_plot.sh $PROJECT $VERSION $ASSEMBLER 4_reordered
sbatch busco_singularity.sh $PROJECT $VERSION $ASSEMBLER 4_reordered