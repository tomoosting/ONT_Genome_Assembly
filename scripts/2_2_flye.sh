#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=50G
#SBATCH --partition=bigmem
#SBATCH --time=10-0:00
#SBATCH --job-name=Flye
#SBATCH -o /nfs/scratch/oostinto/stdout/Flye.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/Flye.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#program
#flye=~/bin/Flye/bin/flye
export PATH=~/bin/Flye/bin/:$PATH

#set parameters
PROJECT=$1
VERSION=$2
ASSEMBLY=$PROJECT'_'$VERSION
#SIZEEST=0.65g #johndory, kahawai
#SIZEEST=0.60g #moki
SIZEEST=0.75g #hapuka
#SIZEEST=0.35g #coral

#paths
READ_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY
FLYE_DIR=$READ_DIR/flye/1_flye
STATS_DIR=$SCRATCH/scripts/genome_assembly/4_assembly_stats
mkdir -p $FLYE_DIR 

#run flye
echo "performing assembly for $ASSEMBLY, using $READ_DIR/$ASSEMBLY.fastq.gz"
flye 	--nano-hq $READ_DIR/$ASSEMBLY.fastq.gz	\
		--out-dir $FLYE_DIR						\
		--genome-size $SIZEEST					\
		--threads 10							\
		--read-error 0.03						\
		--min-overlap 1000						\
		--no-alt-contigs						\
		--scaffold								

#create QC plot for used reads
source ~/python/NanoPlot_env/bin/activate
	NanoPlot -t 10 --fastq $READ_DIR/$ASSEMBLY.fastq.gz -o $READ_DIR/0_NanoPlot
deactivate

#rename output
rename 'assembly' $ASSEMBLY'_raw_flye_assembly' $FLYE_DIR/*

#clean up - never use anything from the subdirectories, so I remove them
rm -r $FLYE_DIR/*/

#run circleplot and BUSCO
cd $STATS_DIR
sbatch circle_plot.sh $PROJECT $VERSION $ASSEMBLER 1_flye
sbatch busco_singularity.sh $PROJECT $VERSION $ASSEMBLER 1_flye
