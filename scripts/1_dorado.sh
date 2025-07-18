#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --exclude=gpu01
#SBATCH --time=0-24:00:00
#SBATCH --job-name=dorado
#SBATCH -o /nfs/scratch/oostinto/stdout/dorado.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/dorado.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#1-1269
# Code to set up GPU
module use /home/software/tools/eb_modulefiles/all/Core
module unuse /home/software/tools/modulefiles
module load GCC/10.2.0
module load CUDA/11.1.1
module load OpenMPI/4.0.5

#program
export PATH=~/bin/dorado/dorado-0.9.0-linux-x64/bin:$PATH
export PATH=/nfs/home/oostinto/bin/samtools-1.18:$PATH
export PATH=/nfs/home/oostinto/bin/htslib-1.18:$PATH
export PATH=/home/software/apps/minimap2/2.16:$PATH

#variables
PROJECT=$1
LIB=$2
FILE=$3

#paths
POD5_DIR=$SCRATCH/projects/$PROJECT/raw_data/pod5/$PROJECT'_'$LIB
BAM_DIR=$SCRATCH/projects/$PROJECT/raw_data/dorado/$PROJECT'_'$LIB/bam

model_dir="/nfs/scratch/oostinto/scripts/genome_assembly/1_base_calling/dorado_models"

#create directories
mkdir -p $BAM_DIR

while read line 
do 
    #run dorado (in simplex)
    echo "running dorado basecaller on $line"
    NAME=$( echo $line  | sed -e 's/\.pod5$//' | xargs -n 1 basename )
    dorado basecaller --models-directory $model_dir --trim all sup,5mCG_5hmCG $line | samtools view --bam - > $BAM_DIR/$NAME'_5mCG_5hmCG.bam'
done < $FILE

