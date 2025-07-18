#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --partition=quicktest
#SBATCH --time=0-05:00
#SBATCH --job-name=FCS_adaptor
#SBATCH -o /nfs/scratch/oostinto/stdout/FSC_adaptor.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/FSC_adaptor.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###load modules
#module load GCC/10.2.0 OpenMPI/4.0.5 Singularity/3.10.2

export PATH=~/singularity/bin:$PATH
export SINGULARITY_BIND="/nfs/scratch/oostinto,/nfs/home/oostinto"
export SINGULARITY_TMPDIR="/nfs/scratch/oostinto/tmp"

#PARAMETERS
PROJECT=$1
VERSION=$2
ASSEMBLY=$PROJECT'_'$VERSION
ASSEMBLER=hifiasm

#FCS
FCS_ADAPTOR_SIF=/nfs/scratch/oostinto/singularity/fcs-adaptor.sif
FCS_GX_SIF=/nfs/scratch/oostinto/singularity/fcs-gx.sif
PY=/nfs/scratch/oostinto/scripts/genome_assembly/4_assembly_stats/fcs.py

#READ_DIR
FASTA=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/hifiasm/2_purge_dups/$ASSEMBLY.asm.bp.p_ctg.reordered_purged.fasta
ASM_NAME=$(basename "$FASTA" | sed 's/\.[^.]*$//')
FCS_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/hifiasm/3_fcs
mkdir -p $FCS_DIR/fcs_adaptor

#run fcs_adaptor
./run_fcsadaptor.sh --fasta-input $FASTA                \
                    --output-dir $FCS_DIR/fcs_adaptor   \
                    --euk                               \
                    --container-engine singularity      \
                    --image $FCS_ADAPTOR_SIF

# clean genome
export FCS_DEFAULT_IMAGE=$FCS_GX_SIF
export NCBI_FCS_REPORT_ANALYTICS=1
cat $FASTA | python3 $PY clean genome                                               \
        --action-report $FCS_DIR/fcs_adaptor/fcs_adaptor_report.txt                 \
        --output $FCS_DIR/fcs_adaptor/$ASM_NAME.fcs_adaptor.clean.fasta             \
        --contam-fasta-out $FCS_DIR/fcs_adaptor/$ASM_NAME.fcs_adaptor.contam.fasta        













































#    singularity run $CONTAINER --bind $FASTA_DIRNAME:/sample-volume/ \
#        --bind $EXPANDED_OUTPUT:/output-volume/ $SINGULARITY_IMAGE \
#        /app/fcs/bin/av_screen_x -o /output-volume/ $TAX /sample-volume/$FASTA_FILENAME $DEBUG
# $TAX = --euk $DEBUG = --debug

#usage: av_screen_x.py [-h] [-o path] [--bin-dir BIN_DIR] [--cwl-dir CWL_DIR] [--blastdb-dir BLASTDB_DIR]
#                      [--term-flex TERM_FLEX] [-q] [-d] [-rt RT] [--cwl_args CWL_ARGS] [--prok | --euk]
#                      [input ...]
#
#Run Adaptor and Vector screens. By default, screen against adaptors only.
#
#positional arguments:
#  input                 Input FASTA file(s).
#
#options:
#  -h, --help            show this help message and exit
#  -o path, --output path
#                        Output directory to be created, which may include a full path
#  --bin-dir BIN_DIR     bin directory
#  --cwl-dir CWL_DIR     CWL directory
#  --blastdb-dir BLASTDB_DIR
#                        blastdb directory
#  --term-flex TERM_FLEX
#                        terminal flexibility
#  -q, --quiet           Quiet mode, for scripts
#  -d, --debug           Debug mode
#  -rt RT, --rt RT       GCP testing, true for real testing, false for force-failure, absent for normal
#                        behavior
#  --cwl_args CWL_ARGS   Extra CWL args.
#  --prok                Prokaryotes: use adaptors_for_proks BLAST db. Also screen against contam_in_prok
#                        (not implemented).
#  --euk                 Eukaryotes: use adaptors_for_euks BLAST db (default). Also screen against
#                        gcontam1 (not implemented).