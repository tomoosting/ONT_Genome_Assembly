#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --mem=512G
#SBATCH --partition=bigmem
#SBATCH --time=1-00:00
#SBATCH --job-name=FCS_gx
#SBATCH -o /nfs/scratch/oostinto/stdout/FSC_gx.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/FSC_gx.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###load modules
#module load GCC/10.2.0 OpenMPI/4.0.5 Singularity/3.10.2

export PATH=~/singularity/bin:$PATH
export SINGULARITY_BIND="/nfs/scratch/oostinto,/nfs/home/oostinto"
export SINGULARITY_TMPDIR="/nfs/scratch/oostinto/tmp"

### BINARIES
export PATH=/nfs/home/oostinto/bin/samtools-1.18:$PATH
export PATH=/nfs/home/oostinto/bin:$PATH

#PARAMETERS
PROJECT=$1
VERSION=$2
ASSEMBLY=$PROJECT'_'$VERSION
ASSEMBLER=hifiasm

#PATHS
REF=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/hifiasm/2_purge_dups/$ASSEMBLY.asm.bp.p_ctg.reordered_purged.fasta
ASM_NAME=$(basename "$REF" | sed 's/\.[^.]*$//')
FASTA=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/hifiasm/3_fcs/fcs_adaptor/$ASM_NAME.fcs_adaptor.clean.fasta
FCS_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/hifiasm/3_fcs
mkdir -p $FCS_DIR/fcs_gx

#FCS
FCS_GX_SIF=/nfs/scratch/oostinto/singularity/fcs-gx.sif
PY=/nfs/scratch/oostinto/scripts/genome_assembly/4_assembly_stats/fcs.py
export FCS_DEFAULT_IMAGE=$FCS_GX_SIF
export NCBI_FCS_REPORT_ANALYTICS=1

#FCS_gx database
LOCAL_DB="$SCRATCH/databases/fcs"
#create database - run these commands only ones to build database
#SOURCE_DB_MANIFEST="https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"
#python3 $PY db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/gxdb"
#check if files are present 
#ls $LOCAL_DB/gxdb

#taxomomic IDs
tax_id=334912 #hapuka
#tax_id=270544 #kahawai
#tax_id=97160  #moki

#screen genome
python3 $PY screen genome --fasta $FASTA --out-dir $FCS_DIR/fcs_gx --gx-db $LOCAL_DB/gxdb --tax-id $tax_id 

#clean fasta
cat $READ_DIR/$FASTA | python3 $PY clean genome                          \
        --action-report $FCS_DIR/fcs_gx/*fcs_gx_report.txt               \
        --output $FCS_DIR/fcs_gx/$ASM_NAME.fcs_gx.clean.fasta            \
        --contam-fasta-out $FCS_DIR/fcs_gx/$ASM_NAME.fcs_gx.contam.fasta      
cp $SCRATCH/stdout/FSC_gx.$SLURM_JOB_ID.err $FCS_DIR/fcs_gx/FSC_gx.$SLURM_JOB_ID.err

# Copy clean assembly to main FCS dir
cp $FCS_DIR/fcs_gx/$ASM_NAME.fcs_gx.clean.fasta $FCS_DIR/$ASSEMBLY'_final_assembly.fa'
seqkit sort --by-length --reverse $FCS_DIR/$ASSEMBLY'_final_assembly.fa' | seqkit replace --pattern '.+' --replacement 'Contig{nr}' > $FCS_DIR/$ASSEMBLY'_final_assembly.fasta'
samtools faidx $FCS_DIR/$ASSEMBLY'_final_assembly.fasta'
rm $FCS_DIR/$ASSEMBLY'_final_assembly.fa'

#run assembly stats
sbatch $SCRATCH/scripts/genome_assembly/4_assembly_stats/assembly_stats.sh $FCS_DIR/$ASSEMBLY'_final_assembly.fasta'


































#trial run withou the use of fcs.py - but py script should get the preference
#singularity exec --bind $GXDB_LOC:/app/db/gxdb/  \
#                 --bind $READ_DIR:/sample-volume/ \
#                 --bind $GX_DIR:/output-volume/  \
#                 $sif python3 /app/bin/run_gx    \
#                 --fasta /sample-volume/$FASTA   \
#                 --out-dir /output-volume        \
#                 --gx-db /app/db/gxdb            \
#                 --tax-id $tax_id
