#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --time=10-0:00:00
#SBATCH --job-name=Mask_Repeats_Hifi
#SBATCH -o /nfs/scratch/oostinto/stdout/Mask_Repeats_Hifi.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/Mask_Repeats_Hifi.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###load modules
module load GCC/10.2.0 OpenMPI/4.0.5 Singularity/3.10.2

#singularity
#singularity=~/singularity/bin/singularity
#export PATH=~/singularity/bin/:$PATH

#bind scratch and home to container
export SINGULARITY_BIND="/nfs/scratch/oostinto,/nfs/home/oostinto"
export SINGULARITY_TMPDIR="/nfs/scratch/oostinto/tmp"

###container
dfam=$SCRATCH/singularity/dfam-tetools-latest.sif

###parameters
PROJECT=$1
VERSION=$2
ASSEMBLER=hifiasm
ASSEMBLY=$PROJECT'_'$VERSION

###PATHS
#draft genome
DRAFT=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/hifiasm/3_fcs/${ASSEMBLY}_final_assembly.fasta

#output paths
REPEAT_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/$ASSEMBLER/4_repeat_masking
LIB_DIR=$REPEAT_DIR/libraries
RM_DB=$REPEAT_DIR/repeatmodeler_db
mkdir -p $LIB_DIR $RM_DB
cd $REPEAT_DIR

################################### build repeat library from draft genome #############################
#build de novo repeat library 
echo "creating RM database"
singularity exec $dfam BuildDatabase -name $RM_DB/$ASSEMBLY'_db' $DRAFT
#Identify repeats including LTRS (long tandem repeats)
echo "running repeatmodeler"
singularity exec $dfam RepeatModeler -threads 32                    \
                                     -LTRStruct                     \
                                     -database $RM_DB/$ASSEMBLY'_db'
mv RM*/ RepeatModelerOutput/
#move output to new folder and removing empty lines
awk 'NF' RepeatModelerOutput/consensi.fa > $LIB_DIR/cleaned_consensi.fa

################################### build repeat library from known repeats ############################
#check if taxonomic group is available
echo "finding family"
singularity exec --bind /nfs/scratch/oostinto/databases/repeatmasker/Libraries:/opt/RepeatMasker/Libraries  \
                 $dfam famdb.py -i /opt/RepeatMasker/Libraries/famdb families                               \
                 -ad Actinopterygii
#extract repeats
echo "extracting from dfam database"
singularity exec --bind /nfs/scratch/oostinto/databases/repeatmasker/Libraries:/opt/RepeatMasker/Libraries  \
              $dfam famdb.py -i /opt/RepeatMasker/Libraries/famdb families                                  \
              -f fasta_name -ad --include-class-in-name --add-reverse-complement                            \
              Actinopterygii > $LIB_DIR/Actinopterygii_library.fa

#merge dfam and repeatmodeler databases
cat $LIB_DIR/cleaned_consensi.fa $LIB_DIR/Actinopterygii_library.fa > $LIB_DIR/combined_library.fa

############################################## run repeat masker #######################################
##run repeatmasker
#echo "running repeatmasker"
#singularity exec $dfam RepeatMasker -pa 8                                   \
#                                    -dir $REPEAT_DIR                        \
#                                    -xsmall                                 \
#                                    -lib $LIB_DIR/combined_library.fa $DRAFT
#rename "${ASSEMBLY}_final_assembly.fasta" "${ASSEMBLY}_final_assembly.softmasked.fasta" $REPEAT_DIR/*
#mv $REPEAT_DIR/${ASSEMBLY}_final_assembly.softmasked.fasta.masked $REPEAT_DIR/${ASSEMBLY}_final_assembly.softmasked.fasta

######################################### create hard masked version ###################################
#convert softmasked to hardmasked
#cp $REPEAT_DIR/${ASSEMBLY}_final_assembly.softmasked.fasta $REPEAT_DIR/${ASSEMBLY}_final_assembly.hardmasked.fasta
#sed -i 's/a/N/g' $REPEAT_DIR/${ASSEMBLY}_final_assembly.hardmasked.fasta
#sed -i 's/c/N/g' $REPEAT_DIR/${ASSEMBLY}_final_assembly.hardmasked.fasta
#sed -i 's/t/N/g' $REPEAT_DIR/${ASSEMBLY}_final_assembly.hardmasked.fasta
#sed -i 's/g/N/g' $REPEAT_DIR/${ASSEMBLY}_final_assembly.hardmasked.fasta
#sed -i 's/ConNiN/Contig/g' $REPEAT_DIR/${ASSEMBLY}_final_assembly.hardmasked.fasta
