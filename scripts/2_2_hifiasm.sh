#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=parallel
#SBATCH --time=0-12:00
#SBATCH --job-name=Hifiasm
#SBATCH -o /nfs/scratch/oostinto/stdout/Hifiasm.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/Hifiasm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#set parameters
PROJECT=$1
VERSION=v14 # v3 .. v14
ASSEMBLY=$PROJECT'_'$VERSION
ASSEMBLER=hifiasm

#binaries
export PATH=/nfs/home/oostinto/bin/samtools-1.18:$PATH
export PATH=/nfs/home/oostinto/bin:$PATH

#paths
READ_DIR=$SCRATCH/projects/$PROJECT/raw_data/dorado/
HIFI_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/hifiasm/1_hifiasm
STATS_DIR=$SCRATCH/scripts/genome_assembly/4_assembly_stats

#read files
#ONT-FASTQ
#ONT_simplex=$READ_DIR/$PROJECT'_simplex_filtered.fastq.gz'
#ONT_1kbp=$READ_DIR/$PROJECT'_simplex_filtered_1kbp_q10.fastq.gz'
#ONT_UL_10kbp=$READ_DIR/$PROJECT'_simplex_filtered_10kbp_q10.fastq.gz'
ONT_UL_20kbp=$READ_DIR/$PROJECT'_simplex_filtered_20kbp_q10.fastq.gz'
ONT_UL_50kbp=$READ_DIR/$PROJECT'_simplex_filtered_50kbp_q10.fastq.gz'
#HIFI-FASTA
HIFI=$READ_DIR/$PROJECT'_herro_hifi.fasta'
HIFI_35x=$READ_DIR/$PROJECT'_herro_hifi_35x.fasta'
HIFI_SPLIT=$READ_DIR/$PROJECT'_herro_hifi_split.fasta'
HIFI_SPLIT_35x=$READ_DIR/$PROJECT'_herro_hifi_split_35x.fasta'

#rm -r $HIFI_DIR #watch out! 
mkdir -p $HIFI_DIR 

#--ul-rate 50kbp
# hapuka  = 0.020
# kahawai = 0.016
# moki    = 0.025
ul_rate=0.020

#run hifiasm - corrected reads and UltraLong ONT simplex reads (50kbp)
source activate /nfs/scratch/oostinto/conda/hifiasm_v0.25.0
    #subsampling split reads and 50kbp UL-reads (v3)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_50kbp --ul-rate $ul_rate $HIFI_SPLIT_35x
    #subsampling  split reads and 20kbp UL-reads (v4)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_20kbp --ul-rate $ul_rate $HIFI_SPLIT_35x
    #no subsampling  split reads and 50kbp UL-reads (v5)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_50kbp --ul-rate $ul_rate $HIFI_SPLIT
    #no subsampling  split reads and 20kbp UL-reads (v6)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_20kbp --ul-rate $ul_rate $HIFI_SPLIT
    #subsampling and 50kbp UL-reads (v7)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_50kbp --ul-rate $ul_rate $HIFI_35x
    #subsampling and 20kbp UL-reads (v8)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_20kbp --ul-rate $ul_rate $HIFI_35x
    #no subsampling and 50kbp UL-reads (v9)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_50kbp --ul-rate $ul_rate $HIFI
    #no subsampling and 20kbp UL-reads (v10)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 --ul $ONT_UL_20kbp --ul-rate $ul_rate $HIFI
    #hifi reads only, split reads + downsampled (v11)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 $HIFI_SPLIT_35x
    #hifi reads only, split reads (v12)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 $HIFI_SPLIT
    #hifi reads only + downsampled (v13)
    #hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 $HIFI_35x
    #hifi reads only (v14)
    hifiasm -o $HIFI_DIR/$ASSEMBLY.asm -t 32 $HIFI
conda deactivate

#convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.gfa > $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.fa

#reorder and rename contigs
seqkit sort --by-length --reverse $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.fa | seqkit replace --pattern '.+' --replacement 'Contig{nr}' > $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.reordered.fasta
samtools faidx $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.reordered.fasta
rm $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.fa

#get stats
sbatch /nfs/scratch/oostinto/scripts/genome_assembly/4_assembly_stats/assembly_stats.sh $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.reordered.fasta

##assembly stats
#source activate /nfs/scratch/oostinto/conda/gfastats
#    gfastats -f $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.reordered.fasta              > $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.reordered.summary.tsv
#    gfastats -f $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.reordered.fasta --seq-report > $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.reordered.seq_summary.tsv
#conda deactivate

#run circleplot and BUSCO
#cd $STATS_DIR
#sbatch circle_plot.sh $PROJECT $VERSION $ASSEMBLER 1_hifiasm
#sbatch busco_singularity.sh $PROJECT $VERSION $ASSEMBLER 1_hifiasm
#sbatch quast.sh $PROJECT $VERSION

#clean up
rm *.asm.re.* *.asm.ovlp.* *.asm.ec.* *.asm.uidx.* *.asm.ul.ovlp.* *_utg.* *.ec.bin *.ovlp.*