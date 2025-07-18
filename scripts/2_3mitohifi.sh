#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=MitoHiFi
#SBATCH -o /nfs/scratch/oostinto/stdout/MitoHiFi.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/MitoHiFi.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#set parameters
PROJECT=$1
VERSION=v1
ASSEMBLY=$PROJECT'_'$VERSION
ASSEMBLER=hifiasm

###load modules
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2

#export PATH to singularity and bind scratch and home to container
#export PATH=~/singularity/bin/:$PATH
#export SINGULARITY_BIND="/nfs/scratch/oostinto,/nfs/home/oostinto"
#export SINGULARITY_TMPDIR="/nfs/scratch/oostinto/tmp"
mitohifi_sif=/nfs/scratch/oostinto/singularity/mitohifi.sif

#binaries
export PATH=/nfs/home/oostinto/bin/samtools-1.18:$PATH
export PATH=/nfs/home/oostinto/bin:$PATH

#paths
READ_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY
HIFI_DIR=$READ_DIR/hifiasm/1_hifiasm
MITO_DIR=$READ_DIR/hifiasm/0_mito
mkdir -p $MITO_DIR 
#data sources
READS_FASTA=$READ_DIR/$ASSEMBLY'_herro_corrected_split.fasta'
REF_GFA=$HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.gfa
REF_FA=$HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.fa

#convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' $REF_GFA > $REF_FA

#find genbenk close match - for some reason this doesn't work for the singnularity container
source ~/python/findMitoReference/bin/activate 
    python ~/MitoHiFi/src/findMitoReference.py --species "Latridopsis ciliaris" --outfolder $MITO_DIR --min_length 15000
deactivate
# hapuka   - Polyprion oxygeneios
# kahawai  - Arripis trutta
# moki     - Latridopsis ciliaris 
# johndory - Zeus faber

### annotation with mito_finder
#run mitohifi on fasta reads
#mkdir $MITO_DIR/reads_mito_finder
#cd $MITO_DIR/reads_mito_finder
#singularity exec --bind /nfs:/nfs $mitohifi_sif mitohifi.py -r $READS_FASTA -f $MITO_DIR/NC_*.fasta -g $MITO_DIR/NC_*.gb -t 4 -o 2 -a animal

#run mitohifi on fasta assembly (assembly)
mkdir $MITO_DIR/assembly_mito_finder
cd $MITO_DIR/assembly_mito_finder
singularity exec --bind /nfs:/nfs $mitohifi_sif mitohifi.py -c $REF_FA     -f $MITO_DIR/NC_*.fasta -g $MITO_DIR/NC_*.gb -t 10 -o 2 -a animal

### annotation with mitos
#run mitohifi on fasta reads
#mkdir $MITO_DIR/reads_mitos
#cd $MITO_DIR/reads_mitos
#singularity exec --bind /nfs:/nfs $mitohifi_sif mitohifi.py -r $READS_FASTA -f $MITO_DIR/NC_*.fasta -g $MITO_DIR/NC_*.gb -t 4 -o 2 -a animal --mitos

#run mitohifi on fasta assembly (assembly)
#mkdir $MITO_DIR/assembly_mitos
#cd $MITO_DIR/assembly_mitos
#singularity exec --bind /nfs:/nfs $mitohifi_sif mitohifi.py -c $REF_FA     -f $MITO_DIR/NC_*.fasta -g $MITO_DIR/NC_*.gb -t 10 -o 2 -a animal --mitos

#clean up 
#rm $REF_FA

#usage: MitoHiFi [-h] (-r <reads>.fasta | -c <contigs>.fasta) -f
#                <relatedMito>.fasta -g <relatedMito>.gbk -t <THREADS> [-d]
#                [-a {animal,plant,fungi}] [-p <PERC>] [-m <BLOOM FILTER>]
#                [--max-read-len MAX_READ_LEN] [--mitos]
#                [--circular-size CIRCULAR_SIZE]
#                [--circular-offset CIRCULAR_OFFSET] [-winSize WINSIZE]
#                [-covMap COVMAP] [-v] [-o <GENETIC CODE>]
#
#required arguments:
#  -r <reads>.fasta      -r: Pacbio Hifi Reads from your species
#  -c <contigs>.fasta    -c: Assembled fasta contigs/scaffolds to be searched
#                        to find mitogenome
#  -f <relatedMito>.fasta
#                        -f: Close-related Mitogenome is fasta format
#  -g <relatedMito>.gbk  -k: Close-related species Mitogenome in genebank
#                        format
#  -t <THREADS>          -t: Number of threads for (i) hifiasm and (ii) the
#                        blast search
#
#optional arguments:
#  -d                    -d: debug mode to output additional info on log
#  -a {animal,plant,fungi}
#                        -a: Choose between animal (default) or plant
#  -p <PERC>             -p: Percentage of query in the blast match with close-
#                        related mito
#  -m <BLOOM FILTER>     -m: Number of bits for HiFiasm bloom filter [it maps
#                        to -f in HiFiasm] (default = 0)
#  --max-read-len MAX_READ_LEN
#                        Maximum lenght of read relative to related mito
#                        (default = 1.0x related mito length)
#  --mitos               Use MITOS2 for annotation (opposed to default
#                        MitoFinder
#  --circular-size CIRCULAR_SIZE
#                        Size to consider when checking for circularization
#  --circular-offset CIRCULAR_OFFSET
#                        Offset from start and finish to consider when looking
#                        for circularization
#  -winSize WINSIZE      Size of windows to calculate coverage over the
#                        final_mitogenom
#  -covMap COVMAP        Minimum mapping quality to filter reads when building
#                        final coverage plot
#  -v, --version         show program's version number and exit
#  -o <GENETIC CODE>     -o: Organism genetic code following NCBI table (for
#                        mitogenome annotation): 1. The Standard Code 2. The
#                        Vertebrate MitochondrialCode 3. The Yeast
#                        Mitochondrial Code 4. The Mold,Protozoan, and
#                        Coelenterate Mitochondrial Code and the
#                        Mycoplasma/Spiroplasma Code 5. The Invertebrate
#                        Mitochondrial Code 6. The Ciliate, Dasycladacean and
#                        Hexamita Nuclear Code 9. The Echinoderm and Flatworm
#                        Mitochondrial Code 10. The Euplotid Nuclear Code 11.
#                        The Bacterial, Archaeal and Plant Plastid Code 12. The
#                        Alternative Yeast Nuclear Code 13. The Ascidian
#                        Mitochondrial Code 14. The Alternative Flatworm
#                        Mitochondrial Code 16. Chlorophycean Mitochondrial
#                        Code 21. Trematode Mitochondrial Code 22. Scenedesmus
#                        obliquus Mitochondrial Code 23. Thraustochytrium
#                        Mitochondrial Code 24. Pterobranchia Mitochondrial
#                        Code 25. Candidate Division SR1 and Gracilibacteria
#                        Code
