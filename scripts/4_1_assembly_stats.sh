#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=assembly_stats
#SBATCH -o /nfs/scratch/oostinto/stdout/assembly_stats.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/assembly_stats.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

### SINGULARITY
singularity=~/singularity/bin/singularity
export SINGULARITY_BIND="/nfs/scratch/oostinto,/nfs/home/oostinto"
export SINGULARITY_TMPDIR="/nfs/scratch/oostinto/tmp"

### PARAMETERS
ASSEMBLY=$1
ASM_NAME=$(basename "$ASSEMBLY" | sed 's/\.[^.]*$//')
STAT_DIR=$(dirname "$ASSEMBLY")/assembly_stats/$ASM_NAME
mkdir -p $STAT_DIR

### GFASTATS
source activate /nfs/scratch/oostinto/conda/gfastats
    gfastats -f $ASSEMBLY              > $STAT_DIR/$ASM_NAME.summary.tsv
    gfastats -f $ASSEMBLY --seq-report > $STAT_DIR/$ASM_NAME.seq_summary.tsv
conda deactivate

### QUAST
source activate /nfs/scratch/oostinto/conda/quast
    quast  --output-dir $STAT_DIR/quast --threads 2  $ASSEMBLY
conda deactivate

### BUSCO
#BUSCO DOWNLOADED LINEAGES
DOWNLOADS=/nfs/scratch/oostinto/scripts/genome_assembly/4_assembly_stats/busco_downloads
#RUN BUSCO
mkdir $STAT_DIR/busco
$singularity run $SCRATCH/singularity/busco_5.7.0.sif busco -i $ASSEMBLY                \
                                                            --out_path $STAT_DIR/busco  \
									   						-o busco_$ASM_NAME          \
															--offline 					\
															--download_path $DOWNLOADS	\
									   						-l actinopterygii_odb10		\
									   						-f -m genome				\
									   						--cpu 10

### CIRCLE PLOT
asm2stats=/nfs/scratch/oostinto/scripts/genome_assembly/4_assembly_stats/circle_plot/pl/asm2stats.pl
#copy template folder to OUTPUT
if [ ! -f $STAT_DIR/circle_plot/assembly-stats.html ]
then
	mkdir $STAT_DIR/circle_plot
    cp -r ./circle_plot/* $STAT_DIR/circle_plot
fi
#create json file
echo "var $ASM_NAME = " > $STAT_DIR/circle_plot/json/$ASM_NAME.json
perl $asm2stats $ASSEMBLY >> $STAT_DIR/circle_plot/json/$ASM_NAME.json
echo "localStorage.setItem('${ASM_NAME}',JSON.stringify(${ASM_NAME}))" >> $STAT_DIR/circle_plot/json/$ASM_NAME.json
#add json to html
sed -i s"%<!--add_jsons_here-->%  <!--add_jsons_here-->\n  <script type=\"text/javascript\" src=\"json/${ASM_NAME}.json\"></script>%"g $STAT_DIR/circle_plot/assembly-stats.html

### TIDK - TELOMERE IDENTIFICATION
source activate /nfs/scratch/oostinto/conda/tidk
    mkdir -p $STAT_DIR/tidk
    cd $STAT_DIR/tidk
    
    echo "running tidk find"
    tidk find --print
    tidk find --clade Perciformes --dir $STAT_DIR/tidk --output $ASM_NAME --log $ASSEMBLY
    
    echo "running tidk plot"
    tidk plot --tsv $ASM_NAME'_telomeric_repeat_windows.tsv' --output $ASM_NAME
conda deactivate

### QUARTET - TELOMERE IDENTIFICATION
source activate /nfs/scratch/oostinto/conda/quartet 
    mkdir -p $STAT_DIR/quartet
    cd $STAT_DIR/quartet
    
    echo "running quarTeT"
    quartet_py=~/bin/quarTeT/quartet.py
    python3 $quartet_py TeloExplorer -i $ASSEMBLY -c animal -p $ASM_NAME
    #options(bitmapType='cairo')
conda deactivate


