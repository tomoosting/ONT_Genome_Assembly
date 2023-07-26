#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=circle_plot
#SBATCH -o /nfs/scratch/oostinto/stdout/circle_plot.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/circle_plot.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

module load singularity/3.7.3

PROJECT=$1
assembly=$SCRATCH/projects/$PROJECT/output/assembly/flye/assembly.fasta

out=$SCRATCH/projects/$PROJECT/output/assembly/flye/circle_plot

perl=$SCRATCH/scripts/genome_assembly/flye/04_assembly_stats/circle_plot/pl/asm2stats.pl
#perl=$SCRATCH/scripts/assembly/assembly_stats/circle_plot/pl/asm2stats.minmaxgc.pl

#copy template folder to output
cp -r $SCRATCH/scripts/genome_assembly/flye/04_assembly_stats/circle_plot/ $out

#create json file
echo "var ${PROJECT}_flye = " > $out/json/$PROJECT'_flye.json'
perl $perl $assembly >> $out/json/$PROJECT'_flye.json'
echo "localStorage.setItem('${PROJECT}_flye',JSON.stringify(${PROJECT}_flye))" >> $out/json/$PROJECT'_flye.json'

#add json to html
sed -i s"%<!--add_jsons_here-->%  <!--add_jsons_here-->\n  <script type=\"text/javascript\" src=\"json/${PROJECT}_flye.json\"></script>%"g $out/assembly-stats.html
