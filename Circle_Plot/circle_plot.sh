#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=circle_plot
#SBATCH -o ~/stdout/circle_plot.%j.out
#SBATCH -e ~/stdout/circle_plot.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@gmail.com

#set variables
NAME="species"
FASTA="PATH/TO/ASSEMBLY.fasta"
OUT="/PATH/TO/OUTPUT/DIR"
DIR="/PATH/TO/Assembly-Stats_template"

#copy template folder to output
cp -r $DIR $OUT

#create json file
echo "var ${$NAME}_flye = " > $OUT/json/$NAME'_flye.json'
perl $DIR/pl/asm2stats.pl $FASTA >> $OUT/json/$NAME'_flye.json'
echo "localStorage.setItem('${$NAME}_flye',JSON.stringify(${NAME}_flye))" >> $out/json/$NAME'_flye.json'

#add json to html
sed -i s"%<!--add_jsons_here-->%  <!--add_jsons_here-->\n  <script type=\"text/javascript\" src=\"json/${NAME}_flye.json\"></script>%"g $OUT/assembly-stats.html
