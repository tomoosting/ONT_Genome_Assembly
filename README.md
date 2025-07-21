# De novo genome assembly using Oxford Nanopore 
!!REPO UNDER CONSTRUCTION!!

This repository will soon contain all documentation for performing a de novo genome assembly using Oxford Nanopore data.

This Repository contains a step by step workflow for performing a de novo assembly using only Oxford Nanopre Technology (ONT) data.
Our goal is to provide people with an easy to adopt protocol for generating a reference genome using a single sequencing technology. 

**IMPORTANT**: 
* The workflow assumes data has been generated using R10.4.1 flowcells with kit14 chemistry (or higher)
* You have some basic coding skills for setting up and executing bash scripts
* You have access or able to install required packages (singularity/apptainer and conda are your friends here)
* you access to HPC that has a GPU partition (dorado runs on GPU)
* I would recommand generating ~50X coverage of raw sequencing data.

**This workflow contains the following steps**
1. Basecalling ([dorado basecaller](https://github.com/nanoporetech/dorado)) 
2. Merge ubam files ([samtools](https://github.com/samtools))
3. Trimming and filtering of simplex reads (porechop, clean, FCX, chopper)
4. Quality control simplex reads (Nanoplot)
5. Haplotype-aware error correction ([dorado correct](https://github.com/nanoporetech/dorado))
6. Genome assembly (Hifiasm)
7. Mitochondiral genome assembly
8. Error correction [Medaka](https://timkahlke.github.io/LongRead_tutorials/ECR_ME.html)
9.  Genome assesment ([BUSCO](https://github.com/WenchaoLin/BUSCO-Mod))
10. Repeat masking
11. Annotation (Eggnog-mapper[https://github.com/eggnogdb/eggnog-mapper])
12. Functional annotation

## 1. basecalling with [dorado basecaller](https://github.com/nanoporetech/dorado)
### i. Convert fast5 (optional)
If your raw data is stored in fast5 files you will have to convert them to pod5. To do so, you can use [pod5](https://pod5-file-format.readthedocs.io/en/latest/).

```
#PATHS
FAST5_DIR=/PATH/TO/FAST5_DIR
POD5_DIR=/PATH/TO/POD5_DIR

#run pod5
pod5 convert fast5 $FAST5_DIR/*.fast5 --output $POD5_DIR/ --one-to-one $FAST5_DIR/
```
If you already have your raw data in pod5 you can skip this step.

### ii. Dorado basecalling
In order to perform basecalling using dorado you will need to know a number of things about how the data was generated.
1. Library preperation protocol was used (e.g. kit14 SQK-LSK114)
2. Flow cell was used (e.g. R10.4.1 FLO-PRO114M)
3. Machine and speed the flow cell was run on (e.g. P2 400 bps)
This information will let you pick the right model --> dna_r10.4.1_e8.2_400bps_sup@v4.1.0

You can downlaod the latest dorado version from [here](https://github.com/nanoporetech/dorado/releases). By aware, from v1.0.0 you can no longer run 4kHz models. If you used older chemistry, download an older version (v0.X.X).


```
POD5_DIR=PATH/TO/POD5_DIR
BAM_DIR=PATH/TO/BAM_DIR

dorado basecaller --trim all sup,5mCG_5hmCG $POD5_DIR > $BAM_DIR/simplex_5mCG_5hmCG.bam
```
By default dorado trims all adapters and barcodes but I prefer to explicitely use the parameter `--trim all`. Dorado had various basecalling models that balance speed and acuracy. The sup (superior) model is the slowest but results in the highest quality. Oxford Nanopore also allows for scoring of basemodification like 5mCG and 5hmCG sites. I recommand scoring all possible base modification available for the chemistry you used. For more information on basecalling models and how to select them see [here](https://github.com/nanoporetech/dorado?tab=readme-ov-file#available-basecalling-models).

If you have issues with running time (like me), you can split up job into multiple smaller ones. See 1_1_submit_array_dorado.sh and 1_dorado.sh in the scritps directory. You will have to merge all resulsting bam files into a single file using `bamtools merge`.

### iii. Read filtering
Now that we have our raw basecalling data we can filter low 



Again, for the dorado basecalling I've written an array script that allows you to basecall all your pod5 files in parallel, drastically cutting back run time. In some cases you might not be able to request enough run time to basecall all your pod5 files in series. Dorado requests a directory is input and performs casecalling on all pod5 files contained in the that directory. The following array script copies a single pod5 file to a temporary directory and performs basecalling.  
```
#!/bin/bash
#SBATCH -a 1-10
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=39G
#SBATCH --time=0-0:30:00
#SBATCH --job-name=dorado

# Code to set up GPU partition **!!!THIS WILL BE SPECIFIC TO YOUR HPC!!!**
module load GCC/10.2.0
module load CUDA/11.1.1
module load OpenMPI/4.0.5

# Program
dorado=~/bin/dorado-0.3.1-linux-x64/bin/dorado

# Download model
$dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.1.0

# Variables
i=${SLURM_ARRAY_TASK_ID}

# Paths
POD5_DIR=PATH/TO/POD5/DIR
BAM_DIR=PATH/TO/BAM/OUTPUT/DIR
mkdir -p $BAM_DIR

# Extract filename from the ith pod5 file in directory
NAME=$( ls $POD5_DIR/*.pod5 | head -n $i | tail -n 1 | sed -e 's/\.pod5$//' | xargs -n 1 basename )

# Create tmp dir and copy pod5 file
TMP_DIR=$BAM_DIR/$NAME/
mkdir -p $TMP_DIR
cp $POD5_DIR/$NAME.pod5 $TMP_DIR

# Run dorado
$dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.1.0 $TMP_DIR > $BAM_DIR/$NAME.bam 

# Remove tmp dir
rm -r $TMP_DIR
```
#### merge bam files and convert bam to fastq
The next step is merging all bam files and convert with samtools, then compress using bgzip.
```
samtools merge $BAM_OUT.bam $BAM_OUT/*.bam
samtools fastq $BAM_OUT.bam > $FASTQ_OUT.fastq
bgzip $FASTQ_OUT.fastq
```
#### create summary file

```
dorado summary $BAM_EXT.bam > $BAM_EXT.sequencing_summary.txt
```

#### remove adapters with [porechop](https://github.com/rrwick/Porechop)
Like other sequencing platforms, ONT reads can have adapter sequences attached. These we'll rempove with porechop. This program is no longers upported but it still seems to be the go to program.
I've had some issues with installing the program due to compatibility issues on my HPC, but I've found a singulatiry container [here](https://forgemia.inra.fr/gafl/singularity/porechop). If you don't have experience working with singulatiry have a look at their [page](https://docs.sylabs.io/guides/3.1/user-guide/index.html). First download the image.
```
singularity pull porechop_v0.2.3.sif oras://registry.forgemia.inra.fr/gafl/singularity/porechop/porechop:latest
```
Then we can run the program as follows:
```
singularity exec porechop_v0.2.3.sif porechop
```

## 2. Quality control

PycQC python3 application which generates summaries based on output from `dorado summary calls.bam` and the `sequencing_summary.txt` file generated by guppy.
The easiest wat to run PycoQC is in a virtual environment. See this [link](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/) on how to install and use virtual python environments.
in short how to install a vitrual environment (don't do this on the login node of your HPC, admins and other users won't like this)
```
# Install package
python3 -m pip install --user virtualenv 
# Create vitrual environment names pycoqc_env
python3 -m venv ~/python/pycoqc_env
# Activate environment
source ~/python/pycoqc_env/bin/activate
# You are not inside your virtunal environment where you can install pycoQC
pip install pycoQC
# Exit virtual environment
deactivate 
```
Now that pycoQC is installed we can run it using the following code
```
source ~/python/pycoqc_env/bin/activate
pycoQC -f sequncing_summary.txt \
       -o output.html           \
       -j output.json           \
       --min_pass_qual 7
```
## 4. Genome assembly
If your FASTQ data looks good we can start creating a draft assembly using [flye](https://github.com/fenderglass/Flye). 
First, let's install Flye on your system:
```
# Navigate to where you would want your installation, for example
cd ~/bin/
# Download Flye, go into the directory and compile 
git clone https://github.com/fenderglass/Flye
cd Flye
make
```
We can now run Flye by typing `~/bin/Flye/bin/Flye`.

If you have done multiple runs on the ONT to generate the data, copy all FASTQ.gz files you want to use for your assembly to a single directory (i.e. FLYE_DATA_DIR)
```
# Point towards your Flye installation
flye=~/bin/Flye/bin/flye
# Run Flye on your data
$flye --nano-hq FLYE_DATA_DIR/*.fastq.gz \
      --out-dir FLYE_OUTPUT              \
      --genome-size 0.7g                 \
      --threads 10                       \
      --read-error 0.03                  \
      --no-alt-contigs                   \
      --scaffold
```
Note that the output directory must not exist, flye will create or give an error if the directory is already present. You also need to provide a rough estimate of your genome size `--genome-size`, here I expect the genome to be around 700Mb or `0.7g`. For Q20 ONT data the [Flye manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) suggests using a `--read-error` of 0.03. I prefer to run Flye with `--no-alt-contigs` and `--scaffold` as this removes alternative alleles (haplotigs) and produces the most contiguous assmebly.

## 5. genome assesment 
To assess how completenes and contiguous the assembly is we can use [BUSCO](https://busco.ezlab.org/busco_userguide.html) and [assembly-stats](https://github.com/rjchallis/assembly-stats).

### BUSCO
BUSCO tries to find a set of genes to determine completeness of the genome.
We will use Singularity to to run BUSCO.
First download the Singularity image (.sif)
```
**!!!THIS WILL BE SPECIFIC TO YOUR CLUSTER!!!**
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.7.3

singularity pull docker://ezlabgva/busco:v5.2.2_cv1
```

```
**!!!THIS WILL BE SPECIFIC TO YOUR CLUSTER!!!**
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.7.3

singularity run busco:v5.2.2_cv1.sif busco -i $ASSEMBLY_DIR/assembly.fasta  \ 
                                           -o busco_assembly                \
                                           -l actinopterygii_odb10          \
                                           -m genome                        \
                                           --cpu 10
```


### Assembly-stats
Assembly-stats procuces a circle plot. To run this using the provides script you'll have to download the circle_plot folder in this repository and place it in the same directory as the script.
```
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=circle_plot

PROJECT=$1 #species name

#directory where your flye assembly.fasta is
ASSEMBLY_DIR=/PATH/TO/FLYE/ASSEMBLY  

#output directory where output is written to
OUT_DIR=/PATH/TO/BUSCO/DIR
mkdir -p $OUT_DIR

#perl script
asm2stats=./circle_plot/pl/asm2stats.pl

#copy template folder to OUT_DIR
cp -r ./circle_plot/* $OUT_DIR

#create json file
echo "var ${PROJECT} = " > $OUT_DIR/json/$PROJECT'.json'
perl $asm2stats $ASSEMBLY_DIR/assembly.fasta >> $OUT_DIR/json/$PROJECT'.json'
echo "localStorage.setItem('${PROJECT}',JSON.stringify(${PROJECT}))" >> $OUT_DIR/json/$PROJECT'.json'

#add json to html
sed -i s"%<!--add_jsons_here-->%  <!--add_jsons_here-->\n  <script type=\"text/javascript\" src=\"json/${PROJECT}.json\"></script>%"g $OUT_DIR/assembly-stats.html
```
Afterwards, copy the entire circle_plot folder this is in `OUT_DIR` to your PC and open up the assemble-stats.html file to view the results.

## 9. Repeat Masking
In order to annotate the genome we first need to mask repeat regions in the genome. To do this we'll use singularity image from [dfam](https://github.com/Dfam-consortium/TETools). Repeats will be masked based on transposabe elements (TE) identified in the assembly and known repeat sequences present in Dfam and RepBase libraries. While the singularity container comes with a Dfam library, there is a much more extentive one (PLEASE NOTE: the extensive dfam library is 700Gb to make sure you have enough space on your system.). There is also an old (but free) RepBase library you can download. We will first download these two liraries and then download and configure the singularity image.  
```
#create a directory to download the dfam and RepBase liraries
mkdir /path/to/Repeat_masker_libs
cd /path/to/Repeat_masker_libs

#download Dfam and check if download is corrupted
wget "https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz"
wget "https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz.md5"
md5sum check Dfam.h5.gz.md5
gzip -d Dfam.h5.gz > Dfam.h5

#download RepBase
#if wget doesn't work you can download the file manually from github and transfer to your system (this file is not big)
wget https://github.com/yjx1217/RMRB/blob/master/RepBaseRepeatMaskerEdition-20181026.tar.gz
tar -xzf RepBaseRepeatMaskerEdition-20181026.tar.gz
mv ./Libraries/* .
rm -r ./Libraries
```
This should produce 3 files: `Dfam.h5`, `RMRBSeqs.embl` and `README.RMRBSeqs`. 
Next we will download the singulatity imamge `dfam-tetools-latest.sif` and update the libraries.
```
#download and create singularity image
singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest

#enter the container in interactive mode and copy the lirary to outside the container
singularity shell /path/to/dfam_container/dfam-tetools-latest.sif
cp -r /opt/RepeatMasker/Libraries/ /path/to/Repeat_masker_libs
exit

#remove the old Dfam lib and move the new libraries to Liraries fodler
rm ./Libraries/Dfam.h5
mv * ./Libraries/

#update the singularity image
singularity exec dfam-tetools-latest.sif addRepBase.pl -libdir Libraries
```
Now that we have created the RepeatMasker singularity image, downloaded the Dfam and RepBase libraries, and updated the singularity image.

Now we have prepaired our singularity image so we can run both [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) and [RepeatMasker](https://www.repeatmasker.org/). This section includes the follwoing steps:
1. build a database from your assembly
2. Identify TEs
3. extract known repeat sequences 
4. merge output from steps 2 & 3 
5. mask repeats
6. convert to hardmasked
```
#export PATH to libraries created in previous section
export LIBDIR=/path/to/Repeat_masker_libs/Libraries

#create output directory for outout
PURGE_DIR=/path/to/purge_haplpotypes/
REPEAT_DIR=/path/to/repeat/output/
mkdir $REPEAT_DIR

#1. build database
mkdir $REPEAT_DIR/libraries
singularity run dfam-tetools-latest.sif BuildDatabase -name $REPEAT_DIR/libraries/repeat_modeler_db \
                                                      -engine ncbi                                  \
                                                      $PURGE_DIR/assembly.fasta

#2. Identify TEs
singularity run dfam-tetools-latest.sif RepeatModeler -threads 32                                 \
                                                      -LTRStruct                                  \
                                                      -database $REPEAT_DIR/libraries/$ASSEMBLY'_db'

#3. check if taxonomic group is available 
singularity exec dfam-tetools-latest.sif famdb.py -i /nfs/scratch/oostinto/RepeatMasker/Libraries/RepeatMaskerLib.h5 lineage -ad 'Actinopterygii'
#extract repeats
singularity exec dfam-tetools-latest.sif famdb.py -i /nfs/scratch/oostinto/RepeatMasker/Libraries/RepeatMaskerLib.h5 families -ad --add-reverse-complement Actinopterygii > $REPEAT_DIR/Actinopterygii_library.fa

#4. merge dfam and repeatmodeler databases
cat $REPEAT_DIR/libraries/repeat_modeler_db-families.fa $REPEAT_DIR/libraries/Actinopterygii_library.fa > $REPEAT_DIR/libraries/combined_library.fa

#5. run repeatmasker
singularity exec dfam-tetools-latest.sif RepeatMasker -pa 32                   \
                                                      -dir $REPEAT_DIR         \
                                                      -xsmall                  \
                                                      -lib $REPEAT_DIR/libraries/combined_library.fa $PURGE_DIR/assembly.fasta
mv $REPEAT_DIR/results/assembly.fasta.masked $REPEAT_DIR/assembly.softmasked.fasta

#6. convert softmasked to hardmasked
sed 's/(acgt)/N/g' $REPEAT_DIRassembly.softmasked.ordered.fasta > $REPEAT_DIR/results/assembly.hardmasked.ordered.fasta
```

## 10. Genome annotation

To annotate the genome we well use [braker3](https://github.com/Gaius-Augustus/BRAKER), a fully automated pipeline that can uitlise RNA-seq + genome assembly and genome assembly only. Because we don't have RNA-seq data we will use the genome assembly only pipeline. Have a look at the braker github to get a detailed description of the steps.

First, we'll download the nessesary files. And again, we'll use a singularity container.
```
#build container
singularity build braker3.sif docker://teambraker/braker3:latest

#Download [ortho-db](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/) tuitable for your species. Here I download a protein database for vertebrates. 
wget "https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Vertebrata.fa.gz"
gunzip Vertebrata.fa.gz > Vertebrata.fa
```
You will also need to obtain the config folder from Augustus and place it at a writable location. You can find the config fodler on the github page from [Augustus](https://github.com/Gaius-Augustus/Augustus). Download the config folder and put it on a an accessable location.

Braker runs on the soft-masked version of your assembly using the following code
```
#Augustus config
AUG_CONFIG=/nfs/scratch/oostinto/scripts/genome_assembly/7_Annotation/config
#ortho_db
ORTHO_DB=/nfs/scratch/oostinto/databases/ortho_db11/Vertebrata.fa
SOFTMASKED_FASTA=$READ_DIR/$ASSEMBLER/4_repeat_masking/...
BRAKER_DIR=$READ_DIR/$ASSEMBLER/6_annotation/

#create output directory
mkdir -p $BRAKER_DIR

#run braker3
singularity exec -B braker3.sif braker.pl --AUGUSTUS_CONFIG_PATH=$AUG_CONFIG   \
                                          --species=species_name               \
                                          --genome=$SOFTMASKED_FASTA           \
                                          --prot_seq=$ORTHO_DB                 \
                                          --workingdir=$BRAKER_DIR             \
                                          --threads=32                         \                                 --busco_lineage=actinopterygii_odb10 \
                                          --gff3
```
note that I use the same lineage for BUSCO as earlier

Finally, to identify genes I would recommand using [blast2go](https://www.blast2go.com/). 
You can run this on either your desktop or cluster. I personally run this locally on my laptop.
I use the swissprot database which is a curated database.
You will also need the output from braker i.e. braker.codingseq which is used is the input.
Be aware, it may run a few days on your laptop.







