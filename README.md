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

If basecalling gets interrupted, you can resume basecalling by providing the `--resume-from` flag.
```
dorado basecaller --trim all sup,5mCG_5hmCG $POD5_DIR --resume-from $BAM_DIR/simplex_5mCG_5hmCG.bam > $BAM_DIR/simplex_5mCG_5hmCG_resumed.bam
mv $BAM_DIR/simplex_5mCG_5hmCG_resumed.bam $BAM_DIR/simplex_5mCG_5hmCG.bam
```

### iii. convert to fastq
The resulting bam file needs to be converted to fastq using [samtools and bgzip from htslib](https://www.htslib.org/download/). Do not remove the bam file! This file stores the base modification scores which will be used later.

To assess the quality of the raw reads you can use [Nanoplot](https://github.com/wdecoster/NanoPlot) to visualise the read quality and read lentghs.  
```
# Convert bam to fastq - note that the base modifications are not written to the fastq file.
samtools fastq simplex_5mCG_5hmCG.bam > simplex.fastq

# Bgzip fastq file
bgzip simplex.fastq

# Run QC using Nanoplot
NanoPlot -t 10 --fastq bgzip simplex.fastq.gz -o PATH/TO/NAOPLOT/OUTPUT/DIRECTORY   
```

## 2. Read filtering
To remove low-quality reads and filter sequencing artifacts we'll use [chopper](https://github.com/wdecoster/chopper) and [clean](https://github.com/rki-mf1/clean). Again, we can run Nanoplot to compare how filtering changed effected our read quality and length distribution.

Remove reads that have a quality score below a phred-score below 8 `-q 8` and shorter than 100 basepairs `-l 100`. While dorado should remove all adapters and primers, I find it a good practice to remove the first 50 and last 20 basepairs to remove any unwanted sequences missed by dorado `--headcrop 50 --tailcrop 20`. Chopper loads the entire fastq file into memory, so make sure you request enough RAM for this job (below I've requested 20 threads and 200Gb of memory).


To run the de

```
# Run chopper
bgzip -cd simplex.fastq.gz | chopper -q 8  -t 20 -l 100 --headcrop 50 --tailcrop 20 | bgzip -c > simplex_porechop_q8_hc50_tc20_l100.fastq.gz

# Run clean
nextflow run rki-mf1/clean -r v1.1.2                  \
 --input_type nano                                    \
 --input simplex_porechop_q8_hc50_tc20_l100.fastq.gz' \
 --control dcs                                        \
 --output $CLEAN_DIR                                  \
 --cores 20                                           \
 --memory '200.GB'                                    \
 -profile conda

# Move clean reads and rename file
mv $CLEAN_DIR/clean/*.fastq.gz simplex_filtered.fastq.gz

# Run Nanoplot
NanoPlot -t 10 --fastq bgzip simplex_filtered.fastq.gz -o PATH/TO/NANOPLOT/OUTPUT/DIRECTORY
```

The final output contains high-quality reads that can be used for assembling a genome. Depending on the assembler, additional filter steps might be required. 

## 3. Haplotype aware error correction (dorado correct)
Dorado provides an error-correction module which is an integration of the [HERRO](https://github.com/lbcb-sci/herro) algorithm. This involves an all-vs-all alignment and subsequent haplotype-aware error correction, producing high accuracy "Hifi-like" reads (similar to those produced by Pac-Bio). We'll follow the protocol presented in the [HERRO paper](https://www.biorxiv.org/content/10.1101/2024.05.18.594796v2) for performing the additional filtering steps required to perform the error correction. 

We'll remove reads shorter than 10kbp (`-l 10000`) and a quality-score below 10 (`-q 10`). We'll also use this opportunity to filter Ultra-Long (UL) reads (> 50 kbp) which will be used for scaffolding during the assembly.
```
# Run chopper for error-correction input
bgzip -cd simplex_filtered.fastq.gz' | chopper -t 10 -l 10000 -q 10 | bgzip -c > simplex_filtered_10kbp_q10.fastq.gz

# Run chopper to extract ULreads
bgzip -cd simplex_filtered_10kbp_q10.fastq.gz | chopper -t 10 -l 50000 | bgzip -c > simplex_filtered_50kbp_q10.fastq.gz

#run Nanoplot for ULreads. you will need to infer error rate (required parameter for Hifiasm)
NanoPlot -t 10 --fastq bgzip simplex_filtered_50kbp_q10.fastq.gz -o PATH/TO/NANOPLOT/OUTPUT/DIRECTORY
```
Then use the resulting file to run dorado correct.
```
# Download model
dorado download --model herro-v1

# Perform correction
dorado correct -t 200 -m herro-v1 simplex_filtered_10kbp_q10.fastq.gz > hifi_reads.fasta
```
Note that the output is written in FASTA (not FASTQ). This implies that there are no quality scores associated and basecalls als assumed to be accurate.

The resulting reads are split into 30kbp fragments to mimic Pac-Bio reads. The assembler Hifiasm was build specifically to handle Pac-Bio data which are approximately 30kbp max in length. Submitting longer reads were observed to reduce the quality of the assembly, hence we'll be doing the same here using [seqkit](https://github.com/shenwei356/seqkit). Any split sequences shorter than 10kbp were subsequently removed.

```
# Split reads in 30kbp sections
seqkit sliding -j 10 -g -s 30000 -W 30000 hifi_reads.fasta | seqkit seq -m 10000 -w 0 > hifi_reads_split.fasta

# Ran Nanoplot
NanoPlot -t 10 --fasta $FASTA'_herro_hifi_split.fasta' -o $NANO/herro_hifi_split
```
You will notice there are strange peaks in the read length distribution. This is an artifact from the binning process of the error-correction algorithm and can be safely ignored ([see](https://github.com/lbcb-sci/herro/issues/67)).

## 4. Genome assembly ([Hifiasm](https://github.com/chhylp123/hifiasm))
We now have all the required input files to perform an assembly using Hifiasm. We'll employ two different strategies to generating a de novo assembly. 

1. Hybrid approach using split hifi reads + UL simplex reads
2. simplex only approach using filtered reads from dorado

Which approach works best for you will depend on the quality of your data and the complexity of the genome. It's worth performing an assembly using both methods and evaluating which approach resulted in the best assembly based on the number scaffolds and completeness (see Assembly statistics section).

### i. Hifiasm - hifi + ULreads
This approach uses the error-corrected hifi reads and ULreads as input. To obtain the error rate of the ULreads you can use average phredscore provided by Nanoplot and convert this to error probability ([online tool](https://jamiemcgowan.ie/bioinf/phred.html)). An error rate of 0.02 is roughly equal to a phredscore of 17.

The output directory will contain multiple assembly files. We'll focus on asm.bp.p_ctg.gfa (i.e. assembly of primary contigs). For a detailed explanation of the different output files see [here](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html).

The graph based output (gfa file) is converted to fasta, indexed using [samtools](https://www.htslib.org/download/), and scaffolds are reordered and renamed based on their length using [seqkit](https://github.com/shenwei356/seqkit).
```
# Output directory
HIFI_DIR=/PATH/HIFIASM/OUTPUT

# Run Hifiasm
hifiasm -o $HIFI_DIR/assembly.asm -t 32 --ul simplex_filtered_50kbp_q10.fastq.gz --ul-rate 0.02 hifi_reads_split.fasta

#convert gfa to fasta
cd $HIFI_DIR
awk '/^S/{print ">"$2;print $3}' assembly.asm.bp.p_ctg.gfa > assembly.asm.bp.p_ctg.fa

#reorder and rename contigs
seqkit sort --by-length --reverse assembly.asm.bp.p_ctg.fa | seqkit replace --pattern '.+' --replacement 'Contig{nr}' > assembly.asm.bp.p_ctg.reordered.fasta

# Index fasta file
samtools faidx assembly.asm.bp.p_ctg.reordered.fasta

# Remove intermediate file
rm $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.fa
```
### ii. Hifiasm - ONT-simplex 
From version v0.21.0 onwards, hifiasm incorporated a new function which allows you to submit ONT simplex reads as the primary input for the assembly.
```
HIFI_DIR=/PATH/HIFIASM/OUTPUT

hifiasm --ont -o $HIFI_DIR/assembly.asm -t32 simplex_filtered.fastq.gz

#convert gfa to fasta
cd $HIFI_DIR
awk '/^S/{print ">"$2;print $3}' assembly.asm.bp.p_ctg.gfa > assembly.asm.bp.p_ctg.fa

#reorder and rename contigs
seqkit sort --by-length --reverse assembly.asm.bp.p_ctg.fa | seqkit replace --pattern '.+' --replacement 'Contig{nr}' > assembly.asm.bp.p_ctg.reordered.fasta

# Index fasta file
samtools faidx assembly.asm.bp.p_ctg.reordered.fasta

# Remove intermediate file
rm $HIFI_DIR/$ASSEMBLY.asm.bp.p_ctg.fa
```
Again, you will need to evaluate yourself which approach resulted in the best assembly for your species. If you want to see how your assembly looks right after your assembly (which I would recommend), you can run gfastats, quast and BUSCO to assess the quality (7. assembly evaluation) 

## 5. Purge haplotigs ([purge_dups](https://github.com/dfguan/purge_dups))

```
add code
```

## 6. Contamination screening ([FCS](https://github.com/ncbi/fcs))

It's possible that DNA from other organisms (e.g. E.coli) or adapter sequences were not filtered out previous steps and represented in your assembly. You can use Foreign Contaminations Screen (FCS) to test for this, and sequently remove these sequences if present. FCS consists of of separate module which test for each type of contamination separately.

### i. [FCS-adapter](https://github.com/ncbi/fcs/wiki/FCS-adaptor-quickstart)

```
# Containers and scripts
FCS_ADAPTOR_SIF=/PATH/TO/FCS-ADAPTER/SINGULARITY/SIF
FCS_GX_SIF=/PATH/TO/FCS-GX/SINGULARITY/SIF
PY=/PATH/FCS/PYTHON/SCRIPT

# PATHS
FASTA=/PATH/TO/PURGED/ASSEMBLY
FCS_DIR=/PATH/TO/FCS/OUTPUT/DIR

# Create output directory
mkdir -p $FCS_DIR/fcs_adaptor

# Run fcs_adaptor
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
        --output $FCS_DIR/fcs_adaptor/assembly.fcs_adaptor.clean.fasta             \
        --contam-fasta-out $FCS_DIR/fcs_adaptor/assembly.fcs_adaptor.contam.fasta        
```

### ii. [FCS-GX](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart)

Running FCS-GX requires pre-building a large database which needs to be loaded into memory. 
first, Build the database. make sure you have enough free space (200Gb)
```
# FCS_gx database
LOCAL_DB="PATH/TO/GCS_GX/DATABASE"
PY=/PATH/FCS/PYTHON/SCRIPT

# Build database
SOURCE_DB_MANIFEST="https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"
python3 $PY db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/gxdb"

# Check if files are present 
#ls $LOCAL_DB/gxdb
```

Now you can run FCS-gx to perform contamination screening.

```
# Containers and scripts
FCS_GX_SIF=/PATH/TO/FCS-GX/SINGULARITY/SIF
PY=/PATH/FCS/PYTHON/SCRIPT

# PATHS
LOCAL_DB="PATH/TO/GCS_GX/DATABASE"
FASTA=/PATH/TO/FCS_ADAPTER/OUTPUT/FASTA
FCS_DIR=/PATH/TO/FCS/OUTPUT/DIR

# Create output directory
mkdir -p $FCS_DIR/fcs_gx

# Export relevant settings
export FCS_DEFAULT_IMAGE=$FCS_GX_SIF
export NCBI_FCS_REPORT_ANALYTICS=1

#taxomomic ID (search for on [ncbi](https://www.ncbi.nlm.nih.gov/taxonomy))
tax_id=334912 #example for Polyprion oxygeneios, hapuka 

#screen genome
python3 $PY screen genome --fasta $FASTA --out-dir $FCS_DIR/fcs_gx --gx-db $LOCAL_DB/gxdb --tax-id $tax_id 

#clean fasta
cat $READ_DIR/$FASTA | python3 $PY clean genome                          \
        --action-report $FCS_DIR/fcs_gx/*fcs_gx_report.txt               \
        --output $FCS_DIR/fcs_gx/assembly.fcs_gx.clean.fasta             \
        --contam-fasta-out $FCS_DIR/fcs_gx/assembly.fcs_gx.contam.fasta      

# Copy clean assembly to main FCS dir
cp $FCS_DIR/fcs_gx/assembly.fcs_gx.clean.fasta $FCS_DIR/assembly_final_assembly.fa

# Reorder and rename assembly contigs
seqkit sort --by-length --reverse $FCS_DIR/assembly_final_assembly.fa | seqkit replace --pattern '.+' --replacement 'Contig{nr}' > $FCS_DIR/assembly_final_assembly.fasta

# Index assembly
samtools faidx $FCS_DIR/assembly_final_assembly.fasta

# Remove intermediate file
rm $FCS_DIR/assembly_final_assembly.fa
```

## 7. Assembly evaluation
There are many tools out there that will allow you to evaluate the quality of your assembly. Below I outline a number of commonly used tools that are quick and easy to run. 

### i. Assembly metrics ([gfastats](https://github.com/vgl-hub/gfastats) & [quast](https://github.com/ablab/quast))
Gfastats and quast are great tools to extract the most important summary metrics about your assembly. The number of scaffolds, average scaffold length, N50-scaffold size ...
They report very similar summary metrics, see which one you find most useful. gfastats `--seq-report` provides metrics for each individuals scaffold. 

```
ASSEMBLY=/PATH/TO/PURGED/ASSEMBLY
STAT_DIR=/PATH/TO/STATS/DIRECTORY

mkdir -p $STAT_DIR

# Gfastats
gfastats -f $ASSEMBLY              > $STAT_DIR/$assembly.summary.tsv
gfastats -f $ASSEMBLY --seq-report > $STAT_DIR/$assembly.seq_summary.tsv

# Quast
quast --output-dir $STAT_DIR/quast --threads 2 $ASSEMBLY
```

### ii. Assembly completeness ([BUSCO](https://busco.ezlab.org/busco_userguide.html))
BUSCO is an incredibly useful tool that evaluates the completeness of your assembly by trying to identify genes present in a reference database. Which database you need to use depends on your organism. Here, I've used the `actinopterygii_odb10` database for fish.

```
# PATHS
ASSEMBLY=/PATH/TO/PURGED/ASSEMBLY
STAT_DIR=/PATH/TO/STATS/DIRECTORY
SIF=/PATH/TO/SINGULARITY/SIF

# Create output directory
mkdir $STAT_DIR/busco

# Run BUSCO
singularity run $SIF busco -i $ASSEMBLY               \
                           --out_path $STAT_DIR/busco \
       			           -o busco_assembly          \
			               -l actinopterygii_odb10	  \
				           -m genome         		  \
				           --cpu 10
```

### iii. Telomere identification ([tidk](https://github.com/tolkit/telomeric-identifier) & [quarTeT](https://github.com/aaranyue/quarTeT))
Because we've used Oxford Nanopore Technology and assume we've generated some ultra-long reads we will want to check if we have assembled entire chromosome. To do this we will try to identify the presence of a telomeric sequence at the ends of all contigs (possibly full chromosomes).

```
add code
```

## 8. Extract mitochondrial genome ([MitoHifi](https://github.com/marcelauliano/MitoHiFi))

```
add code
```

## 9. Repeat masking ([TETools](https://github.com/Dfam-consortium/TETools/))
In onder to perform functional annotation we first need to identify and mask repetetaive regions within the genome. This is a very time consuming process and take the most time to complete by far. I would highly recommand using the [dfam TE Tools container](https://github.com/Dfam-consortium/TETools/). I will describe how to run this container using singularity. 

Repeat masking involved the following steps:
1. Download + build container and dfam database
2. Extract repeart sequences from database
3. Identify repetative sequences within the genome
4. Mask Repeats

Carefully follow the instruntions on the dfam [TETools github page](https://github.com/Dfam-consortium/TETools/) on how to configure the TETools container by adding the nessecary dfam and RepBase libraries to the container. It's a somewhat complecated process which requires accessing the contents within the container and linking the downloaded libraries outside of the container. 

The dfam library is very large and split into multiple partitions so you only have download those sections relevent for your organism. Check the FamDB [page](https://github.com/Dfam-consortium/FamDB) for an explenation.
NOTE: below I extract repeats for `Actinopterygii`. You will have to replace this with the families that is appropriote for your species.

### i. Identifying repeats
```
#PATHS
ASSEMBLY=/PATH/TO/PURGED/ASSEMBLY
TETOOLS_LIB=/PATH/TO/PREBUILD/LIB
DFAM=PATH/TO/DFAM/TETOOLS/SIF


# Output
REPEAT_DIR=/PATH/TO/REPEATMASKING/OUTPUT/DIR
LIB_DIR=$REPEAT_DIR/libraries
RM_DB=$REPEAT_DIR/repeatmodeler_db

# Create and move to output directory
mkdir -p $LIB_DIR $RM_DB
cd $REPEAT_DIR

################################### build repeat library from known repeats ############################
# Extract repeats
singularity exec --bind $TETOOLS_LIB$:/opt/RepeatMasker/Libraries          \
              $DFAM famdb.py -i /opt/RepeatMasker/Libraries/famdb families \
              -f fasta_name -ad --include-class-in-name                    \ 
              --add-reverse-complement                                     \
              Actinopterygii > $LIB_DIR/Actinopterygii_library.fa

################################### build repeat library from draft genome #############################
#build de novo repeat library 
singularity exec $DFAM BuildDatabase -name $RM_DB/assembly_db $ASSEMBLY

#Identify repeats including LTRS (long tandem repeats)
singularity exec $DFAM RepeatModeler -threads 32                    \
                                     -LTRStruct                     \
                                     -database $RM_DB/assembly_db
mv RM*/ RepeatModelerOutput/

#move output to new folder and removing empty lines
awk 'NF' RepeatModelerOutput/consensi.fa > $LIB_DIR/cleaned_consensi.fa

########################  merge the two libraries if known and identified repeats  ######################
cat $LIB_DIR/cleaned_consensi.fa $LIB_DIR/Actinopterygii_library.fa > $LIB_DIR/combined_library.fa
```
Check whether both generated libraries (and thus the merged library as well) contain sequences. If one step fails and produces an empty file the rest will run without giving an error. 

### ii. Masking repeats
Make sure you request enough memory and time for RepeatMasker to run. It is not uncommon for RepeatMaske to run over 10 days. If needed, request if your runtime can be extended. You can check in the logs how far along your run is. This will give you idea of how far along your analysis is. 

```
ASSEMBLY=/PATH/TO/PURGED/ASSEMBLY

# Output paths
REPEAT_DIR=/PATH/TO/REPEATMASKING/OUTPUT/DIR
LIB_DIR=$REPEAT_DIR/libraries
DFAM=PATH/TO/DFAM/TETOOLS/SIF

# Run repeatmasker
singularity exec $DFAM RepeatMasker -pa 8                                   \
                                    -dir $REPEAT_DIR                        \
                                    -xsmall                                 \
                                    -lib $LIB_DIR/combined_library.fa $ASSEMBLY
```
The output contains a .tbl file. This txt file contains a summary of all identified sequences and what percentage of your genome has been masked. 

## 10. Gene Annotation
Gene Annotation involves the identification if protein coding sequences in the (unmasked regions of the) genome. We will perform gene annotation within protein sequences (RNA-seq data). If RNA-seq is available for your species, it highly recommended you include this in your analysis.

### i. Galba

### ii. ACAT
Because we only predicted 


## 11. Functional annotation

### i. EggNog-mapper

### ii. Interproscan








