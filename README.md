# De novo genome assembly using Oxford Nanopore 
## !!REPO UNDER CONSTRUCTION!!
This repository will soon contain all documentation for performing a de novo genome allembly using Oxford Nanopore data.

This Repository contains a step by step workflow for performing a de novo assembly using only Oxford Nanopre Technology (ONT) data.
Our goal is to provide people with an easy to adopt protocol for generating a reference genome.


IMPORTANT: Only use ONT data that has generated using kit12 and R10.4 flowcells or higher. This is ONTs Q20 chemistry which genrates high-quality (~20 phred score) long-read data. Older chemistries will not be compatible for this work flow.

This workflow contains the following steps
1. Basecalling ([dorado](https://github.com/nanoporetech/dorado) or [guppy](https://timkahlke.github.io/LongRead_tutorials/BS_G.html))
2. Quality control ([pycoQC](https://github.com/a-slide/pycoQC))
3. Genome assembly ([Flye](https://github.com/fenderglass/Flye))
4. Genome assesment ([BUSCO](https://github.com/WenchaoLin/BUSCO-Mod), [assembly-stats](https://assembly-stats.readme.io/docs))
5. More steps will be added soon. 

## What you do you need
* raw read data from the oxford nanopore
* settings under which the data was genrated
* rough estimate of genomes size 
* a High Performance Computing cluster with a GPU partition

## 1. basecalling
You can perform basecalling using either **dorado** or **guppy**. Dorado is the preferred basecaller but does not support all the output configurations from ONT platforms. It that case guppy still provides a good alternative to perform basecalling. This is really only an issue if you have used kit12 chemistry with R10.4 flow cells. Dorado does support basecalling for kit10 and 9.4.1 flow cells and kit14 and R10.4.1 flow cells or newer.



### basecaling with [dorado](https://github.com/nanoporetech/dorado)
#### convert with pod5
We'll use dorado to convert the raw output from your sequencing run to fastq data. Often you receive both fast5 and fastq files. In such cases basecalling is performed while output was genereted. However, it's best to redo the basecalling using a hyper acurate model.

First we need to convert the fast5 files for POD5 using [pod5](https://pod5-file-format.readthedocs.io/en/latest/).
To speed things up I've written code for an array.

```
FAST5_DIR=/PATH/TO/FAST5_DIR
POD5_DIR=/PATH/TO/POD5_DIR
```
The following piece of code retrives the name of ith fast5 file in your FAST5_DIR and creates a new file extension for the POD5 file. Then pod5 converts the ith FAST5 file in the array to POD5. This will create the same number of pod5 files as fast5 files.
```
N=${SLURM_ARRAY_TASK_ID}
FAST5_N=$( ls $FAST5_DIR/*.fast5 | head -n $N | tail -n 1 )
POD5_N=$( basename "${FAST5_N%.*}" )
POD5_N=$( basename $POD5_N ) 
echo "converting fast5 to pod for $FAST5_N"
pod5 convert fast5 $FAST5_N --output $POD5_DIR/$POD5_N.pod5
```
If you're not comfotable with arrays or would like a single pod5 per library/run you can run code like this (but will take MUCH longer!)
```
pod5 convert fast5 $FAST5_DIR/*.fast5 --output $POD5_DIR/$output.pod5
```
#### basecall with dorado
[Dorado](https://github.com/nanoporetech/dorado) runs on GPU!

It's imporant you know which 
1. Library preperation protocol was used (e.g. kit14 SQK-LSK114)
2. Flow cell was used (e.g. R10.4.1 FLO-PRO114M)
3. Machine and speed the flow cell was run on (e.g. P2 400 bps)
This information will let chosoe the right model --> dna_r10.4.1_e8.2_400bps_sup@v4.2.0
```
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.2.0
```

Dorado basecaller 



### basecalling with [guppy](https://timkahlke.github.io/LongRead_tutorials/BS_G.html)
basecalling with guppy is a bit more straightforward, FASTQ files are generated straight from FAST5 files.
Guppy has both a GPU and CPU version but I'd recommand using the GPU version unless that is available on your HPC.

The following command will create a FASTQ.gz file for each FAST5 file. `-i` indicates the input directory and `-r` tells guppy to search recursively. `-x` indicates that guppy should use all available GPU cores. `-c` indates which configuration should be used based on how the data has been generated. to get a list of all availble models you can run `guppy_basecaller --print_workflows`. `--min-score 7` filters out any reads with a phrd score below 7 (default).
```
guppy_basecaller -i FAST5_DIR               \
                 -r                         \
                 -s FASTQ_TMP               \
                 -x "cuda:0"                \
                 -c dna_r10.4_e8.1_sup.cfg  \
                 --min_score 7              \
                 --num_callers 32
```
After basecalling all FASTQ can be concatenated using `cat`, and zipped using `bgzip` from [hitlib](https://github.com/samtools/htslib)
```
cat FASTQ_TMP/*fastq > FASTQ_DIR/final.fastq
bgzip FASTQ_DIR/final.fastq
```

## 2. Quality control

PycQC generates summaries based on output from `dorado summary calls.bam` and the `sequencing_summary.txt` file generated by guppy.
The easiest wat to run PycoQC is in a virtual python environment.

