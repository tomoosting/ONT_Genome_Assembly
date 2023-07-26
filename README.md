# De novo genome assembly using Oxford Nanopore 
## !!REPO UNDER CONSTRUCTION!!
This repository will soon contain all documentation for performing a de novo genome allembly using Oxford Nanopore data.

This Repository contains a step by step workflow for performing a de novo assembly using only Oxford Nanopre Technology (ONT) data.
Our goal is to provide people with an easy to adopt protocol for generating a reference genome.


IMPORTANT: Only use ONT data that has generated using kit12 and R10.4 flowcells or higher. This is ONTs Q20 chemistry which genrates high-quality (~20 phred score) long-read data. Older chemistries will not be compatible for this work flow.

This workflow contains the following steps
1. Basecalling ([Dorado](https://github.com/nanoporetech/dorado))
2. Quality control ([pycoQC](https://github.com/a-slide/pycoQC))
3. Genome assembly ([Flye](https://github.com/fenderglass/Flye))
4. Genome assesment ([BUSCO](https://github.com/WenchaoLin/BUSCO-Mod), [assembly-stats](https://assembly-stats.readme.io/docs))
5. More steps will be added soon. 

## What you do you need
1. raw read data from the oxford nanopore
2. settings under which the data was genrated
3. rough estimate of genomes size
4. 

### 1. basecalling

We'll use dorado to convert the raw output from your sequencing run to fastq data. Often you receive both fast5 and fastq files. In such cases basecalling is performed while output was genereted. However, it's best to redo the basecalling using a hyper acurate model.

First we need to convert the fast5 files for POD5 using [pod5](https://pod5-file-format.readthedocs.io/en/latest/).
To speed things up I've written code for an array.

```
FAST5_DIR=
POD5_DIR=
```

```
FAST5_N=$( ls $FAST5_DIR/*.fast5 | head -n $N | tail -n 1 )
POD5_N=$( basename "${FAST5_N%.*}" )
POD5_N=$( basename $POD5_N ) 
echo "converting fast5 to pod for $FAST5_N"
pod5 convert fast5 $FAST5_N --output $POD5_DIR/$POD5_N.pod5
```
This will create the same number of pod5 files as fast5 files.
If you're not comfotable with arrays or would like only a single pod5 file you can run code like this (but can take MUCH longer!)
```
pod5 convert fast5 $FAST5_DIR/*.fast5 --output $POD5_DIR/$output.pod5
```