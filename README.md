# De novo genome assembly using Oxford Nanopore 
## !!REPO UNDER CONSTRUCTION!!
This repository will soon contain all documentation for performing a de novo genome allembly using Oxford Nanopore data.

This Repository contains a step by step workflow for performing a de novo assembly using only Oxford Nanopre Technology (ONT) data.
Our goal is to provide people with an easy to adopt protocol for generating a reference genome.


IMPORTANT: Only use ONT data that has generated using kit12 and R10.4 flowcells or higher. This is ONTs Q20 chemistry which genrates high-quality (~20 phred score) long-read data. Older chemistries will not be compatible for this work flow.

This flow contians the following steps
1. Basecalling (Dorado)
2. Quality control (pycoQC)
3. Genome assembly (Flye)
4. Genome assesment (BUSCO, Circleplot)
5. More steps will be added soon. 

## What you do you need
1. raw read data from the oxford nanopore
2. settings under which the data was genrated
3. rough estimate of genomes size
