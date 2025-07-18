#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --partition=parallel
#SBATCH --time=5-0:00:00
#SBATCH --job-name=Braker3_Hifiasm
#SBATCH -o /nfs/scratch/oostinto/stdout/Braker3_Hifiasm.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/Braker3_Hifiasm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###load modules
module load GCC/10.2.0
module load OpenMPI/4.0.5
module load Singularity/3.10.2

###parameters
PROJECT=$1
VERSION=v3
ASSEMBLER=hifiasm
ASSEMBLY=$PROJECT'_'$VERSION

###build container (only need to do ones)
#singularity build braker3.sif docker://teambraker/braker3:latest
braker3=$SCRATCH/singularity/braker3.0.8.sif
busco=$SCRATCH/singularity/busco_5.7.0.sif
agat=$SCRATCH/singularity/agat_1.4.0--pl5321hdfd78af_0.sif

###Augustus config
#download the config folder from https://github.com/Gaius-Augustus/Augustus
AUG_CONFIG=/nfs/scratch/oostinto/scripts/genome_assembly/7_Annotation/config

###download OrthoDB for vertebrate species and unzip
#wget "https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Vertebrata.fa.gz"
#gunzip Vertebrata.fa.gz > Vertebrata.fa
ORTHO_DB=/nfs/scratch/oostinto/databases/ortho_db11/Vertebrata.fa

###PATHS
#draft genome
FASTA=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/$ASSEMBLER/2_repeat_masking/$ASSEMBLY'.asm.bp.p_ctg.reordered.softmasked.fasta'
#output directory
BRAKER_DIR=$SCRATCH/projects/$PROJECT/output/assemblies/$ASSEMBLY/$ASSEMBLER/4_annotation/braker3
mkdir -p $BRAKER_DIR

###run braker3
#empty working directory
rm -r $BRAKER_DIR/*
#copy augustus config directory to working directory
cp -r ./config $BRAKER_DIR/config
#run container
singularity exec $braker3 braker.pl --AUGUSTUS_CONFIG_PATH=$BRAKER_DIR/config    \
                                    --species=$PROJECT                           \
                                    --workingdir=$BRAKER_DIR                     \
                                    --genome=$FASTA                              \
                                    --prot_seq=$ORTHO_DB                         \
                                    --threads=32                                 \
                                    --gff3

###run busco on annotation
#location downloaded lineages
downloads=/nfs/scratch/oostinto/scripts/genome_assembly/4_assembly_stats/busco_downloads
#create output dir
mkdir $BRAKER_DIR/busco
#run container
singularity run $busco busco    -i $BRAKER_DIR/braker.codingseq \
                                --out_path $BRAKER_DIR/busco	\
        				        -o busco_$ASSEMBLY'_'braker3    \
						        --offline 					    \
						        --download_path $downloads	    \
						        -l actinopterygii_odb10		    \
						        -f -m transcriptome			    \
						        --cpu 32

###run agat on annotation
#create output dir
rm -r $BRAKER_DIR/agat
#run container
singularity run $agat agat_sp_functional_statistics.pl -f $BRAKER_DIR/braker.gtf \
													   -g $FASTA				 \
													   -o $BRAKER_DIR/agat

################################### old code ################################################

##copy Augustus config folder to working directory
#cp -r ./config $BRAKER_DIR/config

#singularity exec -B $braker3 braker.pl --AUGUSTUS_CONFIG_PATH=$BRAKER_DIR/config                   \
#                                       --genome=$SOFTMASKED_FASTA     \
#                                       --prot_seq=/nfs/scratch/oostinto/braker_db/Vertebrata.fa    \
#                                       --workingdir=$BRAKER_DIR                                    \
#                                       --threads 32                                                \
#                                       --gff3

###doesn't run properly anymore with local installment of singularity 
##run braker
#cd $BRAKER_DIR
#$singularity exec $SCRATCH/singularity/braker3.sif braker.pl --genome $SOFTMASKED_FASTA    \
#                                                             --prot_seq /nfs/scratch/oostinto/braker_db/Vertebrata.fa    \
#                                                             --threads 32                                                \
#                                                             --gff3
##move results 
#mv braker/* .
#rm -r braker


##local install of singularity
#singularity=~/singularity/bin/singularity
##bind scratch and home to container
#export SINGULARITY_BIND="/nfs/scratch/oostinto,/nfs/home/oostinto"
#export SINGULARITY_TMPDIR="/nfs/scratch/oostinto/tmp"

#DESCRIPTION
#
#braker.pl   Pipeline for predicting genes with GeneMark-EX and AUGUSTUS with
#            RNA-Seq and/or proteins
#
#SYNOPSIS
#
#braker.pl [OPTIONS] --genome=genome.fa {--bam=rnaseq.bam | --prot_seq=prot.fa}
#
#INPUT FILE OPTIONS
#
#--genome=genome.fa                  fasta file with DNA sequences
#--bam=rnaseq.bam                    bam file with spliced alignments from
#                                    RNA-Seq
#--prot_seq=prot.fa                  A protein sequence file in multi-fasta
#                                    format used to generate protein hints.
#                                    Unless otherwise specified, braker.pl will
#                                    run in "EP mode" or "ETP mode which uses 
#                                    ProtHint to generate protein hints and 
#                                    GeneMark-EP+ or GeneMark-ETP to
#                                    train AUGUSTUS.
#--hints=hints.gff                   Alternatively to calling braker.pl with a
#                                    bam or protein fasta file, it is possible to
#                                    call it with a .gff file that contains
#                                    introns extracted from RNA-Seq and/or
#                                    protein hints (most frequently coming
#                                    from ProtHint). If you wish to use the
#                                    ProtHint hints, use its
#                                    "prothint_augustus.gff" output file.
#                                    This flag also allows the usage of hints
#                                    from additional extrinsic sources for gene
#                                    prediction with AUGUSTUS. To consider such
#                                    additional extrinsic information, you need
#                                    to use the flag --extrinsicCfgFiles to
#                                    specify parameters for all sources in the
#                                    hints file (including the source "E" for
#                                    intron hints from RNA-Seq).
#                                    In ETP mode, this option can be used together
#                                    with --geneMarkGtf and --traingenes to provide
#                                    BRAKER with results of a previous GeneMark-ETP
#                                    run, so that the GeneMark-ETP step can be
#                                    skipped. In this case, specify the hintsfile of
#                                    a previous BRAKER run here, or generate a
#                                    hintsfile from the GeneMark-ETP working
#                                    directory with the script get_etp_hints.py.
#--rnaseq_sets_ids=SRR1111,SRR1115   IDs of RNA-Seq sets that are either in
#                                    one of the directories specified with
#                                    --rnaseq_sets_dir, or that can be downloaded
#                                    from SRA. If you want to use local files, you
#                                    can use unaligned reads in FASTQ format
#                                    (they have to be named ID.fastq if unpaired or
#                                    ID_1.fastq, ID_2.fastq if paired), or aligned reads
#                                    as a BAM file (named ID.bam).
#--rnaseq_sets_dir=/path/to/rna_dir1 Locations where the local files of RNA-Seq data
#                                    reside that were specified with --rnaseq_sets_ids.
#
#FREQUENTLY USED OPTIONS
#
#--species=sname                     Species name. Existing species will not be
#                                    overwritten. Uses Sp_1 etc., if no species
#                                    is assigned
#--AUGUSTUS_ab_initio                output ab initio predictions by AUGUSTUS
#                                    in addition to predictions with hints by
#                                    AUGUSTUS
#--softmasking_off                   Turn off softmasking option (enables by 
#                                    default, discouraged to disable!)
#--esmode                            Run GeneMark-ES (genome sequence only) and
#                                    train AUGUSTUS on long genes predicted by
#                                    GeneMark-ES. Final predictions are ab initio
#--gff3                              Output in GFF3 format (default is gtf
#                                    format)
#--threads                           Specifies the maximum number of threads that
#                                    can be used during computation. Be aware:
#                                    optimize_augustus.pl will use max. 8
#                                    threads; augustus will use max. nContigs in
#                                    --genome=file threads.
#--workingdir=/path/to/wd/           Set path to working directory. In the
#                                    working directory results and temporary
#                                    files are stored
#--nice                              Execute all system calls within braker.pl
#                                    and its submodules with bash "nice"
#                                    (default nice value)
#--alternatives-from-evidence=true   Output alternative transcripts based on
#                                    explicit evidence from hints (default is
#                                    true).
#--fungus                            GeneMark-EX option: run algorithm with
#                                    branch point model (most useful for fungal
#                                    genomes)
#--crf                               Execute CRF training for AUGUSTUS;
#                                    resulting parameters are only kept for
#                                    final predictions if they show higher
#                                    accuracy than HMM parameters.
#--keepCrf                           keep CRF parameters even if they are not
#                                    better than HMM parameters
#--makehub                           Create track data hub with make_hub.py
#                                    for visualizing BRAKER results with the
#                                    UCSC GenomeBrowser
#--busco_lineage=lineage             If you provide a BUSCO lineage, BRAKER will
#                                    run compleasm on genome level to generate hints
#                                    from BUSCO to enhance BUSCO discovery in the
#                                    protein set. Also, if you provide a BUSCO
#                                    lineage, BRAKER will run compleasm to assess
#                                    the protein sets that go into TSEBRA combination,
#                                    and will determine the TSEBRA mode to maximize
#                                    BUSCO. Do not provide a busco_lineage if you
#                                    want to determina natural BUSCO sensivity of
#                                    BRAKER!
#--email                             E-mail address for creating track data hub
#--version                           Print version number of braker.pl
#--help                              Print this help message
#
#CONFIGURATION OPTIONS (TOOLS CALLED BY BRAKER)
#
#--AUGUSTUS_CONFIG_PATH=/path/       Set path to config directory of AUGUSTUS
#                                    (if not specified as environment
#                                    variable). BRAKER1 will assume that the
#                                    directories ../bin and ../scripts of
#                                    AUGUSTUS are located relative to the
#                                    AUGUSTUS_CONFIG_PATH. If this is not the
#                                    case, please specify AUGUSTUS_BIN_PATH
#                                    (and AUGUSTUS_SCRIPTS_PATH if required).
#                                    The braker.pl commandline argument
#                                    --AUGUSTUS_CONFIG_PATH has higher priority
#                                    than the environment variable with the
#                                    same name.
#--AUGUSTUS_BIN_PATH=/path/          Set path to the AUGUSTUS directory that
#                                    contains binaries, i.e. augustus and
#                                    etraining. This variable must only be set
#                                    if AUGUSTUS_CONFIG_PATH does not have
#                                    ../bin and ../scripts of AUGUSTUS relative
#                                     to its location i.e. for global AUGUSTUS
#                                    installations. BRAKER1 will assume that
#                                    the directory ../scripts of AUGUSTUS is
#                                    located relative to the AUGUSTUS_BIN_PATH.
#                                    If this is not the case, please specify
#                                    --AUGUSTUS_SCRIPTS_PATH.
#--AUGUSTUS_SCRIPTS_PATH=/path/      Set path to AUGUSTUS directory that
#                                    contains scripts, i.e. splitMfasta.pl.
#                                    This variable must only be set if
#                                    AUGUSTUS_CONFIG_PATH or AUGUSTUS_BIN_PATH
#                                    do not contains the ../scripts directory
#                                    of AUGUSTUS relative to their location,
#                                    i.e. for special cases of a global
#                                    AUGUSTUS installation.
#--BAMTOOLS_PATH=/path/to/           Set path to bamtools (if not specified as
#                                    environment BAMTOOLS_PATH variable). Has
#                                    higher priority than the environment
#                                    variable.
#--GENEMARK_PATH=/path/to/           Set path to GeneMark-ET (if not specified
#                                    as environment GENEMARK_PATH variable).
#                                    Has higher priority than environment
#                                    variable.
#--SAMTOOLS_PATH=/path/to/           Optionally set path to samtools (if not
#                                    specified as environment SAMTOOLS_PATH
#                                    variable) to fix BAM files automatically,
#                                    if necessary. Has higher priority than
#                                    environment variable.
#--PROTHINT_PATH=/path/to/           Set path to the directory with prothint.py.
#                                    (if not specified as PROTHINT_PATH
#                                    environment variable). Has higher priority
#                                    than environment variable.
#--DIAMOND_PATH=/path/to/diamond     Set path to diamond, this is an alternative
#                                    to NCIB blast; you only need to specify one
#                                    out of DIAMOND_PATH or BLAST_PATH, not both.
#                                    DIAMOND is a lot faster that BLAST and yields
#                                    highly similar results for BRAKER.
#--BLAST_PATH=/path/to/blastall      Set path to NCBI blastall and formatdb
#                                    executables if not specified as
#                                    environment variable. Has higher priority
#                                    than environment variable.
#--COMPLEASM_PATH=/path/to/compleasm Set path to compleasm (if not specified as
#                                    environment variable). Has higher priority
#                                    than environment variable.
#--PYTHON3_PATH=/path/to             Set path to python3 executable (if not
#                                    specified as envirnonment variable and if
#                                    executable is not in your $PATH).
#--JAVA_PATH=/path/to                Set path to java executable (if not
#                                    specified as environment variable and if
#                                    executable is not in your $PATH), only
#                                    required with flags --UTR=on and --addUTR=on
#--GUSHR_PATH=/path/to               Set path to gushr.py exectuable (if not
#                                    specified as an environment variable and if
#                                    executable is not in your $PATH), only required
#                                    with the flags --UTR=on and --addUTR=on
#--MAKEHUB_PATH=/path/to             Set path to make_hub.py (if option --makehub
#                                    is used).
#--CDBTOOLS_PATH=/path/to            cdbfasta/cdbyank are required for running
#                                    fix_in_frame_stop_codon_genes.py. Usage of
#                                    that script can be skipped with option
#                                    '--skip_fixing_broken_genes'.
#
#
#EXPERT OPTIONS
#
#--augustus_args="--some_arg=bla"    One or several command line arguments to
#                                    be passed to AUGUSTUS, if several
#                                    arguments are given, separate them by
#                                    whitespace, i.e.
#                                    "--first_arg=sth --second_arg=sth".
#--skipGeneMark-ES                   Skip GeneMark-ES and use provided
#                                    GeneMark-ES output (e.g. provided with
#                                    --geneMarkGtf=genemark.gtf)
#--skipGeneMark-ET                   Skip GeneMark-ET and use provided
#                                    GeneMark-ET output (e.g. provided with
#                                    --geneMarkGtf=genemark.gtf)
#--skipGeneMark-EP                   Skip GeneMark-EP and use provided
#                                    GeneMark-EP output (e.g. provided with
#                                    --geneMarkGtf=genemark.gtf)
#--skipGeneMark-ETP                  Skip GeneMark-ETP and use provided
#                                    GeneMark-ETP output (e.g. provided with
#                                    --gmetp_results_dir=GeneMark-ETP/)
#--geneMarkGtf=file.gtf              If skipGeneMark-ET is used, braker will by
#                                    default look in the working directory in
#                                    folder GeneMarkET for an already existing
#                                    gtf file. Instead, you may provide such a
#                                    file from another location. If geneMarkGtf
#                                    option is set, skipGeneMark-ES/ET/EP/ETP is
#                                    automatically also set. Note that gene and
#                                    transcript ids in the final output may not
#                                    match the ids in the input genemark.gtf
#                                    because BRAKER internally re-assigns these
#                                    ids.
#                                    In ETP mode, this option hast to be used together
#                                    with --traingenes and --hints to provide BRAKER
#                                    with results of a previous GeneMark-ETP run.
#--gmetp_results_dir                 Location of results from a previous
#                                    GeneMark-ETP run, which will be used to
#                                    skip the GeneMark-ETP step. This option
#                                    can be used instead of --geneMarkGtf,
#                                    --traingenes, and --hints to skip GeneMark.
#--rounds                            The number of optimization rounds used in
#                                    optimize_augustus.pl (default 5)
#--skipAllTraining                   Skip GeneMark-EX (training and
#                                    prediction), skip AUGUSTUS training, only
#                                    runs AUGUSTUS with pre-trained and already
#                                    existing parameters (not recommended).
#                                    Hints from input are still generated.
#                                    This option automatically sets
#                                    --useexisting to true.
#--useexisting                       Use the present config and parameter files
#                                    if they exist for 'species'; will overwrite
#                                    original parameters if BRAKER performs
#                                    an AUGUSTUS training.
#--filterOutShort                    It may happen that a "good" training gene,
#                                    i.e. one that has intron support from
#                                    RNA-Seq in all introns predicted by
#                                    GeneMark-EX, is in fact too short. This flag
#                                    will discard such genes that have
#                                    supported introns and a neighboring
#                                    RNA-Seq supported intron upstream of the
#                                    start codon within the range of the
#                                    maximum CDS size of that gene and with a
#                                    multiplicity that is at least as high as
#                                    20% of the average intron multiplicity of
#                                    that gene.
#--skipOptimize                      Skip optimize parameter step (not
#                                    recommended).
#--skipIterativePrediction           Skip iterative prediction in --epmode (does
#                                    not affect other modes, saves a bit of runtime)
#--skipGetAnnoFromFasta              Skip calling the python3 script
#                                    getAnnoFastaFromJoingenes.py from the
#                                    AUGUSTUS tool suite. This script requires
#                                    python3, biopython and re (regular
#                                    expressions) to be installed. It produces
#                                    coding sequence and protein FASTA files
#                                    from AUGUSTUS gene predictions and provides
#                                    information about genes with in-frame stop
#                                    codons. If you enable this flag, these files
#                                    will not be produced and python3 and
#                                    the required modules will not be necessary
#                                    for running brkaker.pl.
#--skip_fixing_broken_genes          If you do not have python3, you can choose
#                                    to skip the fixing of stop codon including
#                                    genes (not recommended).
#--eval=reference.gtf                Reference set to evaluate predictions
#                                    against (using evaluation scripts from GaTech)
#--eval_pseudo=pseudo.gff3           File with pseudogenes that will be excluded
#                                    from accuracy evaluation (may be empty file)
#--AUGUSTUS_hints_preds=s            File with AUGUSTUS hints predictions; will
#                                    use this file as basis for UTR training;
#                                    only UTR training and prediction is
#                                    performed if this option is given.
#--flanking_DNA=n                    Size of flanking region, must only be
#                                    specified if --AUGUSTUS_hints_preds is given
#                                    (for UTR training in a separate braker.pl
#                                    run that builds on top of an existing run)
#--verbosity=n                       0 -> run braker.pl quiet (no log)
#                                    1 -> only log warnings
#                                    2 -> also log configuration
#                                    3 -> log all major steps
#                                    4 -> very verbose, log also small steps
#--downsampling_lambda=d             The distribution of introns in training
#                                    gene structures generated by GeneMark-EX
#                                    has a huge weight on single-exon and
#                                    few-exon genes. Specifying the lambda
#                                    parameter of a poisson distribution will
#                                    make braker call a script for downsampling
#                                    of training gene structures according to
#                                    their number of introns distribution, i.e.
#                                    genes with none or few exons will be
#                                    downsampled, genes with many exons will be
#                                    kept. Default value is 2.
#                                    If you want to avoid downsampling, you have
#                                    to specify 0.
#--checkSoftware                     Only check whether all required software
#                                    is installed, no execution of BRAKER
#--nocleanup                         Skip deletion of all files that are typically not
#                                    used in an annotation project after
#                                    running braker.pl. (For tracking any
#                                    problems with a braker.pl run, you
#                                    might want to keep these files, therefore
#                                    nocleanup can be activated.)
#
#
#DEVELOPMENT OPTIONS (PROBABLY STILL DYSFUNCTIONAL)
#
#--splice_sites=patterns             list of splice site patterns for UTR
#                                    prediction; default: GTAG, extend like this:
#                                    --splice_sites=GTAG,ATAC,...
#                                    this option only affects UTR training
#                                    example generation, not gene prediction
#                                    by AUGUSTUS
#--overwrite                         Overwrite existing files (except for
#                                    species parameter files) Beware, currently
#                                    not implemented properly!
#--extrinsicCfgFiles=file1,file2,... Depending on the mode in which braker.pl
#                                    is executed, it may require one ore several
#                                    extrinsicCfgFiles. Don't use this option
#                                    unless you know what you are doing!
#--stranded=+,-,+,-,...              If UTRs are trained, i.e.~strand-specific
#                                    bam-files are supplied and coverage
#                                    information is extracted for gene prediction,
#                                    create stranded ep hints. The order of
#                                    strand specifications must correspond to the
#                                    order of bam files. Possible values are
#                                    +, -, .
#                                    If stranded data is provided, ONLY coverage
#                                    data from the stranded data is used to
#                                    generate UTR examples! Coverage data from
#                                    unstranded data is used in the prediction
#                                    step, only.
#                                    The stranded label is applied to coverage
#                                    data, only. Intron hints are generated
#                                    from all libraries treated as "unstranded"
#                                    (because splice site filtering eliminates
#                                    intron hints from the wrong strand, anyway).
#--optCfgFile=ppx.cfg                Optional custom config file for AUGUSTUS
#                                    for running PPX (currently not
#                                    implemented)
#--grass                             Switch this flag on if you are using braker.pl
#                                    for predicting genes in grasses with
#                                    GeneMark-EX. The flag will enable
#                                    GeneMark-EX to handle GC-heterogenicity
#                                    within genes more properly.
#                                    NOTHING IMPLEMENTED FOR GRASS YET!
#--transmasked_fasta=file.fa         Transmasked genome FASTA file for GeneMark-EX
#                                    (to be used instead of the regular genome
#                                    FASTA file).
#--min_contig=INT                    Minimal contig length for GeneMark-EX, could
#                                    for example be set to 10000 if transmasked_fasta
#                                    option is used because transmasking might
#                                    introduce many very short contigs.
#--translation_table=INT             Change translation table from non-standard
#                                    to something else.
#                                    DOES NOT WORK YET BECAUSE BRAKER DOESNT
#                                    SWITCH TRANSLATION TABLE FOR GENEMARK-EX, YET!
#--gc_probability=DECIMAL            Probablity for donor splice site pattern GC
#                                    for gene prediction with GeneMark-EX,
#                                    default value is 0.001
#--gm_max_intergenic=INT             Adjust maximum allowed size of intergenic
#                                    regions in GeneMark-EX. If not used, the value
#                                    is automatically determined by GeneMark-EX.
#--traingenes=file.gtf               Training genes that are used instead of training
#                                    genes generated with GeneMark.
#                                    In ETP mode, this option can be used together
#                                    with --geneMarkGtf and --hints to provide BRAKER
#                                    with results of a previous GeneMark-ETP run, so
#                                    that the GeneMark-ETP step can be skipped.
#                                    In this case, use training.gtf from that run as
#                                    argument.
#--UTR=on                            create UTR training examples from RNA-Seq
#                                    coverage data; requires options
#                                    --bam=rnaseq.bam.
#                                    Alternatively, if UTR parameters already
#                                    exist, training step will be skipped and
#                                    those pre-existing parameters are used.
#                                    DO NOT USE IN CONTAINER!
#                                    TRY NOT TO USE AT ALL!
#--addUTR=on                         Adds UTRs from RNA-Seq coverage data to
#                                    augustus.hints.gtf file. Does not perform
#                                    training of AUGUSTUS or gene prediction with
#                                    AUGUSTUS and UTR parameters.
#                                    DO NOT USE IN CONTAINER!
#                                    TRY NOT TO USE AT ALL!
#
#
#EXAMPLE
#
#To run with RNA-Seq
#
#braker.pl [OPTIONS] --genome=genome.fa --species=speciesname \
#    --bam=accepted_hits.bam
#braker.pl [OPTIONS] --genome=genome.fa --species=speciesname \
#    --hints=rnaseq.gff
#
#To run with protein sequences
#
#braker.pl [OPTIONS] --genome=genome.fa --species=speciesname \
#    --prot_seq=proteins.fa
#braker.pl [OPTIONS] --genome=genome.fa --species=speciesname \
#    --hints=prothint_augustus.gff
#
#To run with RNA-Seq and protein sequences
#
#braker.pl [OPTIONS] --genome=genome.fa --species=speciesname \
#    --prot_seq=proteins.fa --rnaseq_sets_ids=id_rnaseq1,id_rnaseq2 \
#    --rnaseq_sets_dir=/path/to/local/rnaseq/files 
#braker.pl [OPTIONS] --genome=genome.fa --species=speciesname \
#    --prot_seq=proteins.fa --bam=id_rnaseq1.bam,id_rnaseq2.bam 
