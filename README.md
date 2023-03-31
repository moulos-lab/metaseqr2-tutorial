# Transcriptomics practicals

This small tutorial covers a short hands-on on basic transcriptomics (RNA-Seq)
data analysis with the Bioconductor package [metaseqR2](https://www.bioconductor.org/packages/release/bioc/html/metaseqR2.html). 
metaseqR2 implements an RNA-Seq data statistical analysis pipeline by combining 
the p-value outcomes from several individual statistical tests. It has been 
shown that the outcome is more accurate than using individual tests alone. At 
the same time, metaseqR2 offers an interface to several normalization methods 
and implements many gene filtering options. Finally, it presents the results in 
an interactive HTML report with a rich variety of QC and statistical analysis 
visualization graphs.

More information about metaseqR2 can be found in the respective publications:

> Dionysios Fanidis, Panagiotis Moulos: **Integrative, normalization-insusceptible statistical analysis of RNA-Seq data, with improved differential expression and unbiased downstream functional analysis**, *Briefings in Bioinformatics*, 2020, bbaa156, DOI: [10.1093/bib/bbaa156](https://doi.org/10.1093/bib/bbaa156)

> Panagiotis Moulos, Pantelis Hatzis: **Systematic integration of RNA-Seq statistical algorithms for accurate detection of differential gene expression patterns**, *Nucleic Acids Research*, 2015, 43(4):e25, DOI: [10.1093/nar/gku1273](https://doi.org/10.1093/nar/gku1273)

Public examples of metaseqR2 reports can be found [here](http://epigenomics.fleming.gr/~panos/metaseqR2_showcase/metaseqR2_DOX_vs_CON_3I4_76_PANDORA),
[here](http://epigenomics.fleming.gr/~panos/metaseqR2_showcase/ metaseqR2_LiverDevelopment_PANDORA)
and [here](http://epigenomics.fleming.gr/~panos/metaseqR2_showcase/metaseqR2_Smyd3_PANDORA).

## Prerequisites

* A Linux operating system and command-line environment
* [R](https://www.r-project.org/), at least version 4.0
* [Bioconductor](https://www.bioconductor.org/) and [metaseqR2](https://www.bioconductor.org/packages/release/bioc/html/metaseqR2.html) installed

Optionally, some may prefer to work through [RStudio](https://posit.co/download/rstudio-desktop/), 
so this must also be installed, but it is not required.

## RNA-Seq data

Because of time restrictions, we will not examine the step of QC or alignment to
the reference genome. We will use existing aligned BAM files which are part of
the [SeqCVIBE](http://elixir-seqcvibe.hybridstat.gr/seqcvibe) RNA-Seq data 
analysis and exploration tool. If anybody is interested to learn more about 
SeqCVIBE, information can be found in the following publication:

> Efthimios Bothos, Pantelis Hatzis, Panagiotis Moulos: **Interactive Analysis, Exploration, and Visualization of RNA-Seq Data with SeqCVIBE**, *Methods and Protcols*, 2022, 5(2):27, DOI: [10.3390/mps5020027](https://doi.org/10.3390/mps5020027)

The data that we are going to examine have been previously described in 
[this](https://www.frontiersin.org/articles/10.3389/fmicb.2021.760627/full) 
study. Essentially, the authors use RNA-Seq in human brain microvascular 
endothelial cells to study gene expression changes when cells are challenged 
with the bacterium *Borrelia burgdoferi* which is responsible for 
[Lyme disease](https://www.cdc.gov/lyme/index.html). The cells are also treated
with its ligand, Erp23. More information about the data can be found in the
respective [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress)
[entry](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8053). The
analysis protocol mentioned there is not the same as the one mentioned in this
tutorial. For example, the reference genome here is the *Human GRCh37* or *hg19*
version and the reads have been aligned with a 2-step alignment process using
[HiSat2](http://daehwankimlab.github.io/hisat2/) and 
[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/). More information can be
found in the SeqCVIBE publication above.

The dataset has three biological conditions:
- Wild Type (WT)
- Challenged with B.burgdoferi (BB)
- Treated with Erp23 ligand (ERP)

## Working environment

Whether you are working in RStudio or not, it is assumed that all operations
described in this tutorial will be performed in a single directory within your
`HOME` directory in your Linux environment. For example, if your working
directory is `/home/johndoe`, then all the following steps assume that `HOME`
will be replaced with `/home/johndoe`.

# Tutorial

The following sections present step-by-step instructions for data retrieval,
creation of a text files with the targets (the BAM files with conditions) for
metaseqR2 and the execution of the pipeline.

## Working directory

Let's start by creating a working directory in `HOME` where all the operations
will take place.

```
HOME=/home/johndoe # Please replace this accordingly!
cd $HOME
mkdir tutorial && cd tutorial
```

## Data retrieval

We are using the Linux `wget` command to retrieve the BAM files and their 
indexes, it will take some time but not too much. We are in the `tutorial` 
directory.

```
mkdir bam && cd bam

# WT conditions
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367907.bam -O WT_1.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367907.bam.bai -O WT_1.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367908.bam -O WT_2.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367908.bam.bai -O WT_2.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367909.bam -O WT_3.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367909.bam.bai -O WT_3.bam.bai

# BB condition
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367901.bam -O BB_1.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367901.bam.bai -O BB_1.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367902.bam -O BB_2.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367902.bam.bai -O BB_2.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367903.bam -O BB_3.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367903.bam.bai -O BB_3.bam.bai

# ERP condition
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367904.bam -O ERP_1.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367904.bam.bai -O ERP_1.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367905.bam -O ERP_2.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367905.bam.bai -O ERP_2.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367906.bam -O ERP_3.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367906.bam.bai -O ERP_3.bam.bai

cd ..
```

## Targets file

The targets file required by metaseqR2 pipeline is a simple text tab-delimited
file describing the samples and experimental conditions. It has the following
columns (some are required and some are optional).

1. *Sample name*: this column contains a **unique** name for each sample in the
dataset, that is one **unique** name for each BAM file. It must not contain
dashes (-).
2. *BAM file path*: the full path and filename of the BAM file corresponding to
the sample in the first column.
3. *Condition*: The biological condition corresponding to the sample. It must 
not contain dashes (-).
4. *Paired-end reads*: Does the data are in paired-end reads or single-end? If
paired, it should be `paired` otherwise `single`. This column is optional and if 
not provided, reads will be cosnidered single-end.
5. *Strand direction*: If the sequencing protocol is forward-strand based, it
should be `forward` (the default), else `reverse`. This column is optional, if
not provided, it will be treated as forward.

The file *must* have a header, i.e. a line with column titles. Title names are
not important but they must be there. You can find the targets file for this
tutorial in the `files` directory of this repository. You must change the `HOME`
text with your working directory as mentioned above.

```
# Always in the tutorial directory
wget https://github.com/moulos-lab/metaseqr2-tutorial/raw/main/files/targets.txt
```

## Annotation

In order to determine the number of reads accumulated on each genomic region of
interest (in the case of typical RNA-Seq, the exons of the 
[canonical transcript](https://www.ensembl.org/info/genome/genebuild/canonical.html)), 
we need a set of genomic co-ordinates. These co-ordinates are used by counting
algorithms to assign reads falling on them. The metaseqR2 pipeline does this
"automatically" with two possible ways:

* By downloading the required annotation on-the-fly according to the respective
command arguments (organism, gene or exon counting etc.)
* By querying a pre-created database of supported annotations for speed

If we wish to analyze data from a non-supported organism, metaseqR2 can also
accept a [GTF](https://www.ensembl.org/info/website/upload/gff.html) describing
the genome under investigation. We will not cover this case here.

You can download the database as follows:

```
# Always in the tutorial directory
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/annotation.sqlite
```

## Performing differential expression analysis

metaseqR2 offers a complete end-to-end pipeline for 

* Read counting
* Data QC filtering
* Normalization
* Differential expression analysis
* Reporting and visualization

As such, the central command used to run the pipeline has many arguments that
control input options, filtering, differential expression and what will be the
output. Not all of them can be covered here but complete explanations can be
found in the package documentation.

### Preliminaries

Firstly start R and load the metaseqR2 library:

```
library(metaseqR2)
```

All the exported functions (functions that can be used by the R session user)
are documented and help can be found:

```
help(metaseqr2)
# or ?metaseqr2

help(normalizeDeseq2)
# or ?normalizeDeseq2
```

### Define the targets file

We create a variable pointing to the targets file we created:

```
targetsFile <- file.path(HOME,"tutorial","targets.txt")
```

### Define the comparisons

Let's perform three comparisons of the "treatments":

* Challenge with *Borrelia burgdoferi* against Wild Type (`BB_vs_WT`)
* Treatment with its Erp23 ligand (`ERP_vs_WT`)
* Treatment with the Erp23 ligand agains *Borrelia burgdoferi* (`ERP_vs_BB`)

```
theContrasts <- c(
    "BB_vs_WT",
    "ERP_vs_WT",
    "ERP_vs_BB"
)
```

### Run the metaseqR2 pipeline

You can paste the following command in R's environment which will execute a
metaseqR2 pipeline:

```
metaseqr2(
    sampleList=targetsFile,
    contrast=theContrasts,
    org="hg19",
    countType="exon",
    normalization="deseq2",
    statistics="deseq2",
    figFormat="png",
    qcPlots=c(
        "mds","biodetection","countsbio","saturation","readnoise","filtered",
        "correl","boxplot","meandiff","meanvar","deheatmap","volcano","mastat",
        "foldvenn"
    ),
    exportWhere=file.path(HOME,"tutorial","analysis_1"),
    pcut=0.05,
    restrictCores=0.25,
    exportWhat=c("annotation","p_value","adj_p_value","fold_change",
        "counts","flags"),
    exportScale=c("natural","log2","rpgm"),
    exportValues="normalized",
    saveGeneModel=TRUE,
    reportTop=0.1,
    localDb=file.path(HOME,"tutorial","annotation.sqlite")
)
```

After finishing, you should be able to see a new directory named `analysis_1`
in your `HOME/tutorial` directory. Before inspecting the report, let's take a
look at each parameter in the command above:

* `sampleList=targetsFile`: The `sampleList` argument can be a list defining the
conditions and samples in the experiments (please look at metaseqR2 
documentation for its definition). It can also be a file defining these samples
(the targets file). If it's a list, it should be accompanied by a pre-calculated
count matrix and not a targets file with the BAM files. In this case, we are
giving the targets file we defined above.
* `contrast=theContrasts`: This defines the comparisons to be performed. It has
a special syntax, specifically condition names separated by `_vs_`. For pairwise
comparisons, it should be `ConditonA_vs_ConditionB`. Regarding fold change
calculation (ratios of gene expression between two conditions), `ConditionB` is
always the denominator. So, `BB_vs_WT` means that the algorithm will seek 
differentially expressed genes betwee `BB` and `WT` and the fold change (up- or
down-regulation) is `BB/WT`.



# Bonus! BAM stats

##















