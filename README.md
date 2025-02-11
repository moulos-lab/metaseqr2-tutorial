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

If you are running in a Windows system, you should create a folder called `bam`
in the `tutorial` folder and then open R (e.g. RStudio) and navigate to this
folder. Then you can use the `download.file()` R function:

```
# Avoid download timeout
options(timeout=max(3600,getOption("timeout")))

# WT conditions
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367907.bam","WT_1.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367907.bam.bai","WT_1.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367908.bam","WT_2.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367908.bam.bai","WT_2.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367909.bam","WT_3.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/WT/ERR3367909.bam.bai","WT_3.bam.bai")

# BB condition
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367901.bam","BB_1.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367901.bam.bai","BB_1.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367902.bam","BB_2.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367902.bam.bai","BB_2.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367903.bam","BB_3.bam") 
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/B.burgdorferi/ERR3367903.bam.bai","BB_3.bam.bai") 

# ERP condition
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367904.bam","ERP_1.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367904.bam.bai","ERP_1.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367905.bam","ERP_2.bam") 
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367905.bam.bai","ERP_2.bam.bai") 
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367906.bam","ERP_3.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/PRJEB33005/Erp23/ERR3367906.bam.bai","ERP_3.bam.bai")
```

Then, within R, go to up one level to the `tutorial` directory:

```
setwd("../")
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

# or within R for Windows
download.file("https://github.com/moulos-lab/metaseqr2-tutorial/raw/main/files/targets.txt","targets.txt")
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

# or within R
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/annotation.sqlite","annotation.sqlite")
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

#### Basic run

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
* `org="hg19"`: A character denoting which organism should metaseqR2 consider as 
the reference organism (and genome version). See the documentation for the 
supported organisms and how they are provided. In this case it's the human 
genome version "hg19".
* `countType="exon"`: Where should the read counting take place? This depends on
the RNA-Seq protocol. The classic protocol measures the abundance of 
polyadenylated, that is transcribed RNAs (mRNAs). So counting should take place
over exons, hence `countType="exon"`. If the protocol measures total RNA, then
the counting should take place over the whole gene body. See the documentation
for a list of available modes.
* `normalization="deseq2"`: Which normalization algorithm should be applied?
* `statistics="deseq2"`: Which of the supported differential expression tests
should be performed?
* `figFormat="png"`: metaseqR2 creates a directory with the run results, among
which a lot are figures. Apart from their interactive form in the report, they
can be exported in various formats. If we don't specify which, all of them will
be exported (we usually don't want that).
* `qcPlots=c(...)`: Which diagnostic plots to report. Some combinations of plots
and comparisons are not feasible. metaseqR2 detects them, warns the user and
then removes them from the list of plots.
* `exportWhere=...`: A directory to write the pipeline results.
* `pcut=0.05`: A p-value cutoff to report statistically significant results.
These are written in text files in the results directory and are downloadable
from the report. Note that the interactive table in the report contains only
the top scoring genes (controlled by the `reportTop` variable below).
* `restrictCores=0.25`: If the pipeline is running in a system with multiple
cores, the how many to use for parallel computations, when possible (e.g. 
counting). Not available on Windows.
* `exportWhat=c(...)`: The pipeline calculates *a lot* of possible values with
respect to gene expression. The most important ones are annotation elements,
statistical scores (p-values and FDRs) and fold changes. See documentation for
all possibilities.
* `exportScale=c(...)`: What scale to use? Natural (read counts as they are),
log<sub>2</sub> transformed?
* `exportValues="normalized"`: How should all the metrics be calculated? Based
on raw on normalized read counts. Both options are available at the same time.
* `saveGeneModel=TRUE`: This saves the counting and exon summarization results
in an R data file which can be used for additional analyses and re-analyses.
* `reportTop=0.1`: See `pcut` above.
* `localDb=...`: The annotation database. If not provided, metaseqR2 will look
at a default location and if not found, annotation will be downloaded 
on-the-fly.

Detailed lists of differentially expressed genes can be downloaded from the
*Results* section of the report. These are tab-delimited text files which can
be opened in any editor and also Excel. You can download a formatted example
from the `files` directory of this repository or 
[here](https://github.com/moulos-lab/metaseqr2-tutorial/raw/main/files/metaseqr_sig_out_BB_vs_WT.xlsx).

#### Run for a non-pairwise comparison

Apart from pairwise comparisons, metaseqR2 can be used for ANOVA-like 
comparisons where we look for genes that differ in at least one of more than two
conditions. Not all included tests support this. For this reason, metaseqR2 does
not currently support more complex statistical designs.

We execute another run with an ANOVA-like design with less plots (not needed 
since we created them before) and using the gene model we just created to avoid 
many recalculations (counting and summarizing):

```
# Run with an ANOVA-like design (genes that differ in at least one condition)
# with less plots (not needed since we created them before) and using the gene
# model we just created to avoid many recalculations (counting)
theContrasts <- "ERP_vs_BB_vs_WT"

# The gene model is provided through the counts argument
geneModelFile <- file.path(HOME,"tutorial","analysis_1","data",
    "gene_model.RData")

# Run
metaseqr2(
    counts=geneModelFile,
    contrast=theContrasts,
    org="hg19",
    countType="exon",
    normalization="deseq2",
    statistics="deseq2",
    figFormat="png",
    qcPlots="mds",
    exportWhere=file.path(HOME,"tutorial","analysis_2"),
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

#### Remove one sample from the analysis

We saw before that WT_2 clustered a bit differently with respect to correlation
with the other samples from WT condition. Let's see what happens when we remove
it:

```
# Let's see what will happen if we remove one sample (WT_2). This can be done
# by providing a list with the conditions and samples to be excluded. The name
# of the list member is the condition, while the member is the sample
excludeSample <- list(WT="WT_2")

# Define again the first comparisons
theContrasts <- c(
    "BB_vs_WT",
    "ERP_vs_WT",
    "ERP_vs_BB"
)

# Run
metaseqr2(
    counts=geneModelFile,
    excludeList=excludeSample,
    contrast=theContrasts,
    org="hg19",
    countType="exon",
    normalization="edger",
    statistics="edger",
    figFormat="png",
    qcPlots=c("mds","boxplot","volcano"),
    exportWhere=file.path(HOME,"tutorial","analysis_3"),
    pcut=0.05,
    restrictCores=0.25,
    exportWhat=c("annotation","p_value","adj_p_value","fold_change","counts"),
    exportScale=c("natural","log2","rpgm"),
    exportValues="normalized",
    saveGeneModel=TRUE,
    reportTop=0.1,
    localDb=file.path(HOME,"tutorial","annotation.sqlite")
)
```

There is a new argument here:
* `excludeList`: This instructs the pipeline to ignore a sample either from a
precalculated gene model data file, or from a targets file.

We can now inspect the report.

#### The PANDORA algorithm

The power of metaseqR2 is the PANDORA algorithm. Let's try it with 5 out of 9
supported tests. The report adjusts accordingly by reporting the numbers 
returned by each algorithm separately as well as the results obtained with
PANDORA.

PANDORA works with precalculated p-value weights for each organism, available
within the package. Then we instruct metaseqR2 to calculate the PANDORA p-values
as adjusted p-values:

```
# We need to get the pre-calculated p-value weights for human
weights <- getWeights("human")

# The selection about calculating PANDORA p-values is controlled by the metaP
# argument.
panp <- "pandora"
```

Then, the run:

```
metaseqr2(
    counts=geneModelFile,
    contrast=theContrasts,
    org="hg19",
    countType="exon",
    normalization="deseq2",
    statistics=c("deseq2","edger","noiseq","limma","absseq"),
    metaP=panp,
    weight=weights,
    figFormat="png",
    qcPlots=c("mds","boxplot","volcano","statvenn"),
    exportWhere=file.path(HOME,"tutorial","analysis_4"),
    pcut=0.05,
    restrictCores=0.25,
    exportWhat=c("annotation","p_value","adj_p_value","meta_p_value",
        "adj_meta_p_value","fold_change","counts"),
    exportScale=c("natural","log2","rpgm"),
    exportValues="normalized",
    saveGeneModel=TRUE,
    reportTop=0.1,
    localDb=file.path(HOME,"tutorial","annotation.sqlite")
)
```

Notice the two additional values in `exportWhat`:

* `"meta_p_value"`: This will export the PANDORA p-values, as just `p_value`
will only export the p-values from each test performed.
* `"adj_meta_p_value"`: FDR based on the PANDORA p-values.

Let's explore the report.

You can download an Excel file with the PANDORA results from the `files` 
directory of this repository or 
[here](https://github.com/moulos-lab/metaseqr2-tutorial/raw/main/files/metaseqr_sig_out_BB_vs_WT_PANDORA.xlsx)

## Bonus! BAM stats

As a bonus, if time left, we can calculate read count statistics from BAM files
that help assess the overall quality of the mapping. For this, we need the 
targets file we used for the differential expression analysis and a small 
library hosted on GitHub. We also (optionally) need a description of the main
target areas (in this case exons) to calculate stats over these. This
description is essentially a list pointing to a metaseqR2 database. If there is
no database available, the description (`ref` argument) can be a 
[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file (just the first
three columns are enough).

```
library(parallel)
library(Rsamtools)
library(GenomicAlignments)
library(openxlsx)
library(metaseqR2)

# The script
source("https://github.com/moulos-lab/genomics-facility-processes/raw/main/bamstats.R")

# Reference annotation to use for getting stats at areas of interest (exons)
ref <- list(
    org="hg19",
    refdb="ensembl",
    version="auto", # Optional
    countType="exon",
    localDb=file.path(HOME,"tutorial","annotation.sqlite")
)

stats <- getBamStats(targetsFile,ref=ref,rc=0.25)
```

You can also get the output from 
[here](https://github.com/moulos-lab/metaseqr2-tutorial/raw/main/files/alignment_statistics_2023-03-31-19-57-03.xlsx)

