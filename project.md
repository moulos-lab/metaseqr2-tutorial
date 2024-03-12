# RNA-Seq data analysis project

The project you have to complete comprises two parts. The first part is to 
succesfully align one RNA-Seq FASTQ file with [HISAT2](http://daehwankimlab.github.io/hisat2/).
The second part is to perform differential expression analysis in a simple 
RNA-Seq dataset with
[metaseqR2](https://www.bioconductor.org/packages/release/bioc/html/metaseqR2.html).

## Part 1 - Alignment of a FASTQ file

You can retrieve the FASTQ file from [here](http://epigenomics.fleming.gr/~panos/appbio/human.fastq.gz). 
The reference genome is `hg19`.
You have to align the file with [HISAT2](http://daehwankimlab.github.io/hisat2/).
If the respective index is not set up in the environment that you will perform
the alignment, you will have to set it up based on the instructions given by
HISAT2 authors. You will find the required index [here](http://daehwankimlab.github.io/hisat2/download/#h-sapiens).
You would be looking for `hg19`. You need to download and place it in your
location of preference, as long as it can be accessed by the aligner software.
You will also need `samtools` to index the resulting BAM file (or create a BAM
file if you choose to output as SAM). You can find samtools
[here](https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2)

A possible setup is for example:

```
# Create a working directory
mkdir align && cd align

# Get the FASTQ file
wget http://epigenomics.fleming.gr/~panos/appbio/human.fastq.gz

# Get and uncompress the index
wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
tar -xvf hg19_genome.tar.gz

# Get and unzip the aligner itself
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1-Linux_x86_64.zip
unzip hisat2-2.2.1-Linux_x86_64.zip

# Test the installation
./hisat2-2.2.1/hisat2 --help

# Get samtools (if not already in the system you are working at)
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
tar -xvf samtools-1.17.tar.bz2
cd samtools-1.17
./configure && make
cd ..

# Test the installaiton
./samtools-1.17/samtools --help
```

The entry point to HISAT2, supposing we are in the `align` directory is:

```
./hisat2-2.2.1/hisat2
```

while for samtools, the entry point is:

```
./samtools-1.17/samtools
```

You will need to pass any additional parameters according to the manual in order
to perform the spliced alignment with HISAT2. You do *not* have to uncompress 
the FASTQ file. HISAT2 works also with zipped files.

**Hint**: HISAT2 does not produce directly BAM files. You will have to use
samtools to get the BAM. You can do this by:

```
samtools view -bS HISAT_OUT.sam > HISAT_OUT.bam
```

### Tasks

1. Align the file to the reference genome and produce a BAM file.
2. Provide the HISAT2 command used to generate the BAM file.
3. Provide the HISAT2 stats output from the screen. What is the total alignment
rate?

You do **not** have to deliver the BAM file, just the answers to points 2 and 3.

## Part 2 - QC and differential expression analysis

You will work with parts of a dataset described 
[here](https://pubmed.ncbi.nlm.nih.gov/25519703/). It examines the in vivo 
transcriptome profiling of pre- and post-treatment prostatic biopsies from 
patients with advanced hormone-naive prostate cancer treated with docetaxel 
chemotherapy. There are two conditions:

* Pre-docetaxel
* Post-docexatel

The analysis should be performed as explained during the course. There will be
a single comparison, Post vs Pre-docetaxel.

You can download the data using the following commands:

With Linux `wget`:

```
# Pre-docetaxel
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993713.bam  -O Pre_1.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993713.bam.bai -O Pre_1.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993715.bam -O Pre_2.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993715.bam.bai -O Pre_2.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993717.bam -O Pre_3.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993717.bam.bai -O Pre_3.bam.bai

# Post-docetaxel
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993720.bam -O Post_1.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993720.bam.bai -O Post_1.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993722.bam -O Post_2.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993722.bam.bai -O Post_2.bam.bai
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993724.bam -O Post_3.bam
wget http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993724.bam.bai -O Post_3.bam.bai
```

Within R:

```
# Avoid download timeout
options(timeout=max(3600,getOption("timeout")))

# Pre-docetaxel
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993713.bam","Pre_1.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993713.bam.bai","Pre_1.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993715.bam","Pre_2.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993715.bam.bai","Pre_2.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993717.bam","Pre_3.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Pre-Docetaxel/SRR993717.bam.bai","Pre_3.bam.bai")

# Post-docetaxel
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993720.bam","Post_1.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993720.bam.bai","Post_1.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993722.bam","Post_2.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993722.bam.bai","Post_2.bam.bai")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993724.bam","Post_3.bam")
download.file("http://elixir-seqcvibe.hybridstat.gr/seqc_elixir/data/GSE51005/Post-Docetaxel/SRR993724.bam.bai","Post_3.bam.bai")
```

### Tasks

1. Create and deliver the targets file. The BAM files contain **paired-end**
reads. Take this into account in your target file.
2. Perform a differential expression analysis with metaseqR2 for the contrast
Post vs Pre. The normalization method as well as the statistical testing method
should be `"deseq2"`. The following QC plots should be reported: `"mds",
"biodetection","countsbio","saturation","correl","boxplot","meandiff","meanvar",
"volcano","mastat"`. A report should be created. Please deliver the command you
used.
3. Inspect the report. Are there any samples that could be excluded in a second
round of analysis? Why?
4. Report the number of differentially expressed genes based on the report.
5. Create and deliver the alignment metrics using the R scripts and the 
metaseqR2 database as described in the last part of the tutorial, using the
targets file you created for the differential expression analysis.

You do **not** have to deliver the whole report, just what is mentioned in 
points 1-5.

## Deadline

The deadline for the project is **Friday, 15/03/2024, end of day**. The answers 
to the questions should be sent by e-mail to Panagiotis Moulos.
