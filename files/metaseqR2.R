# Load library
library(metaseqR2)

# Define main directory
HOME <- "/home/panos"

# Define targets file
targetsFile <- file.path(HOME,"tutorial","targets.txt")

# Define the statistical comparisons (pairwise)
theContrasts <- c(
    "BB_vs_WT",
    "ERP_vs_WT",
    "ERP_vs_BB"
)

# Run with all pairwise comparisons and many diagnostic plots
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


# Run with an ANOVA-like design (genes that differ in at least one condition)
# with less plots (not needed since we created them before) and using the gene
# model we just created to avoid many recalculations (counting)
theContrasts <- "ERP_vs_BB_vs_WT"

# The gene model is provided through the counts argument
geneModelFile <- file.path(HOME,"tutorial","analysis_1","data",
    "gene_model.RData")

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

# Let's see what will happen if we remove one sample (WT_2). This can be done
# by providing a list with the conditions and samples to be excluded. The name
# of the list member is the condition, while the member is the sample
excludeSample <- list(WT="WT_2")

theContrasts <- c(
    "BB_vs_WT",
    "ERP_vs_WT",
    "ERP_vs_BB"
)

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


# Let's now try the PANDORA algorithm with 5 out of 9 algorithms for just one
# contrast for the shake of time.
theContrasts <- "BB_vs_WT"

# We need to get the pre-calculated p-value weights for human
weights <- getWeights("human")

# The selection about calculating PANDORA p-values is controlled by the metaP
# argument.
panp <- "pandora"

# Let's run
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
