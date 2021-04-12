# m6A-express: uncovering complex and condition-specific m6A regulation of gene expression
m6A-express is a tool to uncover the complex and condition-specific m6A regulation of gene expression

# Pre-installing m6A-express
Before installing m6A-express package, user should have installed python 2.7 software in their platform (Windows systerm or Linux systerm). And some python packages should be installed, such as numpy, statsmodels.api, pandas, scipy. 

# Installation Instructions
Firstly, we need to install exomePeak package to do the peak calling for m6A methylation site
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install(c("SummarizedExperiment","cqn","Rsamtools",
                       "GenomicAlignments","GenomicRanges","GenomicFeatures",
                       "DESeq2","ggplot2","mclust",
                       "genefilter","BSgenome","BiocParallel",
                       "IRanges","S4Vectors","quantreg",
                       "reshape2","rtracklayer","apeglm","RMariaDB"))

> if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
> devtools::install_github("ZW-xjtlu/exomePeak")

Then, we installed the QNB package to detect the differential m6A methylation sites.
> devtools::install_github("NWPU-903PR/QNB")

Installed the reticulate pacakge to call python code in R
> install.packages("reticulate")

The m6A-express package can be installed by the following R commands:
> devtools::install_github("NWPU-903PR/m6Aexpress")
The m6A-express package can be installed by the following R commands:
> devtools::install_github("NWPU-903PR/m6Aexpress")

> library(m6Aexpress)

# Usage Example
The following command code will show how to use this package and output m6A methylation regulated expression gene in excel files
## Basic mode: pooled samples from one or multiple conditions together and identify m6A regulated expression gene in a specific context.
### Input BAM files.
> f1 <- system.file("extdata", "IP1.bam", package="m6Aexpress")

> f2 <- system.file("extdata", "IP2.bam", package="m6Aexpress")

> f3 <- system.file("extdata", "IP3.bam", package="m6Aexpress")

> f4 <- system.file("extdata", "IP4.bam", package="m6Aexpress")

> f5 <- system.file("extdata", "Input1.bam", package="m6Aexpress")
 
> f6 <- system.file("extdata", "Input2.bam", package="m6Aexpress")

> f7 <- system.file("extdata", "Input3.bam", package="m6Aexpress")
 
> f8 <- system.file("extdata", "Input4.bam", package="m6Aexpress")

> IP\_BAM <- c(f1,f2,f3,f4)

> INPUT\_BAM <- c(f5,f6,f7,f8)

### We use GTF file in the following example.
> gtf <- system.file("extdata", "hg19toy.gtf", package="m6Aexpress")

### Predict m6A regulated expression gene by m6A-express model
> m6A_reg\_exp\_gene <- m6Aexpress(express_data=INPUT_BAM, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, annot_type="hg19", GENE_ANNO_GTF=gtf, pvalue=0.05,mode="basic")

## DE-DM mode: In this case, we will detect whether the differential m6A methylation peak sites regulated expression that caused the differential expression genes
### Predict the differential expression genes are regulated by differential methylation peak sites by m6A-express model

> IP\_BAM <- c(f1,f2)

> TREATED\_IP\_BAM <- c(f3,f4)

> INPUT\_BAM <- c(f5,f6)

> TREATED\_INPUT\_BAM <- c(f7,f8)

> m6A\_reg\_exp\_gene <- m6Aexpress(express_data=INPUT_BAM, treated_express_data=TREATED_INPUT_BAM, IP_BAM=IP_BAM, TREATED_IP_BAM=TREATED_IP_BAM, INPUT_BAM=INPUT_BAM, TREATED_INPUT_BAM=TREATED_INPUT_BAM,annot_type="hg19", GENE_ANNO_GTF=gtf,pvalue=0.05,mode="DE-DM")
