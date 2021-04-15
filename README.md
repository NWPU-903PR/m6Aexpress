# m6A-express: uncovering complex and condition-specific m6A regulation of gene expression
m6A-express is a tool to uncover the complex and condition-specific m6A regulation of gene expression

# Pre-installing m6A-express
Before installing m6A-express package, user should have installed python 2.7 software in their platform (Windows systerm or Linux systerm). And some python packages should be installed, such as numpy, statsmodels.api, pandas, scipy. 
In user's python IDE platform, e.g.Spyder or pycharm, load the following python packages.
> import numpy as np
> 
> import re
> 
> import sys
> 
> import numpy as np
> 
> import statsmodels.api as sm
> 
> import pandas as pd
> 
> from scipy.optimize import minimize
> 
> from scipy.optimize import minimize_scalar
> 
> from scipy.stats import nbinom
> 
> from scipy.special import digamma
> 
> import statsmodels.sandbox as sms
> 
# Installation Instructions
The version of current package was supported by R 3.5.3 or new version
Firstly, we need to install exomePeak package to do the peak calling for m6A methylation site

> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
>
> BiocManager::install(c("SummarizedExperiment","Rsamtools",
                       "GenomicAlignments","GenomicRanges","GenomicFeatures",
                       "DESeq2","ggplot2","mclust",
                       "genefilter","BSgenome","BiocParallel",
                       "IRanges","S4Vectors","quantreg",
                       "reshape2","rtracklayer","apeglm"))

> if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
>
> devtools::install_github("ZW-xjtlu/exomePeak")

Then, we installed the QNB package to detect the differential m6A methylation sites.

> install.packages("https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz", repos = NULL, type="source")

Installed the reticulate pacakge to call python code in R

> install.packages("reticulate")
> ##install miniconda to install specific python package
> install_miniconda()
> ##install specific python package in R
> py_install("statsmodels"); py_install("pandas); py_install("scipy"); py_install("numpy")

Before install the m6Aexpress package, you should install the following R package:

> BiocManager::install(c('org.Hs.eg.db','org.Mm.eg.db','Rsubread', 
                         'TxDb.Hsapiens.UCSC.hg19.knownGene','TxDb.Mmusculus.UCSC.mm10.knownGene',
                           'AnnotationDbi'))
                       
The m6A-express package can be installed by the following R commands:
> devtools::install_github("NWPU-903PR/m6Aexpress")

> library(m6Aexpress)

# Usage Example
The following command code will show how to use this package and output m6A methylation regulated expression gene in excel files. The input data for m6Aexpress package includes the INPUT and IP BAM files from MeRIP-seq data. The INPUT BAM files are used to quantify the gene expression under specific context. The IP BAM files with paired INPUT BAM files are used to quantify the methylation intensity for each gene in specific context. m6Aexpress model could detect the correlation between gene expression and methylation and predicated some gene sets, whose gene expression are significantly regulated by methylation in specific context. m6Aexpress can predicate m6A regulated expression gene (m6A-reg-exp) in differential expression and differential methylation context. m6Aexpress can also predicate m6A-re-exp genes in tissue-specific (multiple sub-groups) context, which considered high dynamic peak or high variable peak sites across sub-tissue and predicate m6A-reg-exp genes with high variable peak.
Overall, m6Aexpress can predicate m6A-reg-exp gene set in case-control context and tissue-specific context. The following will introduce how to use m6Aexpress package step by step or only on step to obtain significant m6A-reg-exp gene.
## Step by Step Analysis
### Peak calling for methylation sites in DE-DM context and obtain consisten peak sites
> f1 <- system.file("extdata", "IP1.bam", package="m6Aexpress")

> f2 <- system.file("extdata", "IP2.bam", package="m6Aexpress")

> f3 <- system.file("extdata", "IP3.bam", package="m6Aexpress")

> f4 <- system.file("extdata", "IP4.bam", package="m6Aexpress")

> f5 <- system.file("extdata", "Input1.bam", package="m6Aexpress")
 
> f6 <- system.file("extdata", "Input2.bam", package="m6Aexpress")

> f7 <- system.file("extdata", "Input3.bam", package="m6Aexpress")
 
> f8 <- system.file("extdata", "Input4.bam", package="m6Aexpress")

> IP\_BAM <- c(f1,f2)

> INPUT\_BAM <- c(f5,f6)

> TREATED\_IP\_BAM <- c(f3,f4)

> TREATED\_INPUT\_BAM <- c(f7,f8)

> #Input the gene annotation file 
> gtf <- system.file("extdata", "hg19toy.gtf", package="m6Aexpress")
 
> #Obtain the consistent peak sites
> Get_peak_infor <- Get_peak_sites(IP_BAM, INPUT_BAM,TREATED_IP_BAM, TREATED_INPUT_BAM, GENE_ANNO_GTF=gtf, species="human")
### Differential methylation analysis for the consistent peak sites in case-control context
> DM_sites_infor <- DM_detect(peak_inform=Get_peak_infor,DM_CUTOFF_TYPE="pvalue",num_ctl=2, diff_peak_pvalue=0.05)
### Calculate the methylation intensity for each gene with DM peak sites
> gene_methyintensity <- gene_methy_intensity(peak_inform=DM_sites_infor,txdbinfor=TXDB,GENE_ANNO_GTF=NA, species="human")
### Obtain gene expression for INPUT samples
> get_gene_express <- Get_express_data(INPUT_BAM=c(INPUT_BAM,TREATED\_INPUT\_BAM ), 
                                      isPairedEnd=FALSE,species="human",
                                      GENE_ANNO_GTF = gtf)
### Detect the differential expression gene
> obtain_DEgene <- Select_DEgene(gene_count_infor=get_gene_express,
                               cond1="control", 
                               cond2="treated",
                               num_cond1=2, 
                               num_cond2=2,
                               DIFF_GENE_cutoff_FDR=0.05,
                               DE_CUTOFF_TYPE="padj") 
### Select genes with paired differential expression and differential methylation 
> expr_methy_gene <- match_expr_methy(gene_expre_infor=obtain_DEgene[[1]], 
                                     gene_methy_infor=gene_methyintensity,
                                    OUTPUT_DIR=NA)
### Predicate m6A-reg-exp gene by m6Aexpress model in case-control context
> m6A_Express_model(Input_file=expr_methy_gene,
                  CUTOFF_TYPE="FDR", 
                  FDR=0.05)

## Basic mode: pooled samples from one or multiple conditions together and identify m6A regulated expression gene in a specific context.
### Input BAM files.


### We use GTF file in the following example.
> 

### Predict m6A regulated expression gene by m6A-express model
> m6A_reg\_exp\_gene <- m6Aexpress(express_data=INPUT_BAM, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, annot_type="hg19", GENE_ANNO_GTF=gtf, pvalue=0.05,mode="basic")

## DE-DM mode: In this case, we will detect whether the differential m6A methylation peak sites regulated expression that caused the differential expression genes
### Predict the differential expression genes are regulated by differential methylation peak sites by m6A-express model

> IP\_BAM <- c(f1,f2)

> TREATED\_IP\_BAM <- c(f3,f4)

> INPUT\_BAM <- c(f5,f6)

> TREATED\_INPUT\_BAM <- c(f7,f8)

> m6A\_reg\_exp\_gene <- m6Aexpress(express_data=INPUT_BAM, treated_express_data=TREATED_INPUT_BAM, IP_BAM=IP_BAM, TREATED_IP_BAM=TREATED_IP_BAM, INPUT_BAM=INPUT_BAM, TREATED_INPUT_BAM=TREATED_INPUT_BAM,annot_type="hg19", GENE_ANNO_GTF=gtf,pvalue=0.05,mode="DE-DM")
