# Installation Instructions
The m6A-express package is supported by R 3.5.3 or newer versions. First, you need to install the exomePeak package for m6A peak calling:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools","GenomicAlignments","GenomicRanges",
                       "GenomicFeatures","rtracklayer","DESeq2","apeglm","RMariaDB"))
                       
install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/exomePeak_2.16.0.tar.gz", repos = NULL, type="source")
```

Then, please install the QNB package for identifying differential m6A peaks:
```r
install.packages("locfit")
install.packages("https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz", repos = NULL, type="source")
```
Installed the reticulate pacakge to call python code in R
```r
install.packages("reticulate")
##install miniconda to install specific python package
library(reticulate)
install_miniconda()
##install specific python package in R
py_install("statsmodels"); py_install("pandas"); py_install("scipy"); py_install("numpy")
```
Before installing the m6Aexpress package, you should install the following R package:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
##If the R version R>=4.0, the DESeq package should be installed by
##for Linux installation
install.packages("https://www.bioconductor.org/packages/3.11/bioc/src/contrib/DESeq_1.39.0.tar.gz", repos = NULL, type="source")
##for Windows installation
install.packages("https://www.bioconductor.org/packages/3.11/bioc/bin/windows/contrib/4.0/DESeq_1.39.0.zip", repos = NULL, type="source")
##If the R version R>=3.5.3, the DESeq package should be installed by
BiocManager::install('DESeq')
##Install some annotation packages
BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'TxDb.Hsapiens.UCSC.hg19.knownGene', 
                       'TxDb.Mmusculus.UCSC.mm10.knownGene','AnnotationDbi'))
```                       
Now, the m6A-express package can be installed by the following R commands:
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("NWPU-903PR/m6Aexpress")
library(m6Aexpress)
```
# Usage Example
*m<sup>6</sup>A-express* considers a scenario where transcriptome-wide m<sup>6</sup>A under different conditions (treated/disease vs control; different tissues/infection 
stages) are profiled by MeRIP-seq. Note that m6A-express is not restricted to MeRIP-seq but can be applied to any high-throughput methods such as MAZTER-seq or nanopore 
sequencing that quantify m6A stoichiometry. *m<sup>6</sup>A-express* assumes that for an m<sup>6</sup>A site that regulates mRNA expression, the change of its m<sup>6</sup>A 
level is predictive of the change in the expression level of the methylated gene, where the m<sup>6</sup>A level is quantified by MeRIP-seq IP reads and the expression is 
measured by MeRIP-seq Input reads. *m<sup>6</sup>A-express* is an algorithm designed to assess the degree to which such a predictive relationship exists between 
m<sup>6</sup>A levels and gene expressions for the specific conditions under consideration. 

Before applying *m<sup>6</sup>A-express*, m<sup>6</sup>A peaks are first identified from each MeRIP-seq sample using exomePeak (Figure 1 of the paper). The m6A intensity for 
each gene that harbors m<sup>6</sup>A peaks is computed (**Peak Calling and Quantifying Subsection**). m<sup>6</sup>A-express then selects candidate genes based on the 
following criteria: when two conditions (treated vs. control) are considered, candidate genes are differential expression genes that harbor differential m<sup>6</sup>A peaks 
(or DE-DM genes); otherwise, when there are more than two conditions (multiple tissue types or time points), candidate genes are those that contain highly variable m6A peaks 
(HVPs). Afterward, *m<sup>6</sup>A-express* is applied to all the candidate genes. The candidate genes test significant for FDR<0.05 by the Wald test are termed 
m<sup>6</sup>A-reg-exp genes, whose m<sup>6</sup>A intensities are predicted to regulate their gene expressions. Among the outputs of *m<sup>6</sup>A-express* are a list of 
m<sup>6</sup>A-reg-exp genes, their associated regulatory mode and strength (β<sub>1</sub>), the methylation intensities, and the gene expression levels.  

The following codes show how to use the *m<sup>6</sup>A-express* package (named m6Aexpress) to obtain m<sup>6</sup>A-reg-exp genes genes for two or more conditions. Both step-by-step analysis and one-step prediction are detailed. 


## Step by Step Analysis
### For two (case-control) conditions
#### *Peak calling for methylation sites in DE-DM context and obtain consisten peak sites*
```r
f1 <- system.file("extdata", "IP1.bam", package="m6Aexpress")
f2 <- system.file("extdata", "IP2.bam", package="m6Aexpress")
f3 <- system.file("extdata", "IP3.bam", package="m6Aexpress")
f4 <- system.file("extdata", "IP4.bam", package="m6Aexpress")
f5 <- system.file("extdata", "Input1.bam", package="m6Aexpress")
f6 <- system.file("extdata", "Input2.bam", package="m6Aexpress")
f7 <- system.file("extdata", "Input3.bam", package="m6Aexpress")
f8 <- system.file("extdata", "Input4.bam", package="m6Aexpress")
IP_BAM <- c(f1,f2)
INPUT_BAM <- c(f5,f6)
TREATED_IP_BAM <- c(f3,f4)
TREATED_INPUT_BAM <- c(f7,f8)
#Input the gene annotation file  
gtf <- system.file("extdata", "hg19toy.gtf", package="m6Aexpress")
#Obtain the consistent peak sites
Get_peak_infor <- Get_peakinfor(IP_BAM, INPUT_BAM,TREATED_IP_BAM, TREATED_INPUT_BAM, GENE_ANNO_GTF=gtf)
```
#### *Calling differential methylated (DM) peaks among the consistent peaks* 
```r
DM_sites_infor <- DM_detect(peak_inform=Get_peak_infor,DM_CUTOFF_TYPE="pvalue",num_ctl=2, diff_peak_pvalue=0.05)
```
#### *Calculate the methylation intensity for each gene with DM peaks*
```r
gene_methyintensity <- gene_methy_intensity(peak_inform=DM_sites_infor,txdbinfor=NA,GENE_ANNO_GTF=gtf)
```
#### *Obtain their gene expression for INPUT samples*
```r
get_gene_express <- Get_express_data(INPUT_BAM=c(INPUT_BAM,TREATED_INPUT_BAM ), 
                                      isPairedEnd=FALSE,GENE_ANNO_GTF = gtf,isGTFAnnotationFile=TRUE)
```                                      
#### *Identify the differential expression gene*
```r
obtain_DEgene <- Select_DEgene(gene_count_infor=get_gene_express,
                               cond1="control", 
                               cond2="treated",
                               num_cond1=2, 
                               num_cond2=2,
                               DIFF_GENE_cutoff_FDR=0.05,
                               DE_CUTOFF_TYPE="padj") 
```
#### *Select genes with both differential expression and differential methylation*
```r
expr_methy_gene <- match_expr_methy(gene_count_infor=obtain_DEgene, 
                                     gene_methy_infor=gene_methyintensity,
                                    OUTPUT_DIR=NA)
```                                    
#### *Predict m6A-reg-exp genes* 
```r
m6Aexpress_result <- m6A_Express_model(Input_file=expr_methy_gene,
                                           CUTOFF_TYPE="FDR", 
                                            FDR=0.05)
```                                            
#### *Add differential expression and differential methylation in the result* 
```r
m6A_express_addLFC_DDM <- add_LFC_DDM(expre_methyre=m6Aexpress_result, 
                                    DE_gene=obtain_DEgene, methy_intensity=gene_methyintensity,
                                    num_cond1=2, OUTPUT_DIR=NA)
```
### For more than two conditions
#### *Peak calling*
```r
IP_BAM <- c(f1,f2,f3,f4)
INPUT_BAM <- c(f5,f6,f7,f8)
Get_peak_infor <- Get_peakinfor(IP_BAM, INPUT_BAM, GENE_ANNO_GTF=gtf)
```
#### *Detect highly variable peaks  across multiple conditions*
```r
HVP_infor <- obtain_HVP_sites(peak_inform=Get_peak_infor,CV_values=0.3,
                               num_sample_subgroup=c(2,2))
```                            
#### *Calculate the methylation intensity for each gene with highly variable peaks*
```r
gene_methyintensity <- gene_methy_intensity(peak_inform=HVP_infor,GENE_ANNO_GTF=gtf)
```
#### *Obtain gene expressions from INPUT samples in each condition*
```r
get_gene_express <- Get_express_data(INPUT_BAM=c(INPUT_BAM), 
                                      isPairedEnd=FALSE,
                                      GENE_ANNO_GTF = gtf)
```                                      
#### *Select genes with expression and methylation intensity of highly variable peaks*
```r
expr_methy_gene <- match_expr_methy(gene_count_infor=get_gene_express, 
                                     gene_methy_infor=gene_methyintensity,
                                    OUTPUT_DIR=NA)
```                                    
#### *Predicate m6A-reg-exp gene by m6Aexpress model in tissue-specific context* 
```r
m6A_Express_model(Input_file=expr_methy_gene,
                     CUTOFF_TYPE="FDR", 
                      FDR=0.05)
```
## On step to predicate m6A-reg-exp gene
### For two conditions:
```r
IP_BAM <- c(f1,f2)
TREATED_IP_BAM <- c(f3,f4) 
INPUT_BAM <- c(f5,f6) 
TREATED_INPUT_BAM <- c(f7,f8) 
m6A_reg_exp_gene <- m6Aexpress(express_data=INPUT_BAM, treated_express_data=TREATED_INPUT_BAM, 
                                    IP_BAM=IP_BAM, TREATED_IP_BAM=TREATED_IP_BAM, INPUT_BAM=INPUT_BAM, 
                                    TREATED_INPUT_BAM=TREATED_INPUT_BAM,annot_type="hg19", GENE_ANNO_GTF=gtf,
                                    isGTFAnnotationFile=TRUE, pvalue=0.05,mode="DE-DM")
```                                    
### For more than two conditions
```r
IP_BAM <- c(f1,f2,f3,f4)
INPUT_BAM <- c(f5,f6,f7,f8) 
m6A_reg_exp_gene <- m6Aexpress(express_data=INPUT_BAM, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, 
                                annot_type="hg19", GENE_ANNO_GTF=gtf,isGTFAnnotationFile=TRUE, 
                                pvalue=0.05,mode="HVP", CV_values = 0.3, num_sample_subgroup=c(2,2))
```                                

