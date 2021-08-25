\name{match_expr_methy}
\alias{match_expr_methy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Match gene expression and methylation
}
\description{
Match gene expression and methylation intensity
}
\usage{
match_expr_methy(gene_count_infor, 
                 gene_methy_infor,
                 OUTPUT_DIR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_count_infor}{
A list data format, which includes a dataframe quantified the gene expression with reads count for each gene in each sample under specific context and a vector of size factor to normalize gene expression of INPUT samples in the furture. If the input data are the differential expression gene, it included the reads count of each DE gene, normalized size factor for Input sample and the differential expression information. 
}
  \item{gene_methy_infor}{
A dataframe with the methylation intensity for each gene of each sample.
}
  \item{OUTPUT_DIR}{
A string, which specify the output directory, default: OUTPUT_DIR=NA, the output result will save in the current directory. Otherwise, \code{match_expr_methy} will output genes with paired gene expression and methylation intensity in specific directory given by user.

}
}
\details{
This function is to select genes with paired gene expression (reads count) and methylation intensity in specific context.
}
\value{
\code{match_expr_methy} will output results both:
1. As TAB files on disk (default: "expr_methy.tab") under the specified directory (default crrent working directory), which includes genes with paired reads count and methylation intensity.
2. Retured a list data format uder the R environment:
  \item{gene_expr_methy}{A dataframe including the paired reads count and methylation intensity for each gene in each samples}
  \item{size_factor}{A float vector, which specifices the library size factor for INPUT samples, which is used to normalized the gene expression for INPUT samples}
}

\examples{
\dontrun{
##In differential expression and differential methylation context
##Obtain the DE gene
obtain_DEgene <- Select_DEgene(gene_count_infor=get_gene_express,
                               cond1="control", 
                               cond2="treated",
                               num_cond1=2, 
                               num_cond2=2,
                               DIFF_GENE_CUTOFF_PVALUE=0.05,
                               DIFF_GENE_cutoff_FDR=0.05,
                               DE_CUTOFF_TYPE="padj") 
##Obtain the DM peak sites
DM_sites_infor <- DM_detect(peak_inform=Get_peak_infor,DM_CUTOFF_TYPE=pvalue,num_ctl=2,
                            diff_peak_pvalue=0.05)
##Calculate the methylation intensity for genes with DM peak sites
gene_methyintensity <- gene_methy_intensity(peak_inform=DM_sites_infor,
                                            txdbinfor=TXDB,
                                            GENE_ANNO_GTF=NA, 
                                            species="human")
##match gene expression of DE gene and methylation intensity of DM gene
expr_methy_gene <- match_expr_methy(gene_count_infor=obtain_DEgene, 
                                     gene_methy_infor=gene_methyintensity,
                                    OUTPUT_DIR=NA)
##In high variable peak sites context
##Obtain the reads count for each gene
get_gene_express <- Get_express_data(INPUT_BAM=INPUT_BAM, 
                                      species="human",
                                      annot_file="hg19",
                                      isPairedEnd=FALSE,
                                     GENE_ANNO_GTF = NULL,
                                    isGTFAnnotationFile=FALSE)
##Obtain the high variable peak site
HVP_infor <- obtain_HVP_sites(peak_inform=Get_peak_infor,CV_values=0.3,
                               num_sample_subgroup=c(2,2,2,2))
##Calculate the methylation intensity for genes with high variable peak sites (HVP)
gene_methyintensity <- gene_methy_intensity(peak_inform=HVP_infor,
                                            txdbinfor=TXDB,
                                            GENE_ANNO_GTF=NA, 
                                            species="human")
##match gene expression and methylation intensity of HVP gene
expr_methy_gene <- match_expr_methy(gene_count_infor=get_gene_express, 
                                    gene_methy_infor=gene_methyintensity,
                                     OUTPUT_DIR=NA)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
