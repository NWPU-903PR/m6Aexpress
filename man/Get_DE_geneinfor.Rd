\name{Select_DEgene}
\alias{Select_DEgene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Difference expression analysis
}
\description{
Obtain the differential gene information and sequence library size facor by DESeq2
}
\usage{
Select_DEgene(gene_count_infor,
              cond1, 
              cond2,
              num_cond1, 
              num_cond2,
              DIFF_GENE_CUTOFF_PVALUE=0.05,
              DIFF_GENE_cutoff_FDR=0.05,
              DE_CUTOFF_TYPE="padj")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_count_infor}{
A list data format, which includes a dataframe quantified the gene expression with reads count for each gene in each sample under specific context and a vector of size factor to normalize gene expression of INPUT samples in the furture. 
}
  \item{cond1}{
A string,such as "control", which specifices the name for the one condition.
}
  \item{cond2}{
A integer,such as "treated", which specifices the name for the second condition.
}
  \item{num_cond1}{
A integer,such as "2", which specifices the number of samples in the one condition.
}
  \item{num_cond2}{
A integer,such as "2", which specifices the number of samples in the second condition.
}
  \item{DIFF_GENE_CUTOFF_PVALUE}{
A decimal number, which specifies the p-value cut-off to identify differential expression gene in \code{DESeq2} R-package, default: 0.05.
}
  \item{DIFF_GENE_cutoff_FDR}{
A decimal number, which specifies the fdr cut-off to identify differential expression gene in \code{DESeq2} R-package, default: 0.05.
}
  \item{DE_CUTOFF_TYPE}{
A string, such as "padj", which specifies the tpye of cut-off to identify differential expression gene, default: \code{DE_CUTOFF_TYPE="padj"}.
}
}
\details{
This function is to identify different expresseion gene by DESeq2 and return reads count for different expression (DE) gene, library size factor for further sample normalization and DE gene information, such as \code{"log2foldchange"}, \code{"p-value"},\code{"padj"}.
}
\value{
This function will return a list including reads count for DE gene, size facor and DE gene information.

  \item{DE_gene_count}{A dataframe, which includes the reads count of DE gene for samples under two specific conditions.}
  
  \item{size_factor}{A vector in float type, which are used to normalize samples.}
  
  \item{significant_DE_gene}{A dataframe format, which includes the DE gene information.}
}
\references{
Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. \emph{Genome Biology}, 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
}

\examples{
\dontrun{
##get reads count for each gene
get_gene_express <- Get_express_data(INPUT_BAM=INPUT_BAM, 
                                      isPairedEnd=FALSE,
                                      GENE_ANNO_GTF = gtf)
##do the DE analysis
obtain_DEgene <- Select_DEgene(gene_count_infor=get_gene_express,
                               cond1="control", 
                               cond2="treated",
                               num_cond1=2, 
                               num_cond2=2,
                               DIFF_GENE_CUTOFF_PVALUE=0.05,
                               DIFF_GENE_cutoff_FDR=0.05,
                               DE_CUTOFF_TYPE="padj") 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
