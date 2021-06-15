\name{Get_expression_data}
\alias{Get_expression_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtained the reads count of each gene by Rsubread
}
\description{
This founction is to obtained gene reads count by Rsubread given Gene Anotation File (GTF) or in-built annotation.
}
\usage{
Get_express_data(INPUT_BAM, 
                 TREATED_INPUT_BAM=character(0),
                  species="human",
                  annot_file="hg19",
                  isPairedEnd=FALSE,
                  GENE_ANNO_GTF = NULL,
                  isGTFAnnotationFile=FALSE
                  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{INPUT_BAM}{
a character vector giving names of INPUT samples in MeRIP-seq data, which specifies a path for allignment reads in BAM format.
}
  \item{TREATED_INPUT_BAM}{
a character vector giving names of treated INPUT samples in MeRIP-seq data, which specifies a path for allignment reads in BAM format. These files are only provided in differential methylation context.
}
  \item{annot_file}{
a character string specifying an in-built annotation used to quantify each gene's reads count for expression data. It has four possible values including \code{"mm10"}, \code{"mm9"}, \code{"hg38"} and \code{"hg19"}, corresponding to the NCBI RefSeq annotations for genomes `mm10', `mm9', `hg38' and `hg19', respectively. \code{"hg19"} by default.
}
  \item{species}{
a character string specifying the name of species, which could be \code{"human"}, \code{"mouse"}. This parameter is only appliable when \code{GENE_ANNO_GTF} is \code{NULL}.
}
   \item{GENE_ANNO_GTF}{
A character string giving name of a user-provided annotation file in GTF format. Note that \code{GENE_ANNO_GTF} will override \code{annot_type} if both provided.
}
  \item{isGTFAnnotationFile}{ logical indicating whether the annotation provided via the \code{annot_type} argument is in GTF format. \code{FALSE} by default. This option is only applicable when \code{GENE_ANNO_GTF} is not \code{NULL}.
}
  \item{isPairedEnd}{ logical indicating if counting should be performed on read pairs or reads. FALSE by default. If TRUE, read pairs will be counted instead of individual reads.}
}
\details{
\code{Get_expression_data} is a general-purpose read summarization function that can assign mapped reads from RNA sequencing to genomic features (gene). 

The function takes as input a set of BAM files containing read mapping results.

\code{Get_expression_data} accepts GTF/GFF annotation formats or in-built annotation for human \code{"hg19"}, \code{"hg38"} or mouse muscle \code{"mm9"} or \code{"mm10"}.

}
\value{
\code{Get_expression_data} will return a list data format, which includes a dataframe with reads count for each gene in each INPUT sample and a vector of size factor for each INPUT sample whose columans.

 \item{gene_reads_count}{A dataframe format, which includes the reads count for each gene in each INPUT sample.}
  \item{size_factor}{A vector in float type, which are used to normalize samples.}
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi (2019).
The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads.
\emph{Nucleic Acids Research}, 47(8):e47.
\url{http://www.ncbi.nlm.nih.gov/pubmed/30783653}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{featureCounts}}
}
\examples{
\dontrun{
##Get gene expression data for each INPUT BAM file.
f1 <- system.file("extdata", "Input1.bam", package="m6Aexpress")
f2 <- system.file("extdata", "Input2.bam", package="m6Aexpress")
f3 <- system.file("extdata", "Input3.bam", package="m6Aexpress")
f4 <- system.file("extdata", "Input4.bam", package="m6Aexpress")
INPUT_BAM <- c(f1,f2,f3,f4)
# Input the annotation file
gtf <- system.file("extdata", "hg19toy.gtf", package="m6Aexpress")
# Get reads count for each gene of each sample by given GTF file.
get_gene_express <- Get_express_data(INPUT_BAM=INPUT_BAM, 
                                      isPairedEnd=FALSE,
                                      GENE_ANNO_GTF = gtf)
## Not run:                                      
# Get reads count for each gene of each sample by in-built annotation file.
# get_gene_express <- Get_express_data(INPUT_BAM=INPUT_BAM, 
                                      species="human",
                                      annot_file="hg19",
                                      isPairedEnd=FALSE,
                                     GENE_ANNO_GTF = NULL,
                                    isGTFAnnotationFile=FALSE) 
}
      
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line