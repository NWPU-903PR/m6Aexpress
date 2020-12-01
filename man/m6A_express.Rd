\name{m6Aexpress}
\alias{m6Aexpress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting context-specific m6A regulation of gene expression
}
\description{
The package is designed to predict genes whose expressions are regulaed by m6A methylation based on a negative binomial regression model with pooling prior information across genes to explicitly estimate the sample-to-sample variability between biological replicates and improve the precision and robustness of regression coefficients, which is suitable for MeRIP-seq data in small sample scenario.
}
\usage{
m6Aexpress(# input files in BAM format
            express_data,
            treated_express_data=character(0),
            IP_BAM,
            INPUT_BAM,
            TREATED_IP_BAM=character(0),
            TREATED_INPUT_BAM=character(0),

            # annotation file in GTF format or TXDB format
            annot_type="hg19",
            species="human",
            GENE_ANNO_GTF=NULL,
            isGTFAnnotationFile=FALSE,
            TXDB=NA,

            # specific parameters
            GENOME = NA,
            UCSC_TABLE_NAME = "knownGene",
            isPairedEnd = FALSE,
            model="basic",

            # threshod cutoff
            pvalue=NA,
            FDR=0.05,
            diff_gene_pvalue=NA,
            diff_gene_fdr=0.05,
            diff_peak_pvalue=NA,
            diff_peak_fdr=0.05,

            # output direction
            OUTPUT_DIR= NA
            )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{express_data}{a character vector giving names of input samples containing read mapping results in BAM format. The input files could be the conventional RNA-seq data or the INPUT aligned MeRIP-seq data, which are used to quantify as the gene expression.
}
  \item{treated_express_data}{a character vector giving names of input treated samples containing read mapping results in BAM format. The treated samples could be the conventional RNA-seq data or the INPUT aligned MeRIP-seq data, which are used to quantify as the gene expression under treated condition. These files are only provided in differential expression context.
}
  \item{IP_BAM}{a character vector giving names of IP samples in MeRIP-seq data, which contained read mapping results in BAM format.
}
 \item{INPUT_BAM}{
a character vector giving names of INPUT samples in MeRIP-seq data, which contained read mapping results in BAM format.
}
 \item{TREATED_IP_BAM}{
a character vector giving names of treated IP samples in MeRIP-seq data, which contained read mapping results in BAM format. These files are only provided in differential methylation context.
}
 \item{TREATED_INPUT_BAM}{
a character vector giving names of treated INPUT samples in MeRIP-seq data, which contained read mapping results in BAM format. These files are only provided in differential methylation context.
}
 \item{annot_type}{
a character string specifying an in-built annotation used to quantify each gene's reads count for expression data. It has four possible values including \code{"mm10"}, \code{"mm9"}, \code{"hg38"} and \code{"hg19"}, corresponding to the NCBI RefSeq annotations for genomes `mm10', `mm9', `hg38' and `hg19', respectively. \code{"hg19"} by default.
}
 \item{species}{
a character string specifying the name of species, which could be \code{"human"}, \code{"mouse"} or \code{"yeast"}. This parameter is only appliable when \code{GENE_ANNO_GTF} is \code{NULL}.
}
 \item{GENE_ANNO_GTF}{
A character string giving name of a user-provided annotation file in GTF format. Note that \code{GENE_ANNO_GTF} will override \code{annot_type} if both provided.
}
 \item{isGTFAnnotationFile}{ logical indicating whether the annotation provided via the \code{annot_type} argument is in GTF format. \code{FALSE} by default. This option is only applicable when \code{GENE_ANNO_GTF} is not \code{NULL}.}
 \item{TXDB}{
An optional TxDb object for gene annotation information used in the analysis, default: NA. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
}
 \item{GENOME}{
A string,such as "hg19" or "mm10", which specifies the genome assembly used. If a gene annotation file is provided, this function will not use it directly; otherwise, this function will download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
}
 \item{UCSC_TABLE_NAME}{
A string, which gives the gene annotation used from UCSC, default: "knownGene". Please use function: supportedUCSCtables() to check available tables.
	Some tables may not be available for all genomes, and the "refGene" table does not work correctly due to multiple occuences of the same transcript on the same chromosome.
}
 \item{pvalue}{ a decimal number, which specifies the p-value cut-off in the predicting m6A regulated expression gene of m6A-express model.
}
 \item{FDR}{ a decimal number, which specifies the FDR cut-off in the predicting m6A regulated expression gene of m6A-express model, default: 0.05.
}
 \item{diff_gene_pvalue}{ a decimal number, which specifies the p-value cut-off to identify differential expression gene in \code{DESeq2} R-package.
}
 \item{diff_gene_fdr}{ a decimal number, which specifies the fdr cut-off to identify differential expression gene in \code{DESeq2} R-package, default: 0.05.
}
 \item{diff_peak_pvalue}{ a decimal number, which specifies the p-value cut-off to identify differential methylation peak sites in \code{QNB} R-package.
}
 \item{diff_peak_fdr}{ a decimal number, which specifies the fdr cut-off to identify differential methylation peak sites in \code{QNB} R-package, default: 0.05.
}
 \item{OUTPUT_DIR}{
A string, which specify the output directory, default: OUTPUT_DIR=NA, the output result will save in the current directory. Otherwise, \code{m6A_express} will output the significant m6A regulated expression gene under given cutoff, e.g. \code{FDR<0.05} or \code{pvalue<0.05}.
}
}
\details{
In order to analysis the regulation of m6A methylation on gene expression, we develop the R package m6A-express. The inputs of m6A-express are BAM files of RNA-seq data or INPUT BAM files of MeRIP-seq data, which are used to quantify gene expression, and IP BAM files with paired INPUT BAM files that are used to quantify methylation inensity in m6A-express model. For the origianl gene annotation file, such as GTF or TXDB format file, user could provide it or directly use the default in-built annotation file by given the function parameter \code{species}, \code{annot_type}, \code{GENOME} and \code{UCSC_TABLE_NAME = "knownGene"}.

The \code{m6Aexpress} function is an all-in-one command that performs all the core functions of the m6A-express R-package. The core function of \code{m6A_express} achieved two modes to predict methylation regulated expression genes:

1. Predict m6A regulated expression genes under multiple conditions, such as in mutiple tissues context, and we call it \code{basic} mode.

2. Predict the differential methylation peak sites regulated differential expression genes in specific context and we call it \code{DE-DM} mode.

}
\value{
In the \code{basic} mode, m6A-express will output results in specific directory including:

1. Peak calling result botain by \code{exomePeak} R-package, which are detailedly explained in \code{exomePeak} R-package.

2. The paired reads count and methylation intensity for each gene in tab-delimited text file.

3. The predicting m6A methylation regulated expression (m6A-reg-exp) genes under given cut-off in .XLS (tab-delimited) format file. The .XLS file includes five columns: the name of m6A-reg-exp gene, Beta0, Beta1, pvalue and FDR.
 \item{gene_name}{The name of m6A methylation regulated expression gene in the given annotation file.}
 \item{Beta0}{Beta0 is a gene-specific intercept and models the baseline log gene expression.}
 \item{Beta1}{Beta1 captures the influence of m6A methylation on gene expression.}
 \item{pvalue}{pvalue of the predicted m6A-reg-exp genes.}
 \item{FDR}{FDR of the predicted m6A-reg-exp genes.}

In the \code{DE-DM} mode, m6A-express will output results in specific directory including:

1. Consistent peak sites information between two conditions obtained by \code{exomePeak} R-package, which are detailedly explained in \code{exomePeak} R-package.

2. Differential expression (DE) gene as .XLS file on disk in the specified directory given DE cut-off (e.g. \code{diff_gene_pvalue=0.05} or \code{diff_gene_fdr=0.05}).

3. The differential m6A methylation (DM) peak sites information saved as .XLS file in the specified directory given DM cut-off (e.g. \code{diff_peak_pvalue=0.05} or \code{diff_peak_fdr=0.05}).

4. The paired reads count and methylation intensity for each differential expression gene with differential methylation peaks (DE-DM) in tab-delimited text file.

5. The predicting m6A methylation regulated expression (m6A-reg-exp) genes out from DE-DM genes under given cut-off (e.g. \code{pvalue=0.05} or \code{FDR=0.05}) in .XLS (tab-delimited) format file. The .XLS file includes five columns: the name of m6A-reg-exp gene, Beta0, Beta1, pvalue and FDR.


}

\author{
Teng Zhang <tengzhang156@126.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Basic mode: pooled samples from one or multiple conditions together and identify m6A-reg-exp gene in a specific context.
# Input BAM files
f1 <- system.file("extdata", "IP1.bam", package="m6Aexpress")
f2 <- system.file("extdata", "IP2.bam", package="m6Aexpress")
f3 <- system.file("extdata", "IP3.bam", package="m6Aexpress")
f4 <- system.file("extdata", "IP4.bam", package="m6Aexpress")
f5 <- system.file("extdata", "Input1.bam", package="m6Aexpress")
f6 <- system.file("extdata", "Input2.bam", package="m6Aexpress")
f7 <- system.file("extdata", "Input3.bam", package="m6Aexpress")
f8 <- system.file("extdata", "Input4.bam", package="m6Aexpress")
IP_BAM <- c(f1,f2,f3,f4)
INPUT_BAM <- c(f5,f6,f7,f8)

# Input the annotation file
gtf <- system.file("extdata", "hg19toy.gtf", package="m6Aexpress")

# Predict m6A-reg-exp gene by m6A-express model
m6A_reg_exp_gene <- m6Aexpress(express_data=INPUT_BAM, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, annot_type="hg19", GENE_ANNO_GTF=gtf,isGTFAnnotationFile=TRUE, pvalue=0.05,mode="basic")

# DE-DM mode: Predict the differential methylation peak sites regulated differential expression genes
## Do not to run the following code because of very little numberic DE-DM genes for the input toy BAM files

\dontrun{
# Input BAM files
IP_BAM <- c(f1,f2)
TREATED_IP_BAM <- c(f3,f4)
INPUT_BAM <- c(f5,f6)
TREATED_INPUT_BAM <- c(f7,f8)

# Predict m6A-reg-exp genes out from DE-DM genes
m6A_reg_exp_gene <- m6Aexpress(express_data=INPUT_BAM, treated_express_data=TREATED_INPUT_BAM, IP_BAM=IP_BAM, TREATED_IP_BAM=TREATED_IP_BAM, INPUT_BAM=INPUT_BAM, TREATED_INPUT_BAM=TREATED_INPUT_BAM,annot_type="hg19", GENE_ANNO_GTF=gtf,isGTFAnnotationFile=TRUE, pvalue=0.05,mode="DE-DM")
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
