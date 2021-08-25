\name{Get_peakinfor}
\alias{Get_peakinfor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Peak calling by exomePeak
}
\description{
Peak calling and obtain the reads count for each m6A site using \code{exomePeak} R package.
}
\usage{
Get_peakinfor(IP_BAM, 
               INPUT_BAM,
               TREATED_IP_BAM=character(0),
               TREATED_INPUT_BAM=character(0), 
               species="human",
               GENOME = NA, 
               UCSC_TABLE_NAME = "knownGene", 
               GENE_ANNO_GTF=NULL, 
               TXDB=NA, 
               OUTPUT_DIR= NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{IP_BAM}{
a character vector giving names of IP samples in MeRIP-seq data, which specifies a number of IP samples from the untreted/specific condition in alligned bam format.
}
 \item{INPUT_BAM}{
a character vector giving names of INPUT samples in MeRIP-seq data, which specifies a number of INPUT samples from the untreted/specific condition in alligned bam format.
}
 \item{TREATED_IP_BAM}{
a character vector giving names of treated IP samples in MeRIP-seq data, which specifies a number of IP samples from the treated condition in alligned bam format. These files are only provided in differential methylation context,default: character(0).
}
 \item{TREATED_INPUT_BAM}{
a character vector giving names of treated INPUT samples in MeRIP-seq data, which specifies a number of INPUT samples from the treated condition in alligned bam format. These files are only provided in differential methylation context, default: character(0).
}
 \item{species}{
a character string specifying the name of species, which could be \code{"human"} and \code{"mouse"}. This parameter is only appliable when \code{GENE_ANNO_GTF} is \code{NULL}.
}
 \item{GENOME}{
A string,such as "hg19" or "mm10", which specifies the genome assembly used. If a gene annotation file is provided, this function will not use it directly; otherwise, this function will download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
}
 \item{UCSC_TABLE_NAME}{
A string, which gives the gene annotation used from UCSC, default: "knownGene". Please use function: supportedUCSCtables() to check available tables.
	Some tables may not be available for all genomes, and the "refGene" table does not work correctly due to multiple occuences of the same transcript on the same chromosome.
}
 \item{GENE_ANNO_GTF}{
A character string giving name of a user-provided annotation file in GTF format. Note that \code{GENE_ANNO_GTF} will override \code{annot_type} if both provided.
}
 \item{TXDB}{
An optional TxDb object for gene annotation information used in the analysis, default: NA. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
}
 \item{OUTPUT_DIR}{
A string, which specify the output directory, default: OUTPUT_DIR=NA, the output result will save in the current directory. Otherwise, \code{Get_peak_sites} will output specific dircetory.
}

}
\details{
For peak calling purpose, it requires the IP and input control samples:
An IP sample is the aligned BAM file from the immunoprecipitated sample using RNA modification antibodies such as anti-m6A; 
The input control sample is the aligned BAM file from the total RNAseq shotgun sequencing.
}
\value{
By default, \code{Get_peak_sites} will output results both 

1. A folder named "exomePeak_output" by default included peaks and consist peaks BED/XLS files on diskunder the specified directory (default: current working directory) and a GRangesList object under the R environment, including information generated in the peak calling process can be accessed with the "mcols" command.
2. A list object under current R environment, one includes the consistent peak sites information (e.g. "gene name","seqname","start","end","strand", reads count in each peak of each input samples), another included the libray sizes factor to normalize methylation level of each peak site in further analysis.
}
\references{
Meng, Jia, Xiaodong Cui, Manjeet K. Rao, Yidong Chen, and Yufei Huang. "Exome-based analysis for RNA epigenome sequencing data." Bioinformatics 29, no. 12 (2013): 1565-1567.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{exomepeak}}
}
\examples{
\dontrun{
## Peak calling for case-control context
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
# Get peak sites information
Get_peak_infor <- Get_peakinfor(IP_BAM, INPUT_BAM,GENE_ANNO_GTF=gtf)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
