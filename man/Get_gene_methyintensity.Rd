\name{gene_methy_intensity}
\alias{gene_methy_intensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get methylation intensity
}
\description{
Calculate the methylation intensity of each gene given specific methylation sites.
}
\usage{
gene_methy_intensity(peak_inform,
                     txdbinfor, 
                     GENE_ANNO_GTF, 
                     species="human")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak_inform}{
A list data format, one is a dataframe within specific peak sites (\code{"differential methylation peak (DM)"} or \code{"high variable peak (HVP)"}) information (e.g. "gene name","seqname","start","end","strand", reads count in each peak of each input samples), another is a vector with libray sizes factor to normalize methylation level of each peak site.
}
  \item{txdbinfor}{
An optional TxDb object for gene annotation information used in the analysis, default: NA. The gene_methy_intensity function will 
	first look at txdbinfor, then GENE_ANNO_GTF for gene annnotation information. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
}
  \item{GENE_ANNO_GTF}{ A string, which specifies a gene annotation GTF file if available, default: NA
}
 \item{species}{
A character string specifying the name of species, which could be \code{"human"}, \code{"mouse"}. This parameter is only appliable when \code{GENE_ANNO_GTF} is \code{NULL}.
}

}
\details{
This function is to obtain the genes' methylation intensity given the specific peak sites (\code{"differential methylation peak (DM)"} sites or \code{"high variable peak (HVP)"} sites). The methylation intensity are defined as the mean methylation level of one gene's specific peak sites. For each specific peak site of one gene, its methylation level is decayed by the distance to stop codon. More far away from stop codon, the methylation signial is more weak. 
}
\value{
\code{gene_methy_intensity} will return a dataframe with the methylation intensity for each gene of each sample.
}

\examples{
\dontrun{
gene_methyintensity <- gene_methy_intensity(peak_inform=get_peak_site,txdbinfor=TXDB,GENE_ANNO_GTF=NA, species="human")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
