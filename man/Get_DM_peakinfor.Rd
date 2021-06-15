\name{DM_detect}
\alias{DM_detect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect differential RNA methylation sites.
}
\description{
This function is to detect differential RNA methylation sites between two conditions by \code{qnbtest} from \code{QNB} R package.
}
\usage{
DM_detect(peak_inform,
          DM_CUTOFF_TYPE=pavlue,
          num_ctl,
          diff_peak_pvalue,
          diff_peak_fdr
          )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak_inform}{
A list data format, one is a dataframe within the consistent peak sites information (e.g. "gene name","seqname","start","end","strand", reads count in each peak of each input samples), another is a vector with libray sizes factor to normalize methylation level of each peak site in further analysis.
}
  \item{DM_CUTOFF_TYPE}{
A string, such as "pvalue", which specifies the tpye of cut-off to identify differential expression gene, default: \code{DE_CUTOFF_TYPE="pvalue"}.
}
  \item{num_ctl}{
A string,such as "2", which specifices the number of samples in control condition.
}
  \item{diff_peak_pvalue}{
a decimal number, which specifies the p-value cut-off to identify differential methylation peak sites, default: \code{diff_peak_pvalue=0.05}.
}
  \item{diff_peak_fdr}{
a decimal number, which specifies the fdr cut-off to identify differential methylation peak sites, e.g. \code{diff_peak_fdr=0.05}.
}

}
\value{
Results will return a list including the information differential methylation site in dataframe format and a library size factor in vector format.
  \item{DM_peak_infor}{A dataframe including the information of differential methylation sites}
  \item{size_factor}{A float vector, which specifices the library size factor for IP samples}
}
\references{
Liu, Lian, et al. "QNB: differential RNA methylation analysis for count-based small-sample sequencing data with a quad-negative binomial model." \emph{BMC bioinformatics} 18.1 (2017): 1-12.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{qnbtest}}
}
\examples{
\dontrun{
##Get consist peak site
Get_peak_infor <- Get_peak_sites(IP_BAM, INPUT_BAM,GENE_ANNO_GTF=gtf)
##Detect differential peak sites
DM_sites_infor <- DM_detect(peak_inform=Get_peak_infor,DM_CUTOFF_TYPE=pvalue,num_ctl=2,
                            diff_peak_pvalue=0.05)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
