\name{obtain_HVP_sites}
\alias{obtain_HVP_sites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain High Variable Peak sites
}
\description{
This function is to get the high variable peak sites given RNA methylation peak sites under multiple context.
}
\usage{
obtain_HVP_sites(peak_inform,CV_values=0.3,num_sample_subgroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak_inform}{
A list data format, one is dataframe within the consistent peak sites information (e.g. "gene name","seqname","start","end","strand", reads count in each peak of each input samples) under multiple context (e.g. multiple sub-tissues in Brain or Intestine), another is a vector with libray sizes factor to normalize methylation level of each peak site.
}
 \item{CV_values}{
A decimal number, which is the cutoff for coefficient of variation (CV) used to select high variable peaks. The high variable peaks are defined as CV values are bigger than the cutoff: default \code{CV_values=0.3}.
 }
 \item{num_sample_subgroup}{
A vector of integer, which specifies the number of samples for each sub-group. For example, \code{num_sample_subgroup=c(2,2,3,1)} means four sub-groups have two samples, two samples, three samples and one sample in each sub-group, respectively.    
 }
}
\details{
The high variable peak sites means whose methylation level (\eqn{log2} ratio for normalized reads count by size factor of each peak site) of each sub-group with higher coefficient of variation value than given cutoff: \code{CV_values} across sub-groups.
}
\value{
Results will return a list including the information of high variable peak site in dataframe format and a library size factor in vector format.
  \item{high_variable_peak}{A dataframe including the information of high variable methylation peak sites}
  \item{size_factor}{A float vector, which specifices the library size factor for IP samples}
}
\references{
An, Sanqi, et al. "Integrative network analysis identifies cell-specific trans regulators of m6A." \emph{Nucleic acids research} 48.4 (2020): 1715-1729.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{##Get consist peak site
Get_peak_infor <- Get_peak_sites(IP_BAM, INPUT_BAM,GENE_ANNO_GTF=gtf)
##Detect high variable peak sites
HVP_infor <- obtain_HVP_sites(peak_inform=Get_peak_infor,CV_values=0.3,
                               num_sample_subgroup=c(2,2,2,2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
