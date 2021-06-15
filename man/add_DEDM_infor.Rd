\name{add_LFC_DDM}
\alias{add_LFC_DDM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add DE ad DM information
}
\description{
Add the log2foldchange and the degree of differential methylation intensity for m6A-reg-exp gene
}
\usage{
add_LFC_DDM(expre_methyre, 
            DE_gene, 
            methy_intensity,
            num_cond1,OUTPUT_DIR=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expre_methyre}{
A dataframe format, which specifices a set of significant m6A regulated expression genes predicated by m6Aexpress model.
}
  \item{DE_gene}{
A list including reads count for DE gene, size facor and significant DE gene information. In this function, the log2foldchange of significant DE gene will be used.
}
  \item{methy_intensity}{
A dataframe with the methylation intensity for each gene of each sample under case-control context.
}
  \item{num_cond1}{
A integer, such as "2", which specifices the number of samples in the control condition. It will be used to select the methlation intensity samples in control condition.
}
 \item{OUTPUT_DIR}{
A string, which specify the output directory, default: OUTPUT_DIR=NA, the output result will save in the current directory. Otherwise, \code{add_LFC_DDM} will output the specific directory in \code{.xls} format.
}
}
\details{
This function is only used in the differential expression (DE) and differential methylation (DM) context, which will add the log2foldchange of DE gene and the degree of DM to the m6A-reg-exp gene. The degree of DM is defined as the difference of mean methylation intensity between treated and control group.
}
\value{
A dataframe format, which includes the m6A-reg-exp gene information with addded log2foldchagne and the degree of DM.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Teng Zhang <tengzhang156@126.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
m6A_express_addLFC_DDM <- add_LFC_DDM(expre_methyre=m6Areg_expr_gene, 
                                    DE_gene=DE_gene, methy_distdecay=DM_methy,
                                    num_cond1=2, OUTPUT_DIR=NA)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
