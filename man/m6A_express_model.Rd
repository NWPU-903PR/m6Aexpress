\name{m6A_Express_model}
\alias{m6A_Express_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
m6Aexpress model
}
\description{
Predicate m6A regulated expression gene by m6Aexpress model in specific context.
}
\usage{
m6A_Express_model(Input_file,
                  CUTOFF_TYPE="FDR",pvalue, 
                  FDR,out_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Input_file}{
A list data format, one is a dataframe including the paired reads count and methylation intensity for each gene in each samples;  another is a vector with libray sizes factor for INPUT samples used to normalized the gene expression for INPUT samples.
}
  \item{CUTOFF_TYPE}{
A string, such as "FDR", which specifies the tpye of cut-off to identify significant m6A regulated express gene (m6-reg-exp gene), default: \code{CUTOFF_TYPE="FDR"}.
}
  \item{pvalue}{
a decimal number, which specifies the p-value cut-off in the predicting m6A regulated expression gene of m6A-express model.
}
 \item{FDR}{ a decimal number, which specifies the FDR cut-off in the predicting m6A regulated expression gene of m6A-express model, default: 0.05.
}
 \item{out_dir}{
A string, which specify the output directory, default: out_dir=NA, the output result will save in the current directory. Otherwise, \code{m6A_express} will output the significant m6A regulated expression gene under given cutoff, e.g. \code{FDR<0.05} or \code{pvalue<0.05}.
}
}
\details{
This function is used to predicate m6A regulated expression gene by m6Aexpress model, which is based on a Bayesian Hierical Negative Binomial Regression model. Parts of model is accomplished by calling python codel in the \code{extradata} folder of these package. 
}
\value{
The predicting m6A methylation regulated expression (m6A-reg-exp) genes under given cut-off in .XLS (tab-delimited) format file. The .XLS file includes five columns: the name of m6A-reg-exp gene, Beta0, Beta1, pvalue and FDR.
 \item{gene_name}{The name of m6A methylation regulated expression gene in the given annotation file.}
 \item{Beta0}{Beta0 is a gene-specific intercept and models the baseline log gene expression.}
 \item{Beta1}{Beta1 captures the influence of m6A methylation on gene expression.}
 \item{pvalue}{pvalue of the predicted m6A-reg-exp genes.}
 \item{FDR}{FDR of the predicted m6A-reg-exp genes.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Teng Zhang <tengzhang126@163.com>
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
##match gene expression and methylation intensity of HVP gene
expr_methy_gene <- match_expr_methy(gene_expre_infor=get_gene_express, 
                                     gene_methy_infor=gene_methyintensity,
                                   OUTPUT_DIR=NA)
##predicate m6A-reg-exp gene
m6A_Express_model(Input_file=expr_methy_gene,
                  CUTOFF_TYPE="FDR", 
                  FDR=0.05)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
