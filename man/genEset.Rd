\name{genEset}
\alias{genEset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate an ExpressionSet Object
}
\description{
Generate an ExpressionSet object.
}
\usage{
genEset(ex, pDat, fDat = NULL, annotation = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ex}{
a data matrix stores gene expression data. Rows are genes and columns are subjects.
}
  \item{pDat}{
a data frame stores phenotype data. Rows are subjects and columns are variables describing subjects. Row names of \code{pDat} must be idential to the column names of \code{ex}.
}
  \item{fDat}{
a data frame stores feature data. Rows are genes and columns are variables describing genes. Row names of \code{fDat} must be idential to the row names of \code{ex}.
}
  \item{annotation}{
character. indicates R Bioconductor annotation package that provides annotation for genes.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A BioConductor ExpressionSet object.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Danyang Yu <dyu33@jhu.edu>,
Weiliang Qiu <weiliang.qiu@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ method }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
