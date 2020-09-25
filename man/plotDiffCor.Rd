\name{plotDiffCor}
\alias{plotDiffCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatter Plot of 2 Genes for Cases and Controls
}
\description{
Scatter plot of 2 genes for cases and controls, superimposed with linear regression lines.
}
\usage{
plotDiffCor(x1, 
	    z1, 
	    x0, 
	    z0, 
	    pval = NULL, 
	    xlab = "gene1", 
	    ylab = "gene2", 
	    title = "scatter plots") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
numeric. vector of gene expression for gene 1 for cases.
}
  \item{z1}{
numeric. vector of gene expression for gene 2 for cases.
}
  \item{x0}{
numeric. vector of gene expression for gene 1 for controls.
}
  \item{z0}{
numeric. vector of gene expression for gene 2 for controls.
}
  \item{pval}{
numeric. p-value for testing differential correlation of the 2 genes between cases and controls.
}
  \item{xlab}{
character. label for x-axis.
}
  \item{ylab}{
character. label for y-axis.
}
  \item{title}{
character. plot title.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with 4 elments:
\item{g}{
A ggplot2 object.
}
\item{dat}{
a data frame with 3 variables: x, z, and grp.
}
\item{coef1}{
a vector of length two giving the intercept and slope of linear regression for cases.
}
\item{coef0}{
a vector of length two giving the intercept and slope of linear regression for controls.
}
}
\references{
Danyang Yu, Zeyu Zhang, Kimberly Glass, Jessica Su, Dawn L. DeMeo, Kelan Tantisira, Scott T. Weiss, Weiliang Qiu(corresponding author). New Statistical Methods for Constructing Robust Differential Correlation Networks to characterize the interactions among microRNAs. Scientific Reports 9, Article number: 3499 (2019)
}
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
\examples{
library(Biobase)
set.seed(1234567)
res = generate_data(n1 = 50, n2 = 60, p1 = 5, p2 = 50)
es = res$es
print(es)


# gene expression data
dat = exprs(es)
print(dim(dat))
print(dat[1:2,1:3])

# 3rd gene
x = dat[3,]
# 5th gene
z = dat[5,]

# for cases
x1 = x[which(es$grp == 1)]
z1 = z[which(es$grp == 1)]

# for controls
x0 = x[which(es$grp == 0)]
z0 = z[which(es$grp == 0)]

# st5
res2 = st5(x1 = x1, z1 = z1, x0 = x0, z0 = z0)
pval = res2$pval

plotDiffCor(x1 = x1, 
	    z1 = z1, 
	    x0 = x0, 
	    z0 = z0, 
	    pval = pval, 
	    xlab = "gene3", 
	    ylab = "gene5", 
	    title = "scatter plots" 
	    )


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
