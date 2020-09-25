\name{st6}
\alias{st6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test Differential Correlation Using st6 Method
}
\description{
Compute p-value for the  equal correlation test with another way to combine Spearman corrleaion and percentage bend correlation based on a multiple logistic regression model corresponding to two independent groups

}
\usage{
st6(x1, z1, x0, z0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x1}{
a numeric vector
}
  \item{z1}{
a numeric vector with same length as \code{x1}
}
  \item{x0}{
a numeric vector
}
  \item{z0}{
a numeric vector with same length as \code{x0}
}
}

\value{
A list with 3 elements:
  \item{u}{the vector of Fisher score.}
  \item{pval}{
p-value of test for testing if correlation between \code{x1} and \code{z1} is the same as that between \code{x0} and \code{z0}}
  \item{cov_u}{covariance matrix for \code{u}.}
}
\references{
Danyang Yu, Zeyu Zhang, Kimberly Glass, Jessica Su, Dawn L. DeMeo, Kelan Tantisira, Scott T. Weiss, Weiliang Qiu(corresponding author). New Statistical Methods for Constructing Robust Differential Correlation Networks to characterize the interactions among microRNAs. Scientific Reports 9, Article number: 3499 (2019)
}
\author{
Danyang Yu <dyu33@jhu.edu>,
Weiliang Qiu <weiliang.qiu@gmail.com>
}

\examples{
x1 = ghdist(n = 100, g = 0.2, h = 0.2)
x0 = ghdist(n = 100, g = 0.2, h = 0.2)
z1 = x1 + ghdist(n = 100, g = 0.2, h = 0.2)
z0 = x0 + ghdist(n = 100, g = 0.2, h = 0.2)
res = st6(x1, z1, x0, z0)
print(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ test }% __ONLY ONE__ keyword per line
