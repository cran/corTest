\name{fisher_transfer_test}
\alias{fisher_transfer_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for equal correlation
}
\description{
Compute p-value with Fisher’s Z-transformation test. If biasCorrection is true, the corrected correlation is used. The formula is rho.corrected = rho - rho/(2*(n-1)).
}
\usage{
fisher_transfer_test(x1, z1, x0, z0, biasCorrection = TRUE)
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
  \item{biasCorrection}{
a boolean value
}
}
\value{
p-value of test for testing if correlation between \code{x1} and \code{z1} is the same as that between \code{x0} and \code{z0}
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
p = fisher_transfer_test(x1, z1, x0, z0)
print(p)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ test }% __ONLY ONE__ keyword per line
