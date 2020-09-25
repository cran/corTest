\name{generate_data}
\alias{generate_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Generate Data for Control Subjects and Case Subjects
}
\description{
The function is to generate expression level matrixes of control subjects and case subjects. X matrix is for control subjects with the default sample size n1=50. Z matrix is for case subjects with the default sample size n2=60.
X is generated from multivariate normal distribution N (0, SigmaX), where SigmaX is a block matrix ((SigmaP1, 0), (0, SigmaP2)), sigmaP1 is the p1*p1 matrix and SigmaP2 is the p2*p2 matrix. Z is generated from multivariate normal distribution N (0, SigmaZ), where SigmaZ is a block matrix ((E_P1, 0), (0, SigmaP2)) and E_P1 is p1*p1 identity matrix.
}
\usage{
generate_data(n1, n2, p1, p2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{n1}{
a numeric value representing the sample size of control subjects
}
  \item{n2}{
a numeric value representing the sample size of case subjects
}
  \item{p1}{
a numeric value representing the dimension of sigma_1
}
  \item{p2}{
a numeric value representing the dimension of sigma_0
}
}
\value{
A list with 3 elements:
\item{es}{An ExpressionSet object contains gene expression data.}
\item{covCase}{
  Covariance matrix of all genes for case subjects.
}
\item{covCtrl}{
  Covariance matrix of all genes for control subjects.
}
}
\references{
Danyang Yu, Zeyu Zhang, Kimberly Glass, Jessica Su, Dawn L. DeMeo, Kelan Tantisira, Scott T. Weiss, Weiliang Qiu(corresponding author). New Statistical Methods for Constructing Robust Differential Correlation Networks to characterize the interactions among microRNAs. Scientific Reports 9, Article number: 3499 (2019)
}
\author{
Danyang Yu <dyu33@jhu.edu>,
Weiliang Qiu <weiliang.qiu@gmail.com>
}

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

# phenotype data
pDat = pData(es)
print(dim(pDat))
print(pDat[1:2,])
print(table(pDat$grp, useNA = "ifany"))

# feature data
fDat = fData(es)
print(dim(fDat))
print(fDat[1:2,])
print(table(fDat$memGenes, useNA = "ifany"))

covCtrl = res$covCtrl
covCase = res$covCase

# we expect cov for 1st 5 genes are different between cases and controls
print(round(covCtrl[1:5, 1:5], 2))
print(round(covCase[1:5, 1:5], 2))

# we expect cov for other genes are same between cases and controls
print(round(covCtrl[6:10, 6:10], 2))
print(round(covCase[6:10, 6:10], 2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ test }% __ONLY ONE__ keyword per line
