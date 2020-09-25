\name{construct_network}
\alias{construct_network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct Differential Correlation Network
}
\description{
  Construct differential correlation network with expressionSet,st5 is recommand for testing equal correlation.
}
\usage{
construct_network(es,
                  cor_method = "st5",
                  var.grp,
                  pseudo_adjust_cutoff = FALSE,
                  pAdjMethod = 'fdr',
                  cutoff = 0.05,
                  nPseudo = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{es}{
  an ExpressionSet object of microRNA dataset
}
  \item{cor_method}{
 a string represents the method for equal correlation, 'st5' is recommand.
}
  \item{var.grp}{
character.  phenotype variable name indicating case-control status,0 as control, 1 as case.
}
  \item{pseudo_adjust_cutoff}{
  if the value is TRUE, pseudo probes will be used for setting the cutoff of p-value for differential corrlation test. Otherwise, adjusted p-value will be compared with \code{cutoff}.
}
  \item{pAdjMethod}{
if \code{pAdjMethod='none'}, the function will not do mutiple testing adjustment. If pAdjMethod="fdr"/"BH"/"BY"/"holm"/"hochberg"/"hommel"/"bonferroni"/"BH"/"BY", the specific method will be used for adjusting p-value. \code{pAdjMethod} will only be used when \code{pseudo_adjust_cutoff=FALSE}
}
  \item{cutoff}{
  if p value is smaller than the cutoff, there will be an edge between the two nodes. \code{cutoff} will only be used when \code{pseudo_adjust_cutoff=FALSE}.
}

  \item{nPseudo}{
  if \code{pseudo_adjust_cutoff=TRUE}, then \code{nPseudo} genes will be randomly chosen to produce pseudo genes (i.e., non-differentially correlated with other genes between cases and controls). We will use the minimum of the p-values of these pseudo genes as the cutoff of p-values for differential correlation analysis.
  }
}

\value{
A list with 6 elements:
\item{my_graph}{
  obtained network as igraph object
}
\item{my_dat}{
 obtained netork as data frame with 3 columns: edge id, node_id1,node_id2
}
\item{pvalMat}{
raw p-values for testing differential correlation for each pair of genes
}
\item{pAdjMat}{
adjusted p-values for testing differential correlation for each pair of genes
}
\item{pvalPseudo}{
p-values for testing differential correlation between pseudo genes and other genes
}
\item{alpha1}{
cutoff for p-values for testing differential correlation
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
set.seed(1234567)
res = generate_data(n1 = 50, n2 = 60, p1 = 5, p2 = 50)
es = res$es
print(es)

covCtrl = res$covCtrl
covCase = res$covCase

# we expect cov for 1st 5 genes are different between cases and controls
print(round(covCtrl[1:5, 1:5], 2))
print(round(covCase[1:5, 1:5], 2))

# we expect cov for other genes are same between cases and controls
print(round(covCtrl[6:10, 6:10], 2))
print(round(covCase[6:10, 6:10], 2))

res2 = construct_network(es = es,
                  cor_method = "st5",
                  pseudo_adjust_cutoff = FALSE,
                  var.grp = "grp",
                  pAdjMethod = 'fdr',
                  cutoff = 0.05,
                  nPseudo = 25)

print(res2$graph)
print(res2$network_dat)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ test }% __ONLY ONE__ keyword per line
