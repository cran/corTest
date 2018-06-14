\name{st5}
\alias{st5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
install
}
\description{
#  Compute p-value for the  equal correlation test with combination of Spearman corrleaion and percentage bend correlation based on a logistic regression model corresponding to two independent groups
}
\usage{
st5(x1,z1,x0,z0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x1}{
a numeric vector
}
  \item{z1}{
a numeric vector
}
  \item{x0}{
a numeric vector
}
  \item{z0}{
a numeric vector
}
}

\value{
p-value of test
}
\references{
Danyang Yu, Zeyu Zhang, Kimberly Glass, Jessica Su, Dawn L. DeMeo, Kelan Tantisira, Scott T. Weiss, Weiliang Qiu(2018). New Statistical Methods for Constructing Robust Differential Correlation Networks.(expected to be public)

}
\author{
Danyang Yu <danyangyu@hnu.edu.cn>,
Weiliang Qiu <stwxq@channing.harvard.edu>
}

\examples{
x1=ghdist(n=100,g=0.2,h=0.2)
x0=ghdist(n=100,g=0.2,h=0.2)
z1=x1+ghdist(n=100,g=0.2,h=0.2)
z0=x0+ghdist(n=100,g=0.2,h=0.2)
p=st5(x1,z1,x0,z0)
print(p)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
