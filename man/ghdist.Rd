\name{ghdist}
\alias{ghdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generating varibles from g-and-h distribution
}
\description{
generating varibles from g-and-h distribution
}
\usage{
ghdist(n, g = 0,h = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the number of the varibles you want to generate
}
  \item{g}{
the parameter g of g-and-h distribution
}
  \item{h}{
the parameter h of g-and-h distribution
}
}

\value{
n varibles generated from g-and-h distribution
}

\examples{
x = ghdist(50, 0.2, 0.2)
print(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")
\keyword{ generation }% __ONLY ONE__ keyword per line
