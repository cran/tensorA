\name{mark.tensor}
\alias{mark}
\alias{mark.character}
\alias{mark.numeric}
\alias{mark.tensor}
\title{Marks the names of a tensor with a mark}
\description{
  This modifies the names of the dimensions in a simple and reversible
  way by adding a mark.
}
\usage{
mark(X,mark,...)
\method{mark}{tensor}(X,mark="'",i=1:level.tensor(X),...,by=NULL)
\method{mark}{numeric}(X,mark="'",i=1:length(X),...,by=NULL)
\method{mark}{character}(X,mark="'",i=1:length(X),...,by=NULL)
}
\arguments{
  \item{X}{A tensor or dimension to be marked }
  \item{mark}{a character giving the mark}
  \item{i}{the dimensions to be marked}
  \item{\dots}{generic arguments}
  \item{by}{Dimensions not to be marked. Wins in case of conflicts.}
}
\value{
  A object similar to X but with marked dimensions.
}
\details{
  The concept is very important in tensor algebra since it allows to
  keep dimensions connected without but still
  distinguishable. Eventually later a function for the Riemann summing
  rule will make use of marks to distinguish covariate and contravariate
  dimensions.
}
\seealso{
   \code{\link{delta.tensor}}, \code{\link{diag.tensor}} 
}
\examples{
# The outer product
A <- to.tensor(1:4,c(a=2,b=2))
A %e% mark(A)
}
\author{K. Gerald van den Boogaart}
\keyword{arith}