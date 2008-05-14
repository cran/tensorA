\name{ftable.tensor}
\alias{ftable.tensor}
\title{Pretty printing of tensors}
\description{
  Returns the tensor as (flat) ftable, providing a pretty output.
}
\usage{
\method{ftable}{tensor}(x,...)
}
\arguments{
  \item{x}{the tensor}
  \item{\dots}{additional arguments to ftable}
}
\value{
  an ftable containing the same data as the tensor
}
\details{
  This function is called for a pretty output of a tensor, just try it.
}
\seealso{
  \code{\link{ftable}}
}
\examples{
A <- to.tensor(1:20,c(U=2,V=2,W=5))
A
dim(A)
names(A)
dimnames(A)

ftable(to.tensor(A))
ftable(to.tensor(c(A),dim(A)))

}
\author{K. Gerald van den Boogaart}
\keyword{math}