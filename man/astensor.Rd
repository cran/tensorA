\name{as.tensor}
\alias{as.tensor}
\alias{as.tensor.default}
\alias{as.tensor.tensor}
\title{Coercion to a tensor}
\description{
  Coerces a array to a tensor keeping dimension and names.
}
\usage{
as.tensor(X,...)
\method{as.tensor}{default}(X,...,dims=NULL)
\method{as.tensor}{tensor}(X,...)
}
\arguments{
  \item{X}{a multidimensional array}
  \item{\dots}{further generic arguments}
  \item{dims}{the new dim attribute to be used}
}
\details{
  The main idea is that a multiway array like a vector or a matrix is
  nothing else than a tensor for R, but it still needs the tensor class
  be used with the tensorA library. However this is more a convenience
  function for
  migraters than a proper way construct a tensor, which is done by
  \code{\link{to.tensor}}. 
} \value{
  a tensor containing the same data as X
}
\note{
  You should typically use the \code{\link{to.tensor}} to generate a
  tensor, when you want to write vectorizable functions for tensors.
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- diag(5)
as.tensor(A)
}
\author{K. Gerald van den Boogaart}
\keyword{arith}