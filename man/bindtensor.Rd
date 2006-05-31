\name{bind.tensor}
\alias{bind.tensor}
\title{A cbind/rbind for tensors}
\description{
  Tensors can be put side by side in one dimension if they are of equal
  size in all other dimensions.
}
\usage{
bind.tensor(A,dA=NULL,B,dB=dA)
}
\arguments{
  \item{A}{the first tensor}
  \item{dA}{the dimension of A to be used for binding the tensors}
  \item{B}{the second tensor}
  \item{dB}{the dimension of B to be used for binding the tensors}
}
\value{
  a tensor with the tensors combined to one
}
\details{
}
\references{
}
\note{
  binding does not preserve the sequence of the entries.
}
\seealso{
  \code{\link{base}{cbind}}
}
\examples{
  A <- to.tensor(1:6,c(a=2,b=3))
bind.tensor(A,"a",A)
bind.tensor(A,"b",A)
}
\author{K.Gerald van den Boogaart}
\keyword{arith}