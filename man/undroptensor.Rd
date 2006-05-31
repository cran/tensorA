\name{undrop.tensor}
\alias{undrop.tensor}
\title{Adds a spurious dimension to a tensor}
\description{
  A dimension of length 1 is added a given position to a tensor
}
\usage{
undrop.tensor(A,name,pos=1)
}
\arguments{
  \item{A}{the tensor}
  \item{name}{the name of the dimension to be added}
  \item{pos}{the position, where to insert the new dimension}
}
\value{
  A tensor with one extra dimension of length 1 with name \code{name} at
  position \code{pos}. 
}
\details{
  The function is a pure convenience function.
}
\references{
}
\seealso{
}
\examples{
A <- to.tensor(1:4,c(a=2,b=2))
undrop.tensor(A,"i")
}
\author{K. Gerald van den Boogaart}
\keyword{arith}