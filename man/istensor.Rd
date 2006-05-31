\name{is.tensor}
\alias{is.tensor}
\title{Checking for being a tensor}
\description{
  Checks whether the object has a tensor attribute.
}
\usage{
is.tensor(X)
}
\arguments{
  \item{X}{the objected to be checked}
}
\value{
  boolean 
}
\details{
  This is a simple convenience function to check for the property of
  being a tensor.
}
\references{
  
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- matrix(1:9,nrow=3)
is.tensor(A) # no
A <- to.tensor(A)
is.tensor(A) # yes
}
\author{K. Gerald van den Boogaart}
\keyword{arith}