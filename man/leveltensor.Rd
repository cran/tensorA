\name{level.tensor}
\alias{level.tensor}
\title{The level (number of indices) of a tensor}
\description{
  The level of a tensor is the number of dimensions or subscripts used.
}
\usage{
level.tensor(X,...)
}
\arguments{
  \item{X}{the tensor to be used}
  \item{...}{not used}
}
\value{
  the number of levels
}
\details{
  The level of the tensor is the length of its dim attribute. Objects
  without a dim attribute get level 1 if they are of length > 1 and are
  marked as scalars by 0 level otherwise.
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- to.tensor(1:24,c(a=1,b=2,c=3,d=4))
level.tensor(A)
level.tensor(matrix(1))
level.tensor(1:10)
level.tensor(1)
}
\author{K. Gerald van den Boogaart}
\keyword{arith}