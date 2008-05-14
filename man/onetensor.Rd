\name{one.tensor}
\alias{one.tensor}
\title{Creates a tensor with all entries 1}
\description{
  Creates a tensor with all entries one. 
}
\usage{
one.tensor(d=NULL,dn=NULL)
}
\arguments{
  \item{d}{the dimensions of the new tensor}
  \item{dn}{the dimnames of the new tensor}
}
\value{
  A tensor with dim d and all elements one
}
\details{
  \deqn{E_{i_1\ldots i_n}=1}
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
one.tensor(c(a=3,b=3,c=3))
}
\author{K. Gerald van den Boogaart}
\keyword{arith}