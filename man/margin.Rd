\name{margin.tensor}
\alias{margin.tensor}
\title{Marginalization of tensors}
\description{
  The function removes dimensions from a tensor by summing all entries
  which only differ in these dimensions.
}
\usage{
margin.tensor(X,i=NULL,by=NULL)
}
\arguments{
  \item{X}{the tensor}
  \item{i}{the dimensions to be removed}
  \item{by}{instead of i the dimensions to be kept}

}
\value{
  The tensor with all elements only differing only in the dimensions
  specified added up and only the other dimensions left over.
}
\details{
  This is a tensor multiplication with the \eqn{1_i}{1_i} tensor.
}
\references{
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- diag(1:5)
A
margin.tensor(A,1)

A <- to.tensor(1:30,dim=c(i=3,j=5,k=2))
ftable(A)
margin.tensor(A,"j")

}
\author{K. Gerald van den Boogaart}
\keyword{arith}