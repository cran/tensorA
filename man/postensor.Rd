\name{pos.tensor}
\alias{pos.tensor}
\title{enumeration of index combinations}
\description{
  This gives all combinations of indices of a tensor with dimension d in
  the order of the numbers in the memory.
}
\usage{
pos.tensor(d)
}
\arguments{
  \item{d}{a dim attribute of a tensor}
}
\value{
  a matrix with the same number of rows as the tensor has entries an the
  same number of columns as the tensor has dimensions. Each row
  represents the index combination of a the corresponding element.
}
\details{
  tensors are stored according to the R-convention that the leftmost
  index varies fastest.
}
\seealso{
  \code{\link{reorder.tensor}}
}
\examples{
   (A <- to.tensor(1:20,dim=c(A=2,B=2,C=5)))
   pos.tensor(dim(A))
}
\author{K.Gerald van den Boogaart}
\keyword{array}