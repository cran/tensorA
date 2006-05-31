\name{names.tensor}
\alias{names.tensor}
\alias{names<-.tensor}
\alias{dimnames<-.tensor}
\alias{dimnames.tensor}
\alias{dim<-.tensor}
\title{Getting and setting index and dimensionnames of a tensor}
\description{
  The names of a tensor are the names of its dimension 
}
\usage{
\method{names}{tensor}(x)
\method{names}{tensor}(x) <- value
\method{dimnames}{tensor}(x) 
\method{dimnames}{tensor}(x) <- value
\method{dim}{tensor}(x) <- value
}
\arguments{
  \item{x}{a tensor object}
  \item{value}{The new value. If this is a named list it replaces the
  names of the dimensions. If its an unnamed list it gets the names of
  the dimensions.}
}
\value{
  the names of the dimensions the tensor
}
\details{
  The names of the dimensions of the tensor are very relevant in any
  tensor arithmetic since they are the principle way to specify the
  dimensions to be involved in an operation. The dimnames function is
  here only for convenice to guarantee that the names of the dimnames are
  always the same as the names of the dimensions and to ensure that
  always at least a list with the right length and names. 
}
\references{
  
}
\seealso{
  \code{\link{mul.tensor}} 
}
\examples{
A <- to.tensor(1:20,c(U=2,V=2,W=5))
A
dim(A)
names(A)
names(A) <- c("A","B","C")
A
dim(A)
names(A)

}
\author{K. Gerald van den Boogaart}
\keyword{algebra}