\name{diag.tensor}
\alias{diag.tensor}
\title{Creates a "diagonal" tensor}
\description{
  The diagonal tensor is the tensor equivalent of the diagonal matrix. 
}
\usage{
diag.tensor(X,mark="'",dn=NULL,by=NULL)
}
\arguments{
  \item{X}{ a tensor containing the diagonal entries. }
  \item{mark}{a character to be concatenated to the names of the row
    dimensions to  get the column dimension names}
  \item{dn}{dimnames which are used twice}
  \item{by}{The diagonal tensor is created for each level of the indices
    in \code{by}. 
  }
}
\value{
  a tensor with dimension \code{c(dim(X),mark(dim(X),mark))}
}
\details{
  \deqn{E_{i_1\ldots i_n j_1\ldots
      j_n}=\delta_{i_1j_1}\ldots\delta_{i_nj_n}}{
  E_{i_1...i_n j_1...j_n}=\delta_{i_1j_1}...\delta_{i_n j_n}
  }
}
\references{
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- to.tensor(1:4,c(a=2,b=2))
diag.tensor(A)
diag.tensor(A,by="b")

}
\author{K. Gerald van den Boogaart}
\keyword{arith}