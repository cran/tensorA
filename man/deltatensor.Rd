\name{delta.tensor}
\alias{delta.tensor}
\title{Creates a Kronecker delta tensor}
\description{
  The delta tensor is the tensor equivalent of the identity. 
}
\usage{
delta.tensor(d,mark="'",dn=NULL,by=NULL)
}
\arguments{
  \item{d}{the row dimensions }
  \item{mark}{a character to be concatenated to the names of the row
    dimensions to  get the column dimension names}
  \item{dn}{dimnames for the result}
  \item{by}{the dimensions which should not be duplicated}
}
\value{
  a tensor with dimension \code{c(d,mark(d,mark))}
}
\details{
  \deqn{E_{i_1\ldots i_n j_1\ldots
      j_n}=\delta_{i_1j_1}\ldots\delta_{i_nj_n}}{
  E_{i_1...i_n j_1...j_n}=\delta_{i_1j_1}...\delta_{i_n j_n}}
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
delta.tensor(c(a=2,b=3))
}
\author{K. Gerald van den Boogaart}
\keyword{arith}