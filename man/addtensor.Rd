\name{add.tensor}
\alias{add.tensor}
\alias{-.tensor}
\alias{+.tensor}
\alias{*.tensor}
\alias{/.tensor}
\title{Element-wise arithmetic operations +,-,*,/ with tensors}
\description{
  Adds/subs/multiplies/devides tensors element by element . 
  The luxury difference to a simple \code{+} is that we
  do not need to consider the correct permutation of indices  or rules
  on implicit replication, since all of this is handled automatically.   
}
\usage{
add.tensor(X,Y,op="+",only=NULL)
## Methods for class tensor
# x + y
# x - y
# x * y
# x - y
}
\arguments{
  \item{X}{a tensor}
  \item{Y}{a tensor}
  \item{op}{a binary function used to perform the "addition"}
  \item{only}{a list of dimnames that may be considered as equal. This
    parameter is here to allow parallelization of tensors with only
    partially known structure.}
}
\value{
  A tensor giving the element-wise operation X,Y. If some of the
  indices are missing in one of the tensors they are added by repetition.  
}
\details{
  The tensors are properly reordered such that dimensions of the same
  name are identified. If dimensions are missing in one of the tensors
  it is correspondingly repeated. 
}
\references{
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- to.tensor(1:20,c(U=2,V=2,W=5))
add.tensor(A,A)/2 -A
(A+A)/2
A/A
A * 1/A
norm.tensor(reorder.tensor(A,c(2,3,1)) - A)
}
\author{K. Gerald van den Boogaart}
\keyword{arith}