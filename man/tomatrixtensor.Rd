\name{to.matrix.tensor}
\alias{to.matrix.tensor}
\title{The matrix corresponding to a tensor seen as a linear mapping of tensors.}
\description{
  A tensor can be seen as a linear mapping of a tensor to a tensor. This
  function gives the corresponding matrix of the mapping. 
}
\usage{
to.matrix.tensor(X,i,j,by=NULL)
}
\arguments{
  \item{X}{The tensor }
  \item{i}{The image indices of the linear mapping}
  \item{j}{The domain indices of the linear mapping}
  \item{by}{the operation is done in parallel for these dimensions}
}
\value{
  if no \code{by} is given a matrix. Otherwise a tensor of level
  \code{2+length(dim(X))[by]} 
  giving matrices for each specification of the by dimensions.
}
\details{
  A tensor can be seen as a linear mapping of a tensor to a tensor. This
  function computes the corresponding matrix, mapping the entries of the
  domain tensor to the entries of the image tensor.
}
\seealso{
  \code{\link{to.tensor}}, \code{\link{solve.tensor}},
  \code{\link{inv.tensor}}, \code{\link{svd.tensor}}
}
\examples{
A <- reorder.tensor(to.tensor(1:30,c(a=2,b=3,c=5)),c("c","a","b"))

to.matrix.tensor(A,"a",c("b","c"))              # matrix(1:30,nrow=2)

to.matrix.tensor(A,c("a","b"),c("c"))           # matrix(1:30,nrow=6)

to.matrix.tensor(A,c("a","b"),by=c("c")) # structure(1:30,dim=c(6,1,5)))
to.matrix.tensor(A,c("a"),by=c("c"))     # structure(1:30,dim=c(2,3,5)))
}
\author{K. Gerald van den Boogaart}
\keyword{arith}
