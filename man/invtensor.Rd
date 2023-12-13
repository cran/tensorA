\name{inv.tensor}
\alias{inv.tensor}
\title{Inversion of a tensor as linear mapping from tensors to tensors}
\description{
  A tensor can be seen as a linear mapping of a tensor to a tensor. This
  function computes its (generalized-Moore-Penrose) inverse.
}
\usage{
inv.tensor(X,i,...,allowSingular=FALSE,eps=1E-10,by=NULL)
}
\arguments{
  \item{X}{The tensor to be decomposed}
  \item{i}{The image dimensions of the linear mapping}
%  \item{j}{The coimage dimensions of the linear mapping}
  \item{allowSingular}{A boolean, indicating that a
    Moore-Penrose-Inverse should be computed rather than an error
    generated in case of a numerically singular mapping.}
  \item{\dots}{further arguments for generic use}
  \item{eps}{The limit for condition-number, to select an generalized inverse.}
  \item{by}{the operation is done in parallel for these dimensions}
}
\value{
  a tensor containing the inverse mapping. If allowSingular is given and
  the condition number of the matrix is bellow eps a generalized inverse
  is returned. 
}
\details{
  A tensor can be seen as a linear mapping of a tensor to a tensor. 
  \describe{
    \item{inv.tensor}{ Computes the inverse of the mapping}
  }
}
\seealso{
  \code{\link{to.tensor}}, \code{\link{solve.tensor}}, \code{\link{svd.tensor}}
}
\examples{
# SVD
# inv.tensor
R1  <- matrix(rnorm(9),nrow=3)
R1i <- solve(R1)
R2 <- to.tensor(R1,c(a=3,b=3),what=1:2)
R2i <- to.tensor(R1i,c(b=3,a=3),what=1:2)

inv.tensor(R2,"a","b") - R2i
inv.tensor(R2,"a","b",allowSingular=TRUE) - R2i

inv.tensor(rep(R2,4,1,"K"),"a","b",by="K") - rep(R2i,4,1,"K")
inv.tensor(rep(R2,4,1,"K"),"a","b",by="K",allowSingular=TRUE) - rep(R2i,4,3,"K")


}
\author{K. Gerald van den Boogaart}
\keyword{arith}
