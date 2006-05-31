\name{svd.tensor}
\alias{svd.tensor}
\title{Singular value decomposition of tensors}
\description{
  A tensor can be seen as a linear mapping of a tensor to a tensor. This
  function computes the singular value decomposition of this mapping
}
\usage{
svd.tensor(X,i,j=NULL,...,name="lambda",by=NULL)
}
\arguments{
  \item{X}{The tensor to be decomposed}
  \item{i}{The image dimensions of the linear mapping}
  \item{j}{The coimage dimensions of the linear mapping}
  \item{name}{The name of the eigenspace dimension. This is the
    dimension created by the decompositions, in which the eigenvectors
    are \eqn{e_i}}
  \item{\dots}{further arguments for generic use}
  \item{by}{the operation is done in parallel for these dimensions}
}
\value{
  a tensor or in case of svd a list u,d,v, of tensors like in \code{\link[base]{svd}}.
}
\details{
  A tensor can be seen as a linear mapping of a tensor to a tensor. Let
  denote \eqn{R_i} the space of real tensors with dimensions \eqn{i_1...i_d}.
  \item{svd.tensor}{Computes a singular value decomposition
  \eqn{u_{i_1...i_d\lambda{}}},\eqn{d_\lambda{}}, \eqn{v_{j_1...j_l}\lambda{}} such
  that u and v correspond to orthogonal mappings from \eqn{R_\lambda{}} to
  \eqn{R_i} or \eqn{R_j} respectively.}
}
\references{
}
\seealso{
  \code{\link{to.tensor}}, \code{\link{to.matrix.tensor}},
  \code{\link{inv.tensor}}, \code{\link{solve.tensor}}
}
\examples{
# SVD
A <- to.tensor(rnorm(120),c(a=2,b=2,c=5,d=3,e=2))

SVD <- svd.tensor(A,c("a","d"),c("b","c"),by="e")
dim(SVD$v)
# Property of decomposition
einstein.tensor(SVD$v,diag=SVD$d,SVD$u,by="e") # A
# Property of orthogonality
SVD$v \%e\% SVD$v[[lambda=~"lambda'"]]         # 2*delta.tensor(c(lambda=6))
SVD$u \%e\% SVD$u[[lambda=~"lambda'"]]         # 2*delta.tensor(c(lambda=6)))
SVD$u \%e\% mark(SVD$u,"'",c("a","d"))  # 2*delta.tensor(c(a=2,d=3)))



}
\author{K. Gerald van den Boogaart}
\keyword{arith}
