\name{power.tensor}
\alias{power.tensor}
\title{Compute the power of a symmetric tensor}
\description{
  A tensor can be seen as a linear mapping of a tensor to a tensor. If
  domain and image are the same and the tensor is definite, we can define
  powers. 
}
\usage{
power.tensor(X,i,j,p=0.5,by=NULL)
}
\arguments{
  \item{X}{The tensor to be decomposed}
  \item{i}{The image dimensions of the linear mapping}
  \item{j}{The domain dimensions of the linear mapping}
  \item{p}{the power of the tensor to be computed}
  \item{by}{the operation is done in parallel for these dimensions}
}
\value{
  a tensor
}
\details{
  A tensor can be seen as a linear mapping of a tensor to a tensor. Let
  denote \eqn{R_i} the space of real tensors with dimensions
  \eqn{i_1...i_d}.
  \cr
  To compute a power \code{dim(X)[i]} and \code{dim(X)[j]} need to be
  equal and the tensor symmetric between these dimension. Some exponents
  are only valid with positive definite mappings. None of these
  conditions is checked. 
}
\note{symmetry of the matrix is assumed but not checked.}
\seealso{
  \code{\link{svd.tensor}},
}
\examples{
A <- to.tensor(rnorm(120),c(a=2,b=2,c=5,d=3,e=2))
AAt <- A \%e\% mark(A,"'",c("a","b"))
AAt

power.tensor(AAt,c("a","b"),c("a'","b'"),-1)

inv.tensor(AAt,c("a","b"))

power.tensor(AAt,c("a","b"),c("a'","b'"),2)
mul.tensor(AAt,c("a","b"),AAt,c("a'","b'"))

power.tensor(power.tensor(AAt,c("a","b"),c("a'","b'"),1/pi),
                     c("a","b"),c("a'","b'"),pi)


AAt <- einstein.tensor(A , mark(A,"'",c("a","b")),by="e")

power.tensor(AAt,c("a","b"),c("a'","b'"),-1,by="e")

inv.tensor(AAt,c("a","b"),by="e")

power.tensor(AAt,c("a","b"),c("a'","b'"),2,by="e")
mul.tensor(AAt,c("a","b"),AAt,c("a'","b'"),by="e")

power.tensor(power.tensor(AAt,c("a","b"),c("a'","b'"),1/pi,by="e"),
c("a","b"),c("a'","b'"),pi,by="e")


}
\author{K. Gerald van den Boogaart}
\keyword{arith}
