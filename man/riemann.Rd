\name{riemann.tensor}
\alias{riemann.tensor}
\alias{\%r\%}
\alias{\%r\%.tensor}
\title{Tensor multiplication with Riemann's convention}
\description{
  Multiplies tensors by multiplying over all pairs with one covariate
  and one contravariate variable with the same name according
  to Riemann's summing convention.
}
\usage{
riemann.tensor(...,only=NULL,by=NULL)
## Methods for class tensor
# x \%r\% y
## Default method
# x \%r\% y
}
\arguments{
  \item{\dots}{some tensors, or a renaming code}
  \item{only}{an optional list of the dimension names to be recognized
    for duplication to allow parallel processing on lists of tensors}
  \item{x}{a tensor}
  \item{y}{a tensor}
  \item{by}{Riemannian summing is done in parallel in these dimensions.}
}
\value{
  the tensor product of all the tensors along all duplicate dimensions.
}
\details{
  see \code{\link{mul.tensor}} on details on tensor
  multiplication. In \code{einstein.tensor} complex operations can be
  performed by command and renaming code: The arguments are processed
  from left to right and multiplied. Unnamed attributes are regarded as
  tensors or scalars and
  multiplied with the current result by the Riemann summing convention,
  which means an inner product over all pairs of covariate and
  contravariate indices with the same
  name. Named attributes can either have the name \code{diag}, which performs a
  \code{diagmul} according to the same-name convention or be of the form
  \code{A="B"} or \code{"A"="B"}, for which we have two cases. Typically
  both are given covariate. The first specifies the covariate to be used
  in the multiplication and the second the contravariate.
  If both
  names are
  present in the current result, an inner multiplication (trace) of on
  these two dimensions is
  performed. If only the covariate or the contravariate is present up to
  this point, the specific
  dimension is renamed to the second name, but keeps its type. This
  renaming might be
  visible in the result or inducing a multiplication according to the
  Riemann convention later if the other shows up. 
}
\seealso{
  \code{\link{mul.tensor}}, \code{\link{to.tensor}}, \code{\link{riemann.tensor}}
}
\examples{
A <- to.tensor(1:20,c(U=2,"^V"=2,W=5))
B <- to.tensor(1:20,c("^U"=2,V=2,Q=5))
riemann.tensor(A,B)
A \%r\% B 

}
\author{K. Gerald van den Boogaart}
\keyword{arith}