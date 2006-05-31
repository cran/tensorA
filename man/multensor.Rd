\name{mul.tensor}
\alias{mul.tensor}
\title{Tensor multiplication for the tensor class}
\description{
  Performs a tensor multiplication like tensor(), but with named indices,
  keeping dimnames, and vectorized.
}
\usage{
mul.tensor(X,i=c(),Y,j=i,by=NULL)
}
\arguments{
  \item{X}{a tensor to be multiplied}
  \item{i}{numeric or character vector specifying the dimension to be
    used in the multiplication for X}
  \item{Y}{a tensor to be multiplied}
  \item{j}{numeric or character vector specifying the dimension to be
    used in the multiplication for Y}
  \item{by}{the by dimensions if present and not mentioned in i or j are
    used as sequence dimensions. tensors in these dimensions are
    processed in parallel. So in this dimension the product is neither
    inner nor outer but parallel like \code{a*b}, rather than
    \code{a\%*\%b} or \code{a\%o\%b}. Unmentioned dimensions get an
    outer product. Mentioned dimensions an inner. 
  }
}
\value{
  The tensor product of X and Y with respect to the regarding
  dimensions. 
}
\details{
  Say \deqn{X_{i_1\ldots i_n h_1 \ldots h_l}}
  and \deqn{Y_{j_1\ldots j_n k_1 \ldots k_m}}
  the the result is:
  \deqn{E_{h_1\ldots h_l k_1 \ldots k_m}= \sum_{i_1,\ldots,i_n} X_{i_1\ldots i_n h_1 \ldots h_l}Y_{j_1\ldots j_n k_1 \ldots k_m}}

  This is an full outer product with i,j not given and a full inner product
  product of i=dim(X)
  
}
\seealso{
  \code{\link{to.tensor}},  \code{\link{\%e\%}},
  \code{\link{\%r\%}}, \code{\link{diagmul.tensor}},
  \code{\link{einstein.tensor}}, \code{\link{riemann.tensor}},
  \code{\link{solve.tensor}} 
}
\examples{
A <- to.tensor(1:20,c(A=2,B=2,C=5))
B <- to.tensor(1:20,c(D=2,B=2,E=5))
mul.tensor(A,"A",A,"B")


}
\author{K. Gerald van den Boogaart}
\keyword{arith}