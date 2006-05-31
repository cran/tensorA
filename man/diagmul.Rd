\name{diagmul.tensor}
\alias{diagmul.tensor}
\title{Multiplication of a tensor with a tensor given by its diagonal}
\description{
  This is a convenience function for scaling elements of a tensor with
  different numbers based on their position in the tensor.
}
\usage{
diagmul.tensor(X,i=names(D),D,j=i,by=NULL)
}
\arguments{
  \item{X}{The tensor to be scaled}
  \item{D}{A tensor containing scaling constants}
  \item{i}{numeric of character vector giving the dimensions of X to be
    used for the product.}
  \item{j}{numeric of character vector giving the dimensions of D to be
    used for the product.}
  \item{by}{Every operation is parallel for all levels of by in X and/or D. }
}
\value{
  A tensor with the shape   and dimensions as X with entries
  \eqn{X_{ik}}{X_ik} scaled by \eqn{D_im}{D_im}, where  \eqn{i}{i} and
  \eqn{k}{k} can represent multi-indices.
}
\details{
  Let \deqn{X_{i_1\ldots i_d k_1 \ldots k_d}} and
  \deqn{D_{j_1\ldots j_d}} than the result is:
  \deqn{E_{i_1\ldots i_d k_1 \ldots k_d}=X_{i_1\ldots i_d k_1 \ldots k_d}D_{j_1\ldots j_d}}
}
\references{
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
(A <- matrix(rep(1:3,each=3),nrow=3))
(b <- to.tensor(c(1,1/2,1/3)))
diagmul.tensor(as.tensor(A),2,as.tensor(c(1,1/2,1/3)),1)
diagmul.tensor(as.tensor(A),1,as.tensor(c(1,1/2,1/3)),1)
A \%*\% diag(b)
diag(b) \%*\% A
}
\author{K. Gerald van den Boogaart}
\keyword{arith}