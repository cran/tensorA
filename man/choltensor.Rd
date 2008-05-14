\name{chol.tensor}
\alias{chol.tensor}
\title{Cholesky decomposition of a tensor}
\description{
  A tensor can be seen as a linear mapping of a tensor to a tensor. This
  function computes its Cholesky decomposition.
}
\usage{
chol.tensor(X,i,j,...,name="lambda")
}
\arguments{
  \item{X}{The tensor to be decomposed}
  \item{i}{The image dimensions of the linear mapping}
  \item{j}{The coimage dimensions of the linear mapping}
  \item{name}{The name of the eigenspace dimension. This is the
    dimension created by the decompositions, in which the eigenvectors
    are \eqn{e_i}}
  \item{\dots}{for generic use only}
}
\value{
  a tensor
}
\details{
  A tensor can be seen as a linear mapping of a tensor to a tensor. Let
  denote \eqn{R_i} the space of real tensors with dimensions
  \eqn{i_1...i_d}.
  \itemize{
  \item{chol.tensor}{Computes for a  tensor
  \eqn{ a_{i_1 \ldots i_dj_1 \ldots j_d} }{ a_{i_1...i_dj_1...j_d} } representing a positive definit mapping
  form \eqn{R_j}{R_j} to \eqn{R_i}{R_j} with equal dimension structure in \eqn{i}{i}
  and \eqn{j}{j} its "Cholesky" decomposition
  \eqn{L_{i_1 \ldots i_d \lambda{}}}{L_{i_1...i_d lambda}} such that
  \deqn{
    a_{i_1...i_dj_1...j_d}=\sum_{\lambda{}} L_{i_1...i_d \lambda{}}L_{j_1...j_d \lambda{}}
  }
}
}
}
\note{
  A \code{by} argument is not necessary, since both processing
  dimensions have to be given.
}
\seealso{
  \code{\link{to.tensor}}, \code{\link{svd.tensor}}
}
\examples{


A <- to.tensor(rnorm(15),c(a=3,b=5))
AAt <- einstein.tensor(A,mark(A,i="a"))
ch <- chol.tensor(AAt,"a","a'",name="lambda")
#names(ch)[1]<-"lambda"
einstein.tensor(ch,mark(ch,i="a")) # AAt

A <- to.tensor(rnorm(30),c(a=3,b=5,c=2))
AAt <- einstein.tensor(A,mark(A,i="a"),by="c")
ch <- chol.tensor(AAt,"a","a'",name="lambda")
einstein.tensor(ch,mark(ch,i="a"),by="c") #AAt

	     

}
\author{K. Gerald van den Boogaart}
\keyword{arith}