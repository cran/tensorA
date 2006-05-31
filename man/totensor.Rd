\name{to.tensor}
\alias{to.tensor}
\alias{to.tensor.default}
\title{Creates a tensor object}
\description{
  Constructs a "tensor". 
  A tensor is the generalization of vectors and matrices to multi-index
  arrays. 
}
\usage{
to.tensor(X,...)
\method{to.tensor}{default}(X,dims=NULL,ndimnames=NULL,what=1,addIndex=FALSE,...)
}
\arguments{
  \item{X}{the numeric data with the entries of the tensor. If the
    object is already a tensor only the subtensors given by the dimension
    \code{what} are converted }
    \item{dims}{These dimensions to be added for the new tensor.
      If the object is to big or \code{addIndex} an extra dimension is
      added. }
    \item{ndimnames}{The new dimnames to be used}
    \item{what}{a numeric or character vector specifying the dimensions
      to be removed.}
    \item{addIndex}{boolean or character, FALSE says  no additional
      dimension, or string to give the name of the dimension}
    \item{\dots}{further arguments to other instances of the generic function}
}
\value{
  a tensor of the specified shape
}
\details{
  This package provides a class \code{"tensor"}
  allowing easy computation regarding tensorial computation in the
  Einstein convention and allows an easier control of the computation
  than aperm and tensor. The package is made to work with things like
  matrices of matrices and linear mapping of matrices to matrices, etc.
  
  A tensor is a multidimensional array, with specific mathematical
  meaning, generalizing vectors and matrices. Tensors can be added,
  subtracted and
  multiplied and used in linear equations. While two matrices A,B are
  commonly only multiplied in two ways \code{A\%*\%B} or \code{B\%*\%A}
  and have some more \code{t(A)\%*\%B},
  \code{B\%*\%t(A)}, \code{sum(A*B)},
  \code{sum(A*t(B))},\code{kronecker(A,B)}
  the tensor calculus brings all of them into a organized system. 
  \cr
  An important aspect for that is  the name of its dimensions. Thus we
  are not bound to work with rows and columns, but can name the
  dimensions to be multiplied. This leads to much more organized
  computation of linear mappings of matrices or datasets of matrices or
  other genuine tensor arithmetic gets involved.
  \cr
  The package provides a full linear algebra support of tensors
  including tensor addition, tensor multiplication, norms, deltatensors,
  binding,
  inversion, normalization,
  Einstein summing
  convention, trace, , dimension renaming, smart display of tensors,
renaming and reshaping, solving
  equation system and giving decompositions  and parallelized data
  processing ,
  
}
\note{This constructor is not called tensor() according to the general
  convention of constructors to avoid conflicts with the tensor
  multiplication routine
  in the \code{tensor} package
}
\seealso{
  \link{tensorA},  \code{\link{level.tensor}},
  \code{\link{diag.tensor}}, \code{\link{norm.tensor}}
  \code{\link{drag.tensor}}, \code{\link{one.tensor}}, 
  \code{\link{mul.tensor}}, \code{\link{\%e\%}}, \code{\link{\%r\%}}, ,
  \code{\link{drag.tensor}}, , \code{\link{trace.tensor}}, 
  \code{\link{solve.tensor}}, \code{\link{svd.tensor}},
  \code{\link{mean.tensor}}
}
\examples{
A <- to.tensor(1:20,c(U=2,V=2,W=5))
B <- to.tensor(1:20,c(U=2,VV=2,WW=5))
A \%e\% B

}
\author{K. Gerald van den Boogaart}
\keyword{algebra}