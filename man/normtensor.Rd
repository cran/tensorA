\name{norm.tensor}
\alias{norm}
\alias{norm.tensor}
\alias{opnorm}
\alias{opnorm.tensor}
\title{Calculate the Euclidean norm or Euclidean operator norm of a tensor or its subtensors}
\description{
  Calculates the Euclidean norm of a tensor or its subtensors. 
}
\usage{
norm(X,...)
\method{norm}{tensor}(X,i=NULL,...,by=NULL)
opnorm(X,...)
\method{opnorm}{tensor}(X,i=NULL,...,by=NULL)
}
\arguments{
  \item{X}{The tensor}
  \item{i}{For norm the dimensions to of the subtensors to be used. If missing
    the norm of the whole tensor is computed. For opnorm the dimensions
    of the image.}
  \item{\dots}{unused}
  \item{by}{the list dimension, if i is not specified the norm is
    calculated for each of these
    in parallel.}
}
\value{
  
  \item{norm}{either a single number giving the norm of the tensor or a tensors with
  the dimensions i removed
  containing the individual norms in each entry.}
\item{opnorm}{a tensor of dimension \code{dim(X)[by]} giving the
  Euclidean operator norm of the tensor (i.e. its largest singular
  value)
}
}
\details{
  
  \describe{
  \item{norm}{
    The function computes the Euclidean norm, which is the square root
    over the sum of all entries and not the operator norm. }
  \item{opnorm}{
    The function computes the Euclidean operator norm, which is largest
    factor in changing the Euclidean norm, when mapped with the linear
    mapping corresponding to the tensor. }
  }
}

\seealso{
  \code{\link{to.tensor}}
}
\examples{
C <- to.tensor(1:20,c(A=4,B=5))
norm(C,"A")
norm(C,2)
norm(C,c("A","B"))
opnorm(C,"A")
}
\author{K. Gerald van den Boogaart}
\keyword{algebra}
