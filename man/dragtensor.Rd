\name{drag.tensor}
\alias{drag.tensor}
\alias{is.covariate}
\alias{is.covariate.tensor}
\alias{is.covariate.numeric}
\alias{is.covariate.character}
\alias{is.contravariate}
\alias{is.contravariate.tensor}
\alias{is.contravariate.numeric}
\alias{is.contravariate.character}
\alias{contraname}
\alias{as.covariate}
\alias{as.contravariate}
\alias{as.covariate.character}
\alias{as.contravariate.character}
\title{Managing covariate and contravariate indices}
\description{
  Each index of a tensor can be covariate or contravariate. The \code{is.*}
  routines check the state of the individual indices based on the
  tensor, its dimension or its index names. \code{drag.tensor} can
  change the state for the tensor and contraname for the names of the tensor.
}
\usage{
drag.tensor(x,g,d)
contraname(x)
is.covariate(x,...)
\method{is.covariate}{tensor}(x,...)
\method{is.covariate}{numeric}(x,...)
\method{is.covariate}{character}(x,...)
as.covariate(x,...)
\method{as.covariate}{character}(x,...)
is.contravariate(x,...)
\method{is.contravariate}{numeric}(x,...)
\method{is.contravariate}{character}(x,...)
as.contravariate(x,...)
\method{as.contravariate}{character}(x,...)
}
\arguments{
  \item{x}{the tensor, its dimension (for \code{*.numeric}) or its
    index-names (for \code{*.character} and \code{contraname})
  }
  \item{g}{The geometry tensor \eqn{g_{ij}}{g_ij} giving the
    transformation between covariate and contravariate. It needs to have
    either  covariate and or contravariate indices.
  }
  \item{d}{a vector (or list) of indices that should be dragged,
    i.e. multiplied with \eqn{g_i^{\;j}}{g_i^j} in the right way such
    that it changes from covariate to contravariate or vice versa. The
    name of the index is kept, only its state changes. The index is thus
    dragged from one state to the other. Indices can given in covariate
    or contravariate form.
  }
  \item{\dots}{only for generic use}
}
\value{
  \item{drag.tensor}{returns a tensor like x but with the dimension }
  \item{is.covariate}{returns a boolean vector giving true for every
    covariate index}
  \item{is.contravariate}{returns a boolean vector giving true for every
    contravariate index}
  \item{as.*}{changes the state of the indices}
  \item{contraname}{returns the names with opposite the opposite covariate
    and contravariate state}
}
\details{
  The covariate and contravariate state of a dimension corresponds to
  column and row  vectors. The transformation between these type is done
  by a linear mapping give by the geometry tensor g, which is the
  identity matrix if the enclosing the geometry is represented by the
  orthonormal basis and ordinary scalar product.  
}
\references{
}
\seealso{
  \code{\link{riemann.tensor}}, \code{\link{to.tensor}}, \code{\link{Tensor}}
}
\examples{
g <- to.tensor(c(1,2,0,1),c(i=2,j=2))
A <- to.tensor(rnorm(8),c(a=2,b=2,c=2))
A2 <- drag.tensor(A,g,c("b","c"))
A2
names(A2)
as.covariate(names(A2))
as.contravariate(names(A2))
is.covariate(A2)
is.contravariate(A2)
riemann.tensor(A2,g)

}
\author{K. Gerald van den Boogaart}
\keyword{arith}
