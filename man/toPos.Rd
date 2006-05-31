\name{toPos.tensor}
\alias{toPos.tensor}
\title{get the position of an index of tensor }
\description{
  Calculates the position of a tensor index, which specified in any
  possible way.
}
\usage{
toPos.tensor(M,l=NULL,mnames=names(dim(M)),by=NULL,...,both=FALSE,missing.ok=FALSE)
}
\arguments{
  \item{M}{a tensor}
  \item{l}{a vector specifying the indices as positions or names}
  \item{mnames}{The names of the indices of the tensor. This can be
    specified instead of M. }
  \item{both}{Matches the index in its covariate and contravariate
    form. }
  \item{by}{the list dimension, all operations are done in parallel for
    all levels of these dimensions. Thus in the case of toPos all other
    dimensions are returned if they are not specified.}
  \item{...}{not used}
  \item{missing.ok}{If TRUE does give an error on missing dimension. Rather
    returns NA in that place.}
}
\value{
  a numeric vector giving the positions of the dimensions selected.
}
\details{
  The function is only here to provide a consistent interface which
  provides the same functionality for positions and characters.
}
\references{
}
\seealso{
}
\examples{
A <- to.tensor(1:30,c(a=2,b=3,c=5))
toPos.tensor(A,c("b","c"))
toPos.tensor(A,c(2,1))     # only returns the values
toPos.tensor(A,c("^a"),both=TRUE)
}
\author{K. Gerald van den Boogaart}
\keyword{}