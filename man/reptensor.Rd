\name{reptensor}
\alias{rep.tensor}
\title{Repeats a tensor}
\description{
  The tensor is repeated like a number is repeated by rep and an
  additional dimension is added to select the different tensors. 
}
\usage{
\method{rep}{tensor}(x,times,pos=1,name="i",...)
}
\arguments{
  \item{x}{the tensor to be repeated}
  \item{times}{the number of copies that should be created. If \code{times} is
    a vector, x is seen as a sequence of tensors in dimension \code{pos} and
    each of the tensors is repeated according to the corresponding
    entry of times. }
  \item{name}{the name of the additional dimension. if NA no additional
    dimension is used.}
  \item{pos}{the position where the extra dimension should be added}
  \item{\dots}{not used, only here for generic consistency}
}
\value{
  A tensor with one  additional dimensions of length times. 
}
\details{
  This function is modeled as much as possible to mimic rep, by
  repeating tensors rather than numbers. The
  \code{each} argument is not necessary, since sequence of the dimensions
  can more precisely be controlled by pos. Another problem is the a
  ambiguity between \code{rep(x,3)} and \code{rep(x,c(3))} as a special
  case of \code{rep(x,c(3,2))}. If the second is wanted it can be forced by
  \code{rep(x,c(3),NA)} through setting the name argument to NA. 
}
\seealso{
  \code{\link[base]{rep}}
}
\examples{
A <- to.tensor(1:4,c(A=2,B=2))
rep(A,3)
rep(A,3,3,"u")
rep(A,c(2,3))
A <- to.tensor(1:4,c(A=1,B=4))
rep(A,5,pos="A",name=NA)
}
\author{K. Gerald van den Boogaart}
\keyword{arith}
