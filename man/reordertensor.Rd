\name{reorder.tensor}
\alias{reorder.tensor}
%\alias{reorder.tidx}
\title{Permutation of indices and storage sequence of a tensor}
\description{
  This permutes tensor dimensions like aperm. However the interface is
  more flexible since not all dimensions have to given and names can be
  used instead of numbers.
}
\usage{
\method{reorder}{tensor}(x,i=NULL,...,by=NULL)
%reorder.tidx(x,i,...)
}
\arguments{
  \item{x}{the tensor}
  \item{i}{numeric or character giving dimensions intended to come
    first}
  \item{\dots}{further arguments to other instances of the generic function}
  \item{by}{the complement of i, if i is not given}
}
\value{
  \code{reorder.tensor} returns a tensor equal to x but stored with a different
  sequence of dimensions.
  %\code{reorder.tidx} returns the the
  %permutation of the memory sequence to achieve that as a numeric
  %vector. 
}
\details{
  the remaining dimensions keep their relative sequence and follow at
  the end of the dimension attribute.
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- to.tensor(1:20,c(A=2,B=2,C=5))
A
reorder(A,"C")
reorder(A,"B")
}
\author{K.Gerald v.d. Boogaart}
\keyword{array}