\name{slice.tensor}
\alias{slice.tensor}
\alias{[[.tensor}
\alias{[.tensor}
\alias{[[<-.tensor}
\alias{[<-.tensor}
\title{Working with the indices of a tensor (accessing, slicing, renaming, ...)}
\description{
  Indexing of tensors allows beside the ordinary selection of ranges
  of indices the renaming of indices. The functions are mainly here to
  keep the the tensor property of the results.
}
\usage{
slice.tensor(X,i,what,drop=FALSE)
## Methods for class tensor
# X[...,drop=TRUE]
# X[...,drop=TRUE] <- value
# X[[...,drop=TRUE]]
# X[[...,drop=TRUE]] <- value
}
\arguments{
  \item{X}{A tensor}
  \item{i}{an index given as number or character}
  \item{what}{levels of the index, a number or a character from
    dimnames}
  \item{drop}{a boolean, if true, indices with only a single level are
    removed}
  \item{\dots}{arguments of the form \code{name=}indices, and for the
    \code{[[ ]]} functions it also allowed to give names from the
    corresponding dimnames \code{name=c("a","b")} to select indices by
    names or \code{name=~newname} to rename dimensions, the first use
    makes a usual array access in the given dimension, where \code{[[ ]]}
    only supports a single index, while \code{[ ]} allows vectors. The
    other type changes the names.} 
}
\value{
  a new tensor with dimensions renamed or individual levels selected
}
\details{
  The functions allow to rename dimensions and to take select a part of
  the tensor.
}
\seealso{
  \code{\link{einstein.tensor}}
}
\examples{
A <- to.tensor(1:20,c(A=2,B=2,C=5))
A[C=1]
A[C=1:3]
A[[B=~b]]                  # renaming dimensions
A[[B=~b,A=~aaa]]      
A[[B=~b,A=~aaa,aaa=1]]      
A[[A=1,B=~gamma]][C=1:2]
A %e% A[[C="D"]]
}
\author{K. Gerald van den Boogaart}
\keyword{arith}