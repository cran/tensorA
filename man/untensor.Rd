\name{untensor}
\alias{untensor}
\title{Removes indices/dimensions from a tensor}
\description{
  untensor is more or less the inverse of to.tensor. It flattens
  tensorial dimensions. However the result is still a tensor.
}
\usage{
untensor(X,i=NULL,name=NULL,pos=1,by=NULL)
}
\arguments{
  \item{X}{the tensor}
  \item{i}{the names of the dimensions to be removed and combined to a
    single new 
    one as a character vector or a named list of character vectors if
    the remove should be done in multiple chunks.
    \code{pos} and \code{name} is in this case ignored.}
  \item{name}{the name of the new dimension to replace the others}
  \item{pos}{where to insert the the new dimension}
  \item{by}{if i not given the dimensions to be kept}
}
\value{
  a tensor with the dimensions i removed. 
}
\details{
  The dimensions to be removed are gathered and 
}
\seealso{
  \code{\link{to.tensor}}
}
\examples{
A <- to.tensor(1:64,c(a=2,b=2,c=2,d=2,e=2,f=2))
untensor(A,list(c(1,5),c(2,4)),name=c("i","j"))
untensor(A,by=c("c","f"))
untensor(A,c("a","d"))
}
\author{K.Gerald van den Boogaart}
\keyword{arith}