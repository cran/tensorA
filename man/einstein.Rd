\name{einstein.tensor}
\alias{einstein.tensor}
\alias{\%e\%}
\alias{\%e\%.tensor}

\title{Tensor multiplication with Einstein's convention, by summing over
all equally named indices.}
\description{
  Multiplies tensors by multiplying over all duplicate names according
  to Einsteins summing convention by doing an implicit inner product
  over all dimensions with the same name.
}
\usage{
einstein.tensor(...,only=NULL,by=NULL)
## Methods for class tensor
# x \%e\% y
## Default method
# x \%e\% y
}
\arguments{
  \item{\dots}{some tensors, or a renaming code}
  \item{only}{optional list, if given only names in this list are
    automatically processed}
  \item{x}{a tensor}
  \item{y}{a tensor}
  \item{by}{the parallel dimensions}
}
\value{
  the tensor product of all the tensors along all duplicate dimensions.
}
\details{
  see \code{\link{mul.tensor}} on details on tensor
  multiplication. In \code{einstein.tensor} complex operations can be
  performed by command and renaming code: The arguments are processed
  from left to right and multiplied. Unnamed attributes are regarded as
  tensors or scalars and
  multiplied with the current result by the Einstein summing convention,
  which means an inner product over all dimensions with the same
  name. Named attributes can either have the name diag, which performs a
  diagmul according to the same-name convention or be of the form
  \code{A="B"} or \code{"A"="B"}, for which we have two cases. If both
  names are
  present in the current result, an inner multiplication (trace) of on
  these two dimensions is
  performed. If only the first is a name up to this point, the specific
  dimension is renamed to the second name. This renaming might be
  visible in the result or inducing a multiplication according to the
  Einstein convention later. 
}
\references{
}
\seealso{
  \code{\link{mul.tensor}}, \code{\link{to.tensor}}, \code{\link{riemann.tensor}}
}
\examples{
A <- to.tensor(1:20,c(U=2,V=2,W=5))
B <- to.tensor(1:30,list(U=c("a","b","c"),V=c("B1","B2"),W=1:5))
einstein.tensor(A,U="U'",B)
einstein.tensor(A,U="U'",mark(B,"k"))
einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk")
einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",1/10)
einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",diag=to.tensor(c(1,1/10,1/100),c(Uk=3)))

ftable(einstein.tensor(A,U="U'",B))
ftable(einstein.tensor(A,U="U'",mark(B,"k")))
ftable(einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk"))
ftable(einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",1/10))
ftable(einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",diag=to.tensor(c(1,1/10,1/100),c(Uk=3))))

dim(A[[U=~M]])
A[[U=~M]] %e% B
A[[U=~M,V=~"L"]] %e% B

}
\author{K. Gerald van den Boogaart}
\keyword{arith}