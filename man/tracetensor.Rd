\name{trace.tensor}
\alias{trace.tensor}
\title{Collapse a tensor}
\description{
  Collapses the tensor over dimensions i and j. This is like a trace for
  matrices or like an inner product of the dimensions i and j.
}
\usage{
trace.tensor(X,i,j)
}
\arguments{
  \item{X}{the tensor}
  \item{i}{a numeric or character vector of dimensions of \code{X},
    used for the inner product.
  }
  \item{j}{a numeric or character vector of dimensions of \code{X} with
    the same length but other elements than i.  }
}
\value{
  A tensor like X with the i and j dimensions removed.
}
\details{
  Let be
  \deqn{X_{i_1\ldots i_n j_1\ldots j_n k_1 \ldots k_d}}{
    X_{i_1...i_n j_1...j_n k_1...k_d}
  } the tensor. Then the result is given by
  \deqn{E_{k_1 \ldots k_d}=sum_{i_1\ldots i_n} X_{i_1\ldots i_n
      i_1\ldots i_n k_1 \ldots k_d} }{
    E_{k_1...k_d}=sum_{i_1...i_n} X_{i_1...i_n
      i_1...i_n k_1...k_d} 
  }
  With the Einstein summing convention we would write:
  \deqn{
    E_{k_1 \ldots k_d}=X_{i_1\ldots i_n j_1\ldots j_n k_1 \ldots
      k_d}\delta_{i_1j_1}\ldots \delta_{i_nj_n}{
      E_{k_1...k_d}=X_{i_1...i_n j_1...j_n k_1 ...
      k_d}\delta_{i_1j_1} ... \delta_{i_nj_n
      }
  }}
}
\seealso{
  \code{\link{mul.tensor}}, \code{\link{to.tensor}}
}
\examples{
A <- to.tensor(1:20,c(i=2,j=2,k=5))
A
trace.tensor(A,"i","j")
}
\author{K. Gerald van den Boogaart}
\keyword{arith}