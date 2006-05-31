\name{tripledelta.tensor}
\alias{tripledelta.tensor}
\title{A tensor with entry 1 if and only if three indices are equal}
\description{
  The tensor mapping a tensor of dimension d to its corresponding
  diagonal tensor of dimension c(d',d*)
}
\usage{
tripledelta.tensor(d,mark1="'",mark2="*",dn=NULL)
}
\arguments{
  \item{d}{the first of three dimension vectors}
  \item{mark1}{the mark for the second dimension vectors}
  \item{mark2}{the mark for the third dimension vectors}
  \item{dn}{list of character vectors, optional dimnames}
}
\value{
  The tensor given by:
  \deqn{E_{i_1\ldots i_n j_1\ldots
      j_n k_1\ldots k_n}=\delta_{i_1j_1}\delta_{i_1k_1}
    \ldots \delta_{i_nj_n}\delta_{i_nk_1}}{
    E_{i_1...i_n j_1...
      j_n k_1...k_n}=\delta_{i_1j_1}\delta_{i_1k_1}
    ...\delta_{i_nj_n}\delta_{i_nk_1}
  }
}
\details{
  The tripledelta is the tensor mapping a tensor to a corresponding
  diagonal tensor. 
}
\references{
  
}
\seealso{
  \code{\link{delta.tensor}}, \code{\link{diag.tensor}} 
}
\examples{
tripledelta.tensor(3)
}
\author{K. Gerald van den Boogaart}
\keyword{arith}