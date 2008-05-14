\name{solve.tensor}
\alias{solve.tensor}
\title{Solving linear equations with tensors}
\description{
	We can formulate linear equation systems with tensors. This functions
	solves these systems or gives a least squares fit of minimal
	norm.
}
\usage{
solve.tensor(a,b,i,j=i,...,allowSingular=FALSE,eps=1E-10,by=NULL)
}
\arguments{
  \item{a}{The a of ax=b}
  \item{b}{The a of ax=b}
  \item{i}{The dimensions of the equation in a}
  \item{j}{The dimensions of the equation in b}
  \item{allowSingular}{A boolean, indicating the that a least squares
    fit should be generated with singular equations systems.}
  \item{\dots}{further arguments for generic use}
  \item{eps}{The limit for the smallest singular value in inversion}
  \item{by}{the operation is done in parallel for these dimensions}
}
\value{
  a tensor such that ax=b as good as possible for each combination of by
  values.
}
\details{
  A tensor can be seen as a linear mapping of a tensor to a tensor. Let
  denote \eqn{R_i} the space of real tensors with dimensions
  \eqn{i_1...i_d}.
  \itemize{
  \item{solve.tensor}{Solves the equation for
  \eqn{a_{i_1...i_dk_1...k_p}}, \eqn{b_{j_1...j_dl_1...l_q}} and
  \eqn{x_{k_1...k_pl_1...l_q}} the equation  
  \deqn{ \sum_{k_1,...,k_p}
    a_{i_1...i_dk_1...k_p}x_{k_1...k_pl_1...l_q}=
    b_{j_1...j_dl_1...l_q}}.
  }
  }
}
\seealso{
  \code{\link{to.tensor}}, \code{\link{svd.tensor}},
  \code{\link{inv.tensor}}, \code{\link{chol.tensor}},
  \code{\link{power.tensor}} 
}
\examples{
R1  <- matrix(rnorm(9),nrow=3)
R1i <- solve(R1)
R2 <- to.tensor(R1,c(a=3,b=3),what=1:2)
R2i <- to.tensor(R1i,c(b=3,a=3),what=1:2)

inv.tensor(R2,"a","b") - R2i
inv.tensor(R2,"a","b",allowSingular=TRUE) - R2i

inv.tensor(rep(R2,4,1,"K"),"a","b",by="K") - rep(R2i,4,1,"K")
inv.tensor(rep(R2,4,1,"K"),"a","b",by="K",allowSingular=TRUE) - rep(R2i,4,3,"K")

R3 <- to.tensor(rnorm(15),c(a=3,z=5))

mul.tensor(R2i,"b",mul.tensor(R2,"a",R3)) # R3

solve.tensor(R2i,R3[[z=1]],"a")
mul.tensor(R2,"a",R3[[z=1]])

solve.tensor(R2i,R3,"a")
mul.tensor(R2,"a",R3)

solve.tensor(R2i,R3[[z=1]],"a",allowSingular=TRUE)
mul.tensor(R2,"a",R3[[z=1]])

solve.tensor(R2i,R3,"a",allowSingular=TRUE)
mul.tensor(R2,"a",R3)

solve.tensor(rep(R2i,4,1,"K"),R3[[z=1]],"a",by="K")
rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K")

solve.tensor(rep(R2i,4,1,"K"),rep(R3[[z=1]],4,1,"K"),"a",by="K")
rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K")

}
\author{K. Gerald van den Boogaart}
\keyword{arith}
