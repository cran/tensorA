\name{mean.tensor}
\alias{mean.tensor}
\alias{var.tensor}
\title{Mean and variance of tensors}
\description{
  Mean and variance of tensors again tensors. 
}
\usage{
 \method{mean}{tensor}(x,along,...,na.rm=FALSE)
 \method{var}{tensor}(x,y=NULL,...,along,by=NULL,na.rm=FALSE,mark="'")
}
\arguments{
  \item{x}{(set of) dataset(s) of tensors represented by a tensor}
  \item{y}{a second dataset of connected tensors represented by a tensor}
  \item{along}{the indices indexing the datasets}
  \item{\dots}{here for generic compatibility with the compositions package}
  \item{by}{the indices indexing the set of datasets}
  \item{na.rm}{a boolean, if FALSE and missings are in the dataset a
    error is given. If TRUE pairwise exclusion is used.}
  \item{mark}{the to mark the second instance of indices in var(x,...)}
}
\value{
  \item{mean}{gives a tensor like x without the along dimensions
    representing the a mean over all tensors in the dataset. It is
    not necessary to have a by dimension since everything not in along
    is automatically treated parallel}
  \item{var(x,...) }{Gives the covariate tensor representing the covariance
    of x and y. The data tensor indices of x any y should be different,
    since otherwise duplicated names exist in the result.}
  \item{var(x,...) }{Gives the covariate representation of the variance
    of x. All data indices (i.e. all indices neither in by nor in along
  are duplicated. One with and one without the given mark.} 
}
\details{
  Let denote \eqn{a} the along dimension,
  \eqn{i_1,\ldots,i_k}{i_1,...,i_k} and \eqn{j_1,\ldots,j_l}{j_1,...,j_l}
  the data dimension, and
  b the by dimension, then the mean is given by:
  \deqn{M^x_{bi_1,\ldots,i_k}=\frac1{n}\sum_{a}x_{abi_1,\ldots,i_k}}{
    M^x_{bi_1,...,i_k}=\frac1{n}\sum_{a}x_{abi_1,...,i_k}}
  the covariance by
  \deqn{C_{ab i_1,\ldots,i_kj_1,\ldots,j_l}=\frac1{n-1}\sum_{a}
    (x_{abi_1,\ldots,i_k}-M^x_{bi_1,\ldots,i_k})(y_{abj_1,\ldots,j_l}-M^y_{bj_1,\ldots,j_l})}{
C_{ab i_1,...,i_kj_1,...,j_l}=\frac1{n-1}\sum_{a}
    (x_{abi_1,...,i_k}-x_{bi_1,...,i_k})(y_{abj_1,...,j_l}-M^y_{bj_1,...,j_l})}
  and the variance by
  \deqn{V_{ab i_1,\ldots,i_ki'_1,\ldots,i'_l}=\frac1{n-1}\sum_{a}
    (x_{abi_1,\ldots,i_k}-M^x_{bi_1,\ldots,i_k})(x_{abi'_1,\ldots,i'_k}-M^x_{bi'_1,\ldots,i'_l})}{
V_{ab i_1,...,i_kj_1,...,j_l}=\frac1{n-1}\sum_{a}
    (x_{abi_1,...,i_k}-M^x_{bi_1,...,i_k})(x_{abi'_1,...,i'_k}-M^x_{bi'_1,...,i'_l})}
  

}
\seealso{
  \link{tensorA}
}
\examples{
 d1 <- c(a=2,b=2)
 d2 <- c("a'"=2,"b'"=2)
 # a mean tensor:
 m <- to.tensor(1:4,d1)             
 # a positive definite variance tensor:
 V <- delta.tensor(d1)+one.tensor(c(d1,d2))
 V
 # Simulate Normally distributed tensors with these moments:
 X <- (power.tensor(V,c("a","b"),c("a'","b'"),p=1/2)  \%e\%
      to.tensor(rnorm(1000*2*2),c(i=1000,d2))) + m
 # The mean
 mean.tensor(X,along="i")
 # Full tensorial covariance:
 var.tensor(X,along="i")
 # Variance of the slices  X[[b=1]] and X[[b=2]] :
 var.tensor(X,along="i",by="b")
 # Covariance of the slices X[[b=1]] and X[[b=2]] :
 var.tensor(X[[b=1]],X[[a=~"a'",b=2]],along="i")

}
\author{K.Gerald van den Boogaart}
\keyword{multivariate}