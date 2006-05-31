\name{tensorA-package}
\alias{tensorA-package}
\alias{tensor}
\alias{Tensor}
\alias{tensorA}
\docType{package}
\title{The tensorA package for tensor arithmetic}
\description{
  tensorA stands for "tensor arithmetic". A tensor is a mathematical
  generalization of vector and matrix with
  many applications in physics, geometry and in the statistics of
  vectors valued data. However the package is also useful in any case,
  where computations on sequences of matrices, vectors or even tensors
  is involved. 
}
\details{
  \tabular{ll}{
    Package: \tab tensorA\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2006-06-08\cr
    License: \tab GPL Version 2 or newer\cr
  }

  The tensorA package is made to allow programming for tensors in R on
  the same level of abstraction as we know from matrices. It
  provides many of the
  mathematical operations common in tensor arithmetics including the
  whole tensor calculus of covariate and contravariate indices, naming
  of indices, sequence of indices, decompositions of tensors, Einstein
  and Riemann summing conventions and vectorized computations on datasets
  of tensors just like the well vectorization of numbers in R. It
  provides tools to write tensor formulae very close to there paper form
  and to handle tensors of arbitrary level with simple programs.
  \cr
  The whole documentation of the package is best read in
  pdf or dvi format since it contains complicated mathematical formulae
  with multi-indices.
  \cr
  
  Simply speaking a tensor (see \code{\link{to.tensor}}) is just a
  multidimensional array
  \code{A[,,]}. The number
  of indices (i.e. \code{length(dim(A))} is called the level of the
  tensor (see \code{\link{level.tensor}}). A tensor is mathematically it
  is denoted
  by a core symbol (e.g. A) with multiple indices:e.g.
  \deqn{A_{ijk}}
  The indices \eqn{i,j,k} can be seen as names for the dimensions
  and as integer numbers giving the respective index into the array.
  However the tensor is an algebraical object with many algebraical
  operations defined on it, which are also of relevancy for programming,
  e.g. in the parallel treatment of multiple linear equation systems.

  To understand the package we need to understand tensors including
  their mathematical origin, the corresponding calculus, notation and
  basic
  operations.
  \cr
  One mathematical interpretation of a tensor, which is most relevant
  for physics, that of a multi-linear
  form of \eqn{level(A)} vectors, i.e. a function of \eqn{level(A)} many
  vectors
  to the real or complex numbers, which is linear with respect to each
  of its arguments. E.g. the two vectors "plane face direction" and "force
  direction" are mapped to the actual force by the stress tensor.
  \cr
  Row vectors are a special case of that and likewise column vectors as
  linear forms for row vectors. Matrices are bilinear forms of a row
  vector and a column vector. Thus Vectors and Matrices are examples of
  tensors of level 1 and 2. 
  
  Another interpretation of a tensor is the that of a linear mapping,
  quite like a matrix, but from a tensor space (e.g. the space of
  matrices or vectors seen as tensor) to another tensor space
  (e.g. again a space of matrices). An
  example for that is the Hook elasticity tensor mapping the strain
  tensor (i.e. a
  matrix describing the local deformation) to the stress tensor (i.e. a
  matrix describing the local forces). The Hook tensor is a tensor of
  level 4. Statistically relevant tensors of level 4 are
  e.g. covariances of matrices mapping two linear forms (i.e. 2 level 2
  tensors) on observed
  matrices to there covariance.  The mapping is performed with the
  tensor product, which is not unlike a matrix product, however more
  general. Let denote \eqn{A} a matrix and \eqn{v} a vector, we
  would write \eqn{r=Ab} for the matrix product and \code{r <- A\%*\%b} in R,
  which is defined as:
  \deqn{r_i = \sum_{j=1}^{j_{\max}} A_{ij}b_j }{r_i = \sum_{j=1}^{j_{\max}} A_{ij}b_j }
  We know that we have to use the \(j\)-dimension in the summing, since
  the matrix multiplication rule says "row times column". 
  Since a tensor can have more than two indices there is no row or
  column specified and we need to specify the inner product
  differently. To do this in
  the Einstein-Notation writing the tensor always with indices
  \eqn{r_i=A_{ij}b_j} and according to the Einstein summing rule the
  entries of \(r_i\) are given by an implicit sum over all indices which
  show
  up twice in this notation:
  \deqn{ r_i=\sum_{j=1}^{j_{\max}} A_{ij}b_j }{ r_i=\sum_{j=1}^{j_{\max}}
    A_{ij}b_j }
  
  This notation allows for a multitude of other products:
  \eqn{ A_{ij}b_i=t(A)b }, \eqn{ A_{ij}b_k=outer(A,b)
  }{ A_{ij}b_k=outer(A,b) } ,
  \eqn{ A_{ii}b_j=trace(A)b }{ A_{ii}b_j=trace(A)b } with equal
  simplicity and without any
  additional functions. 
  More complicated products involving more than tensors of level two
  can not even be formulated in pure matrix algebra without
  re-dimensioning of arrays e.g. \eqn{b_ib_jb_k}{b_ib_jb_k},
\eqn{A_{ijk}b_j}{A_{ijk}b_j}. The
  Einstein summing rule is implemented in
  \code{\link{einstein.tensor}} and supported by the index sequencing
  functions \code{\link{$.tensor}} and \code{\link{|.tensor}}. A general
  multiplication allowing to
  identify and sum over any two indices is implemented in
  \code{\link{trace.tensor}}, when the indices are in the same tensor
  and in \code{\link{mul.tensor}}, when the indices to sum over are in
  different tensors. 
  \cr
  Tensors with the same level and dimensions (identified by name
  and dimension) can also be added like matrices according to the rule
  that the values with the same combination of index values are added
  (see \code{\link{add.tensor}}). The implementation takes care of the
  sequence of the indices and rearranges them accordingly to match
  dimensions with the same name. E.g. the tensor addition
  \deqn{E_ijk=A_{ijk}+B_{kji}}{E_ijk=A_ijk+B_kji} 
  has the effect, which is expressed by the same formula read in
  entries, which is also true for the more surprising
  \deqn{E_ijk=A_{ij}+B_{kj}}{E_ij=A_ijk+B_kj} 
  \cr
  Like a matrix a tensor can also be seen as a mapping from one tensor
  space to another:
  \deqn{A_{i_1\ldots i_d j_1 \ldots j_e}x_{j_1 \ldots j_e}=b_{i_1\ldots
      i_d}}{A_{i_1...i_d j_1...j_e}x_{j_1...j_e}=b_{i_1...i_d}
  }
  In this reading all the standard matrix computations and
  decompositions get a tensorial interpretation and generalization. The
  package provides some of these (see \code{\link{svd.tensor}}). 
  \cr
  Another interpretation of tensors is as a sequence of tensors of lower
  level. E.g. a data matrix is seen as a sequence of vectors in
  multivariate dataset. The tensorA library provides means to do
  computation on these in parallel on these sequences of tensors like we
  can do parallel computation on sequences of numbers. This is typically
  done by the \code{by=} argument present in most functions and giving
  the index enumerating the elements of the sequence.\cr
  E.g. If we have
  sequence \eqn{V_{ijd}} of variance matrices \eqn{V_{ij}} of some sequence
  \eqn{v_{id}} of vectors and we would like to transform the vectors
  with some Matrix \eqn{M_{i'i}} we would get the sequence of
  transformed variances
  by \eqn{V_{ijd} M_{i'i}M_{j'j}}.
  However if the \eqn{M_{ki}} are different for each of the elements
  in sequence we would have stored them in a tensor \eqn{M_{kid}} and
  would have to replace  \eqn{M_{kid}} with \eqn{M_{kidd'}=M_{kid}} if
  \eqn{d=d'} and zero otherwise. We can than get our result by
  \deqn{V_{ijd}M_{i'id'd}M_{j'jd'd''}} and we would have a by dimension
  of
  \code{by="d"}. These operations are not strictly
  mathematical tensor operation, but generalizations of the
  vectorization approach of R. This is also closely related to
  \code{\link{diagmul.tensor}} or \code{\link{diag.tensor}}.   
  \cr
  To complicate things the Einstein rule is only valid in case of
  tensors represented with respect to a orthogonal basis. Otherwise
  tensors get lower and upper indices like
  \deqn{A_{i\cdot k}^{\cdot j \cdot}}
  for representation in the covariate and contravariate form of the
  basis. In this case the Riemann summing rule applies which only sums
  over pairs of the same index, where one is in the lower and one is in
  the upper position. The contravariate form is represented with indices
  prefixed by \code{^}. 
  \cr 
  The state of being covariate or contravariate can be changed by the
  dragging rule, which allows to switch between both state through the
  multiplication with the geometry tensors \eqn{g_i^{\;j}}{g_i^j}. This
  can be done through \code{\link{drag.tensor}}.
}
\references{
}
\seealso{
  \code{\link{to.tensor}}, \code{\link{mul.tensor}} ,
  \code{\link{einstein.tensor}}, \code{\link{add.tensor}},
  \code{\link{[[.tensor}}, \code{\link{|.tensor}} 
}
\examples{
A <- to.tensor( 1:20, c(a=2,b=2,c=5) )
A
ftable(A)
B <- to.tensor( c(0,1,1,0) , c(a=2,"a'"=2))
A \%e\% B
drag.tensor( A , B, c("a","b"))
A \%e\% one.tensor(c(c=5))/5     # a mean of matrices
reorder.tensor(A,c("c","b","a"))
A -  reorder.tensor(A,c("c","b","a"))  # =0 since sequence is irrelevant
inv.tensor(A,"a",by="c")  

}
\author{K.Gerald van den Boogaart <boogaart@uni-greifswald.de}
\keyword{package}