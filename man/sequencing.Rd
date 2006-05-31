\name{sequencing}
\alias{$.tensor}
\alias{^.tensor}
\alias{|.tensor}
\alias{renamefirst.tensor}
\title{Working with index sequences}
\description{
  In typical tensor notation the indices are not identified by names but
  by positions. The operators allow to identify names and positions
  transparently during calculation. 
}
\usage{
## Methods for class tensor
# x $ y
# x ^ y
# x | y
renamefirst.tensor(x,y)
}
\arguments{
  \item{x}{A tensor }
  \item{y}{
    Typically a character vector specifying a sequence of names for the
    tensor. The names can be specified in various ways:\cr
    The following specifications are equal and specify a sequence of the
    names i,j and k:\cr
    \code{x$ijk}, \code{x$i.j.k}, \code{i.j.k.}, \code{x"$ijk"},
    \code{x^"i.j.k"}, \code{x^c("i","j","k")},\code{x^c("i.j","k")},
    \code{x^c("$i.j","k")},\code{x^c("$ij","k")},
    \code{x^c("$","ijk")}\cr
    In general names are separated by dots.
    All notations with \code{\$} either as operator or as the first
    character of the first string allow to omit the dots assuming that
    all names are single character. If any dot is present all dots must
    be given. The difference of \code{\$} and \code{\^} is that the
    first accepts a name and the second an character valued expression.
    \cr
    Multi letter indices like "alpha","beta","gamma"
    can only be given in the dot-free version of the notation
    making the following
    specifications equal:
    \code{x$alpha.beta.gamma}, \code{alpha.beta.gamma.},
    \code{x^"$alpha.beta.gamma"},
    \code{x^"alpha.beta.gamma"}, \code{x^c("alpha","beta","gamma")},
    \code{x^c("alpha.beta","gamma")},
    \code{x^c("$alpha.beta","k")},
    \code{x^c("$","alpha.beta.gammak")}\cr
    The specification for \code{|} is equal to that for \code{^}.
  }
}
\value{
  A tensor of the same shape as x but with reordered dimensions (for
  \code{|}) or renamed dimensions (for the others) 
}
\details{
  These functions are used to mimic the mathematical notation in tensor
  analysis.
  Formulae of the form (with Einstein convention):
  \deqn{E_{ijk}= A_{ihl}C_{hj}C_{lk}}{E_ijk= A_ihl C_hj C_lk }
  with defined tensors \eqn{A_{ijk}}{A_ijk} and \eqn{C_{ij}}{C_ij} can
be given the
  simple
  form \cr
  \code{ E <- A$ihl \%e\% C$hj \%e\% C$lk |"$ijk"}\cr
  or alternatively for multi letter names:\cr
  \code{ E <- A$i.h.l \%e\% C$h.j \%e\% C$l.k |"i.j.k"}\cr
  or more flexible in computation with arguments I,J,K:\cr
  \code{ E <- A^c(I,"h.l") \%e\% C^c("h",J) \%e\% C^c("l",K) | c(I,J,K)}\cr
  The \code{$} or \code{^} binds to the tensors with high precedence
  and renames the first elements. The \code{|} binds with very low
  precedence and reorders the tensor according to the
  assumed index sequence of the result afterwards.

}
\references{
}
\seealso{
  \code{\link{reorder.tensor}}, \code{\link{names<-.tensor}}, \code{\link{[[.tensor}}
}
\examples{
A <- to.tensor(1:20,c(i=5,j=2,k=2))
C <- to.tensor(1:4,c(i=2,j=2))
E <- A$ihl \%e\% C$hj \%e\% C$lk |"$ijk"
E
# Same as:
E2 <- reorder.tensor(A[[j=~h,k=~l]] \%e\% C[[i=~h]] \%e\% C[[i=~l,j=~k]],c("i","j","k"))
E-E2
E <- A$i.h.l \%e\% C$h.j \%e\% C$l.k |"i.j.k"
E
E-E2
E <- A^"i.h.l" \%e\% C^"h.j" \%e\% C^"l.k" |"i.j.k"
E
E-E2
}
\author{K. Gerald van den Boogaart}
\keyword{arith}