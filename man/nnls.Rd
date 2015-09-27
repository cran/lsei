\name{nnls}
\alias{nnls}
\alias{pnnls}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Least Squares under Nonnegativity Constraints}
\description{
  
  The \code{nnls} function solves the least squares problem under
  nonnegativity (NN) constraints. It is an R interface to the NNLS
  function that is described in Lawson and Hanson (1974, 1995). Its
  Fortran implementation is public domain and available at
  \url{http://www.netlib.org/lawson-hanson}.

  The \code{pnnls} function solves the least squares problem when
  possibly only some coefficients are subject to nonnegativity
  constraints. It also works when the NN-restricted coefficients are
  further restricted to have a fixed positive sum.
}

\usage{
nnls(a, b) 
pnnls(a, b, k=0, sum=NULL) 
}

\arguments{
  
  \item{a}{Design matrix.}
  
  \item{b}{Response vector.} 
  
  \item{k}{Integer, meaning that the first \code{k} coefficients are not
  NN-restricted.} 
  
  \item{sum}{= NULL, if NN-restricted coefficients are not further
    restricted to have a fixed sum; 

    = a positive value, if NN-restricted coefficients are further
    restricted to have a fixed positive sum.}
  
}

\details{
  Given matrix \code{a} and vector \code{b}, \code{nnls} solves the
  nonnegativity least squares problem:
  
  \deqn{\mathrm{minimize \ \ } || a x - b ||^2,}{minimize  || a x - b ||^2,} \deqn{\mathrm{\ \ \ subject\ to\ \ } x \ge 0.}{ subject to  x >= 0.}
  
  Function \code{pnnls} can also solve the above nonnegativity least squares
  problem when \code{k=0}, but it may leave the first \code{k} values
  unrestricted. The returned value of \code{k} can be different from the
  input, if \code{a} has linearly dependent columns. If \code{sum} is a
  positive value, \code{pnnls} solves the problem by further restricting that
  the NN-restricted coefficients must sum to the given value.
  
}

\value{
  
  \item{x}{Solution}
  
  \item{rnorm}{Euclidean norm of the residual vector.}
  
  \item{index}{Indices of the columns; those in the positive set are
    first given, and then those in the zero set.}
  
  \item{mode}{= 1, successful computation;

    = 2, bad dimensions of the problem;

    = 3, iteration count exceeded (more than 3 times the number
    of variables iterations).}
  
  \item{r}{The upper-triangular matrix \code{Q*a}, pivoted by variables
    in the order of \code{index}.}
  
  \item{b}{The vector \code{Q*b}, pivoted by variables in the order of
    \code{index}.}
  
}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{

  Lawson and Hanson (1974, 1995). Solving least squares
  problems. Englewood Cliffs, N.J., Prentice-Hall.

  Dax (1990). The smallest point of a polytope. Journal of Optimization
  Theory and Applications, 64, pp. 429-432.
  
  Wang (2010). Fisher scoring: an interpolation family and its Monte
  Carlo implementations. Computational Statistics and Data Analysis, 54,
  pp. 1744-1755.
  
}

\seealso{ \code{\link{lsei}}, \code{\link{hfti}}. }

\examples{
a = matrix(rnorm(40), nrow=10)
b = a \%*\% c(0,1,-1,1) + rnorm(10)
nnls(a, b)$x                   # constraint x >= 0,
pnnls(a, b, k=0)$x             # same as nnls(a, b)
pnnls(a, b, k=2)$x             # first two coeffs are not NN-constrained
pnnls(a, b, k=2, sum=1)$x      # NN-constrained coeffs must sum to 1
pnnls(a, b, k=2, sum=2)$x      # NN-constrained coeffs must sum to 2
}
\keyword{ array }   % at least one, from doc/KEYWORDS
\keyword{ algebra }   % at least one, from doc/KEYWORDS
