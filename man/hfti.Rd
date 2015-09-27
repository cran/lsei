\name{hfti}
\alias{hfti}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Householder Least Squares Solution under Nonnegativity Constraints}
\description{
  Solves the least squares problem using Householder forward
  triangulation with column interchanges. It is an R interface to the
  HFTI function that is described in Lawson and Hanson
  (1974, 1995). Its Fortran implementation is public domain and is
  available from \url{http://www.netlib.org/lawson-hanson}.
}
\usage{
hfti(a, b, tol=1e-7) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Design matrix.}
  \item{b}{Response vector.} 
  \item{tol}{Tolerance for determining the pseudorank.} 
}
\details{
  Given matrix \code{a} and vector \code{b}, \code{hfti} solves the
  least squares problem:
  
  \deqn{\mathrm{minimize\ \ } || a x - b ||^2.}{minimize  || a x - b ||^2.}
}
\value{
  \item{b}{First \code{krank} elements contains the solution}
  \item{krank}{Psuedo-rank}
  \item{rnorm}{Euclidean norm of the residual vector.}
}
\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{Lawson and Hanson (1974, 1995). Solving least squares
  problems. Englewood Cliffs, N.J., Prentice-Hall.
}

\seealso{ \code{\link{lsei}}, \code{\link{nnls}}. }
\examples{
a = matrix(rnorm(10*4), nrow=10)
b = a \%*\% c(0,1,-1,1) + rnorm(10)
hfti(a, b)
}
\keyword{ array }   % at least one, from doc/KEYWORDS
\keyword{ algebra }   % at least one, from doc/KEYWORDS
