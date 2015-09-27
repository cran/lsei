\name{lsei}
\alias{lsei}
\alias{lsi}
\alias{ldp}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Least Squares Solution under Equality and Inequality Constraints}
\description{
  The \code{lsei} function solves the least squares problem under
  both equality and inequality constraints and implements the LSEI
  algorithm described in Lawson and Hanson (1974,
  1995).

  The \code{lsi} function solves the least squares problem under
  inequality constraints and implements the LSI algorithm described in
  Lawson and Hanson (1974, 1995).
  
  The \code{ldp} function solves the least distance programming problem
  under inequality constraints and implements the LDP algorithm
  described in Lawson and Hanson (1974, 1995).
  
  The NNLS Fortran implementation used internally is downloaded from
  \url{http://www.netlib.org/lawson-hanson}.
}

\usage{
lsei(a, b, c, d, e, f)
lsi(a, b, e, f)
ldp(e, f, tol=1e-15)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Design matrix.}
  \item{b}{Response vector.} 
  \item{c}{Matrix for equality constraint.}
  \item{d}{Vector for equality constraint.} 
  \item{e}{Matrix for inequality constraint.}
  \item{f}{Vector for inequality constraint.} 
  \item{tol}{Tolerance for checking compatibility of inequalities.} 
}

\details{

  Given matrices \code{a}, \code{c} and \code{e}, and vectors \code{b},
  \code{d} and \code{f}, function \code{lsei} solves the least squares
  problem under both equality and inequality constraints:
  
  \deqn{\mathrm{minimize\ \ } || a x - b ||^2,}{minimize  || a x - b ||^2,} \deqn{\mathrm{subject\ to\ \ } c x = d, e x \ge f.}{subject to  c x = d, e x >= f.}
  
  Function \code{lsi} solves the least squares problem under inequality
  constraints:

  \deqn{\mathrm{minimize\ \ } || a x - b ||^2,}{minimize  || a x - b ||^2,} \deqn{\mathrm{\ \ \ subject\ to\ \ } e x \ge f.}{subject to  e x >= f.} 
  
  Function \code{ldp} solves the least distance programming problem
  under inequality constraints:

  \deqn{\mathrm{minimize\ \ } || x ||^2,}{minimize  || x ||^2,} \deqn{\mathrm{\ \ \ subject\ to\ \ } e x \ge f.}{subject to  e x >= f.}
  
}
\value{
  A vector of the solution values
}
\author{ Yong Wang <yongwang@auckland.ac.nz>}

\references{Lawson and Hanson (1974, 1995). Solving least squares
  problems. Englewood Cliffs, N.J., Prentice-Hall.
}

\seealso{ \code{\link{nnls}},\code{\link{hfti}}. }
\examples{
beta = c(rnorm(2), 1)
beta[beta<0] = 0
beta = beta/sum(beta)
a = matrix(rnorm(18), ncol=3)
b = a \%*\% beta + rnorm(3,sd=.1)
c = matrix(rep(1, 3), nrow=1)
d = 1
e = diag(rep(1,3))
f = rep(0,3)
lsei(a, b, c, d, e, f)

## Example from Lawson and Hanson (1974), p.170
a = cbind(c(.25,.5,.5,.8),rep(1,4))
b = c(.5,.6,.7,1.2)
e = cbind(c(1,0,-1),c(0,1,-1))
f = c(0,0,-1)
lsi(a, b, e, f)      # Solution: 0.6213152 0.3786848

## Example from Lawson and Hanson (1974), p.171:
e = cbind(c(-.207,-.392,.599), c(2.558, -1.351, -1.206))
f = c(-1.3,-.084,.384)
ldp(e, f)            # Solution: 0.1268538 -0.2554018
}
\keyword{ array }   % at least one, from doc/KEYWORDS
\keyword{ algebra }   % at least one, from doc/KEYWORDS
