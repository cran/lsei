\name{lsei}
\alias{lsei}
\alias{lsi}
\alias{ldp}
\alias{qp}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Least Squares Solution under Equality and Inequality Constraints}
\description{
  The \code{lsei} function solves a least squares problem under
  both equality and inequality constraints and implements the LSEI
  algorithm described in Lawson and Hanson (1974,
  1995).

  The \code{lsi} function solves a least squares problem under
  inequality constraints and implements the LSI algorithm described in
  Lawson and Hanson (1974, 1995).
  
  The \code{ldp} function solves a least distance programming problem
  under inequality constraints and implements the LDP algorithm
  described in Lawson and Hanson (1974, 1995).
  
  The \code{qp} function solves a quadratic programming problem, by
  transforming the problem into a least squares one under the same
  equality and inequality constraints, which is then solved by function
  \code{lsei}.
  
  The NNLS Fortran implementation used internally is downloaded from
  \url{http://www.netlib.org/lawson-hanson}.
}

\usage{
lsei(a, b, c=NULL, d=NULL, e=NULL, f=NULL, lower=-Inf, upper=Inf)
lsi(a, b, e=NULL, f=NULL, lower=-Inf, upper=Inf)
ldp(e, f, tol=1e-15)
qp(q, p, c=NULL, d=NULL, e=NULL, f=NULL, lower=-Inf, upper=Inf, tol=1e-15)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Design matrix.}
  
  \item{b}{Response vector.}
  
  \item{c}{Matrix of numeric coefficients on the left-hand sides of
    equality constraints. If it is NULL, \code{c} and \code{d} are
    ignored.} 
  
  \item{d}{Vector of numeric values on the right-hand sides of equality
    constraints.}
  
  \item{e}{Matrix of numeric coefficients on the left-hand sides of
    inequality constraints. If it is NULL, \code{e} and \code{f} are
    ignored.} 
  
  \item{f}{Vector of numeric values on the right-hand sides of
    inequality constraints.}
  
  \item{q}{Matrix of numeric values for the quadratic term of a
    quadratic programming problem.}
  
  \item{p}{Vector of numeric values for the linear term of a quadratic
    programming problem.}
  
  \item{lower, upper}{Bounds on the solutions, as a way to specify such
    simple inequality constraints.}
  
  \item{tol}{Tolerance, for checking compatibility of inequalities in
    \code{lsi} and for calculating pseudo-rank in \code{qp}.}

}

\details{

  Given matrices \code{a}, \code{c} and \code{e}, and vectors \code{b},
  \code{d} and \code{f}, function \code{lsei} solves the least squares
  problem under both equality and inequality constraints:
  
  \deqn{\mathrm{minimize\ \ } || a x - b ||,}{minimize  || a x - b ||,} \deqn{\mathrm{subject\ to\ \ } c x = d, e x \ge f.}{subject to  c x = d, e x >= f.}
  
  Function \code{lsi} solves the least squares problem under inequality
  constraints:

  \deqn{\mathrm{minimize\ \ } || a x - b ||,}{minimize  || a x - b ||,} \deqn{\mathrm{\ \ \ subject\ to\ \ } e x \ge f.}{subject to  e x >= f.} 
  
  Function \code{ldp} solves the least distance programming problem
  under inequality constraints:

  \deqn{\mathrm{minimize\ \ } || x ||,}{minimize  || x ||,} \deqn{\mathrm{\ \ \ subject\ to\ \ } e x \ge f.}{subject to  e x >= f.}
  
  Function \code{qp} solves the quadratic programming problem:
  
  \deqn{\mathrm{minimize\ \ } \frac12 x^T q x + p^T x,}{minimize  0.5 x^T q x + p^T x,} \deqn{\mathrm{subject\ to\ \ } c x = d, e x \ge f.}{subject to  c x = d, e x >= f.}

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
beta = beta / sum(beta)
a = matrix(rnorm(18), ncol=3)
b = a \%*\% beta + rnorm(3,sd=.1)
c = t(rep(1, 3))
d = 1
e = diag(1,3)
f = rep(0,3)
lsei(a, b)                        # under no constraint
lsei(a, b, c, d)                  # under eq. constraints
lsei(a, b, e=e, f=f)              # under ineq. constraints
lsei(a, b, c, d, e, f)            # under eq. and ineq. constraints
lsei(a, b, rep(1,3), 1, lower=0)  # same solution
q = crossprod(a)
p = -drop(crossprod(b, a))
qp(q, p, rep(1,3), 1, lower=0)    # same solution

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
