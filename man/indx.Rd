\name{indx}
\alias{indx}

\title{Index-finding in a Sorted Vector}

\description{

  For each of given values, \code{indx} finds the index of the
  value in a vector sorted in ascending order that the given value is
  barely greater than or equal to.

}

\usage{
indx(x, v)
}

\arguments{
  
  \item{x}{vector of numeric values, the indices of which are to be found.}
  
  \item{v}{vector of numeric values sorted in ascending order.}
 
  
}

\details{
  For each x[i], the function returns integer j such that 
  \deqn{v_j \le x_i < v_{j+1}}{v[j] <= x[i] < v[j+1]}
  where \eqn{v_0 = - \infty \mathrm{ and } v_{n+1} = \infty}{v[0] = -Inf and v[n+1] = Inf}.

}

\value{
  Returns a vector of integers, that are indices of x-values in vector v.

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\examples{
indx(0:6,c(1:5,5))
indx(sort(rnorm(5)),-2:2)
}
\keyword{ array }   % at least one, from doc/KEYWORDS
\keyword{ algebra }   % at least one, from doc/KEYWORDS
