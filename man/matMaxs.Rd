\name{matMaxs}
\alias{matMaxs}

\title{Row or Column Maximum Values of a Matrix}

\description{

  Finds either row or column maximum values of a matrix.

}

\usage{
matMaxs(x, dim=1)
}

\arguments{
  
  \item{x}{numeric matrix.}
  
  \item{dim}{\code{=1}, for row maximum values; \code{=2}, for column
    maximum values. }
 
  
}

\details{
  Matrix \code{x} may contain \code{Inf} or \code{-Inf}, but not
  \code{NA} or \code{NaN}. 
}

\value{
  Returns a numeric vector with row or column maximum values.

  The function is very much the same as using \code{apply(x, 1, max)} or
  \code{apply(x, 2, max)}, but faster.

}

\author{ Yong Wang <yongwang@auckland.ac.nz>}

\examples{
x = cbind(c(1:4,Inf), 5:1)
matMaxs(x)
matMaxs(x, 2)
}
\keyword{ array }   % at least one, from doc/KEYWORDS
\keyword{ algebra }   % at least one, from doc/KEYWORDS
