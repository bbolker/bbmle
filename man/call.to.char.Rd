\name{call.to.char}
\alias{call.to.char}
\title{Convert calls to character}
\description{
  Utility function (hack) to convert calls such
  as y~x to their character equivalent
}
\usage{
call.to.char(x)
}
\arguments{
  \item{x}{a formula (call)}
}
\details{
  It would be nice if \code{as.character(y~x)}
  gave "y~x", but it doesn't, so this hack achieves
  the same goal
}
\value{
  a character vector of length 1
}
\author{Ben Bolker}
\examples{
as.character(y~x)
call.to.char(y~x)
}
\keyword{misc}
