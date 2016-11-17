\name{relist2}
\alias{relist2}
\title{reconstruct the structure of a list}
\description{
  reshapes a vector according to a list template
}
\usage{
relist2(v, l)
}
\arguments{
  \item{v}{vector, probably numeric, of values to reshape}
  \item{l}{template list giving structure}
}
\details{
  attempts to coerce \code{v} into a list with the same
  structure and names as \code{l}
}
\value{
  a list with values corresponding to v and structure corresponding to l
}
\author{Ben Bolker}
\examples{
  l = list(b=1,c=2:5,d=matrix(1:4,nrow=2))
  relist2(1:9,l)
}
\keyword{misc}
