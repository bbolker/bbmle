\name{namedrop}
\alias{namedrop}
\title{drop unneeded names from list elements}
\description{
  goes through a list (containing a combination
  of single- and multiple-element vectors) and
  removes redundant names that will make trouble
  for mle
}
\usage{
namedrop(x)
}
\arguments{
  \item{x}{a list of named or unnamed, typically numeric, vectors}
}
\details{
  examines each element of \code{x}.  If the element has length
  one and is a named vector, the name is removed; if \code{length(x)}
  is greater than 1, but all the names are the same, the vector
  is renamed 
}
\value{
  the original list, with names removed/added
}
\author{Ben Bolker}
\examples{
x = list(a=c(a=1),b=c(d=1,d=2),c=c(a=1,b=2,c=3))
names(unlist(namedrop(x)))
names(unlist(namedrop(x)))
}
\keyword{misc}
