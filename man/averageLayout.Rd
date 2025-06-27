\name{averageLayout}
\alias{averageLayout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes an average layout over several graphs
}
\description{
This function can be used to compute a joint layout over several graphs.
}
\usage{
averageLayout(..., layout = "spring", repulsion = 1, layout.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Multiple graph objects such as \code{cigraph} objects or weights matrices.
}
  \item{layout}{
Same as in \code{\link{cigraph}}
}
  \item{repulsion}{
    The repulsion parameter as used in  \code{\link{cigraph}}.
  }
  \item{layout.par}{
Same as in \code{\link{cigraph}}
}
}

\value{
A layout matrix
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

