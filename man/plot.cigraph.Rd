\name{plot.cigraph}
\alias{plot.cigraph}

\title{
Plot method for "cigraph"
}
\description{
Plots an object created by \code{\link{cigraph}}.
}
\usage{
\method{plot}{cigraph}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A "cigraph" object
}
  \item{\dots}{
Not used
}
}

\details{
If the result of \code{\link{cigraph}} is stored, such as \code{Graph <- cigraph(...)}, the plot can be recreated in two ways. \code{cigraph(Graph, ...))} reruns \code{cigraph} with the same arguments used in the origina call except those restated in the dots. For example \code{cigraph(Graph, shape = "square")} will recreate the same plot but now use square nodes instead of circular. \code{plot(Graph)} will NOT rerun \code{cigraph} but simply plot the cigraph object. This means that now specific graph attributes can be changed before plotting.

More specific, \code{cigraph(Graph)} will base the new plot on the \code{Arguments} element of the \code{cigraph} object and \code{plot(cigraph)} will base the new plot on the \code{graphAttributes} element of cigraph.

}

\author{
Sacha Epskamp (mail@sachaepskamp.com)
}
