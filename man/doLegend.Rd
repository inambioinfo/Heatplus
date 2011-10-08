\name{doLegend}
\alias{doLegend}
\title{A simple legend}
\description{
Add a simple legend in form of a color bar to a plot.
}
\usage{
doLegend(breaks, col, side)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{breaks}{a vector of breaks defining a set of intervals for the data
}
  \item{col}{a vector of colors corresponding to the intervals.
}
  \item{side}{integer between 1 and 4, indicating on which side of the main plot the legend is supposed to be drawn. Standard interpretation: 1 = below, continuing clock-wise.
}
}
\details{
    This is an extremely simple way of giving a visual impression of what numerical values correspond to a given color. The actual plot is done via a call to \code{\link{image}} and \code{\link{axis}}.
}
\value{The locations of the ticks returned by the call to \code{\link{axis}}
}
\seealso{
    \code{\link{plot.annHeatmap}}, \code{\link{niceBreaks}}, \code{\link{g2r.colors}}
}
\examples{
    ## Set up data
    doLegend(1:9, g2r.colors(8), 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

