\name{doLegend}
\alias{doLegend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add a legend to the heatmap}
\description{
Add a simple legend in form of a color bar to a plot.
}
\usage{
doLegend(r, col, side)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{r}{a vector of length two, giving the minimum and maximum value for the legend (think range).
}
  \item{col}{a vector of colors to be shown.
}
  \item{side}{integer between 1 and 4, indicating on which side of the main plot the legend is supposed to be drawn. Standard interpretation: 1 = below, continuing clock-wise.
}
}
\details{
    This is an extremely simple way of giving a visual impression of what numerical values correspond to a given color. The range of values given by \code{r} is split into as many classes of equal length as there are colors in \code{col}. The actual plot is done via a call to \code{\link{image}} and \code{\link{axis}}.
}
\value{The locations of the ticks returned by the call to \code{\link{axis}}
}
\author{Alexander Ploner \email{Alexander.Ploner@ki.se}
}
\seealso{
    \code{\link{plot.annHeatmap}}
}
\examples{
    ## Set up data
    x = matrix(rnorm(100), nrow=10)
    r = range(x)
    ## Set up breaks and colors
    brks = seq(r[1], r[2], 11)
    cols = heat.colors(10)
    
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

