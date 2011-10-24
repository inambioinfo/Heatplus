\name{breakColors}
\alias{breakColors}
\title{Color palette for (symmetric) breaks}
\description{
Given a vector of breaks specifying a set of intervals, this function provides a vector of colors for the indicating the intervals graphically. If the intervals are arranged symmetrically around a specified value, the colors try to reflect this.
}
\usage{
breakColors(breaks, colors, center = 0, tol = 0.001)
}
\arguments{
  \item{breaks}{a vector of breaks}
  \item{colors}{either an explicit vector of colors, or a palette function that takes a number and returns a vector of colors}
  \item{center}{optional center around which to check for symmetry}
  \item{tol}{tolerance (as relative error) for deviation from mathematically exact symmetry}
}
\details{
    The meaning of symmetrical is rather generous here: it is enough that the intervals specified by \code{breaks} are of equal length and that \code{center} is one of the breaks. This means we allow for more or less intervals on one side of \code{center}.
    
    This really only works well if \code{colors} is specified as \code{g2r.colors}, which returns a symmetrical color vector (from green to red) if an even number of colors is requested. The whole point is then that if there are more classes to one side of \code{center} than to the other, this will be reflected by deeper shades of red or green on the appropriate side.
}
\value{
    A vector of colors, of length one less than the number of breaks.
}
\seealso{
    \code{\link{g2r.colors}}
}
\examples{
    ## Fully symmetrical breaks
    br1 = (-3) : 3
    co1 = breakColors(br1, g2r.colors)
    co1
    doLegend(br1, co1, 1)
    
    ## Truncated on one side
    br2 = (-2) : 4
    co2 = breakColors(br2, g2r.colors)
    co2
    doLegend(br2, co2, 1)
    
    ## Does not work with other color schemes
    co3 = breakColors(br2, heat.colors)
    co3
    doLegend(br2, co3, 1)    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

