\name{g2r.colors}
\alias{g2r.colors}
\title{Palette from green to red via black
}
\description{
Returns a color vector of the requested length, ranging from pure red to pure green via slighlty tinted black.
}
\usage{
g2r.colors(n = 12, min.tinge = 0.33)
}
\arguments{
    \item{n}{the number of requested colors}
    \item{min.tinge}{the proportion of red/green added to black to make it recognizably green or red}
}
\details{
    If \code{n} is even, the colors range from pure green to green-tinted black to red-tinted black to pure red. If \code{n} is odd, the colors range from pure red to pure green, with full black for the median class.
}
\value{
    A vector of (RGB-) colors of the specified length
}
\seealso{
    \code{\link{breakColors}}
}
\examples{
    ## Even number: residual tint shows left/right of center
    co_even = g2r.colors(10)
    co_even
    doLegend(1:11, co_even, 1)
    
    ## Odd number: central class all black
    co_odd = g2r.colors(9)
    co_odd
    doLegend(1:10, co_odd, 1)

    ## Lighter tint in the middle
    co_light = g2r.colors(10, min.tinge=0.50)
    co_light
    doLegend(1:11, co_light, 1)
}
\keyword{utilities}

