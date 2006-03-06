\name{RGBColVec}
\alias{RGBColVec}
\alias{RainbowPastel}
\title{Alternative color schemes}
\description{
\code{RGBColVec} returns a vector of colors that is equally spaced from red through black to green, suitable for heatmaps.

\code{RainbowPastel} returns a vector of colors like \code{rainbow}, but more pastelly.
}
\usage{
RGBColVec(nrgcols = 12)

RainbowPastel(n, blanche = 200, ...)
}
\arguments{
  \item{nrgcols, n}{desired number of colors}
  \item{blanche}{the amount of whiteness added; value between 0 and 255}
  \item{\dots}{extra arguments to \code{rainbow}}
}
\value{
A character vector of length \code{nrgcols} or \code{n} giving the RGB codes for the colors.}
\author{\code{RGBColVec} is based on function \code{rgcolors.func} in package \code{sma} by Sandrine Dudoit and  Jane Fridlyand.

\code{RGBColVec} as documented and \code{RainbowPastel} by Alexander Ploner 
}
\seealso{\code{\link{heat.colors}}}
\examples{
 # A Color Wheel
pie(rep(1,12), col=RGBColVec(12))

# A color wheel in the original rainbow
pie(rep(1,6), col=rainbow(6))

# Pastel
pie(rep(1,6), col=RainbowPastel(6))

# Less whiteness
pie(rep(1,6), col=RainbowPastel(6, blanche=127))

# More steps require less whiteness
pie(rep(1,12), col=RainbowPastel(12, blanche=60))

# Test your screen & eyes: any differences?
pie(rep(1,12), col=RainbowPastel(12, blanche=80))
}
\keyword{color}% at least one, from doc/KEYWORDS

