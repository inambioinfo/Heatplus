\name{niceBreaks}
\alias{niceBreaks}
\title{Get nice (symmetric) breaks for an interval}
\description{
Given a minimum and a maximum, this function returns a vector of equidistant breaks that covers this interval, and has a pretty interval length (1, 2, or 5 times a power of 10). If the interval contains zero, it will be one of the breaks, so that the intervals are arranged somewhat symmetrically around it.
}
\usage{
niceBreaks(xr, breaks)
}
\arguments{
    \item{xr}{the range to be covered, as \code{c(min, max)}}
    \item{breaks}{either the desired number of breaks, or a pre-specified vector of breaks}
}
\details{
The number of desired breaks is honored as far as possible, which is not actually that often in practice. However, major deviations of three or more are reasonably rare.

The functiona allows the specification of a set of breaks instead of the desired number of breaks, somewhat like in \code{cut}. However, if the length of \code{breaks} is greater than one, the function just returns \code{breaks}.
}
\value{A vector of pretty breaks covering the specified interval, more or less of the desired length. 
}
\seealso{
    \code{\link{pretty}}
}
\examples{
    
    ## Niceness overrules specified number
    niceBreaks(c(-1,1), 5)
    niceBreaks(c(-1,1), 6)

    ## Zero appears always as break
    niceBreaks(c(-2.75, 1.12), 8)

    ## Not invariant to translation (of course)
    niceBreaks(3.27 + c(-2.75, 1.12), 8)

    
}
\keyword{utilities}

