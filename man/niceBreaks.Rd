\name{niceBreaks}
\alias{niceBreaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
niceBreaks(xr, breaks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xr}{
%%     ~~Describe \code{xr} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xr, breaks) 
{
    if (length(breaks) > 1) {
        return(breaks)
    }
    if ((xr[1] < 0) & (xr[2] > 0)) {
        xminAbs = abs(xr[1])
        xmax = xr[2]
        nneg = max(round(breaks * xminAbs/(xmax + xminAbs)), 
            1)
        npos = max(round(breaks * xmax/(xmax + xminAbs)), 1)
        nbr = pretty(c(xr[1], 0), nneg)
        pbr = pretty(c(0, xr[2]), npos)
        diff = ((nbr[2] - nbr[1]) + (pbr[2] - pbr[1]))/2
        nbr = diff * ((xr[1]\%/\%diff):0)
        pbr = diff * (1:(xr[2]\%/\%diff + 1))
        breaks = c(nbr, pbr)
    }
    else {
        breaks = pretty(xr, breaks)
    }
    breaks
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
