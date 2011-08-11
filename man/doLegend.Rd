\name{doLegend}
\alias{doLegend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add a legend to the heatmap}
\description{
Add a legend to the heatmap
}
\usage{
doLegend(r, col, side)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{side}{
%%     ~~Describe \code{side} here~~
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
function(r, col, side)
{
    nc = length(col)
    zval = seq(r[1], r[2], length=nc)
    z  = matrix(zval, ncol=1)
    if (side \%in\% c(1,3)) {
        image(x=zval, y=1, z=z, xaxt="n", yaxt="n", col=col)
    } else {
        image(x=1, y=zval, z=t(z), xaxt="n", yaxt="n", col=col)
    }        
    axis(side, las=1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
