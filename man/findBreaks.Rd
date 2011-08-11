\name{findBreaks}
\alias{findBreaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specify breaks for a range of values}
\description{
Specify breaks for a range of values
}
\usage{
findBreaks(x, n, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function(x, n, type)
{
    if (missing(type)) type = "length"
    if (missing(n))    n = nclass.Sturges(x)
    if (type=="length") {
        rx = range(x, na.rm=TRUE)    
        dx = diff(rx)
        if (dx==0) dx = abs(rx[1])
        brks = seq.int(from=rx[1]-dx/1000, to=rx[2]+x/1000, length.out=n+1)
    } else if (type=="count") {
        brks = quantile(x, (0:n)/n, na.rm=TRUE)
    } else stop("Invalid type argument")
    attr(brks, "type") = type
    brks
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
