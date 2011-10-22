\name{g2r.colors}
\alias{g2r.colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
g2r.colors(n = 12, min.tinge = 0.33)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{min.tinge}{
%%     ~~Describe \code{min.tinge} here~~
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
function (n = 12, min.tinge = 0.33) 
{
    k <- trunc(n/2)
    if (2 * k == n) {
        g <- c(rev(seq(min.tinge, 1, length = k)), rep(0, k))
        r <- c(rep(0, k), seq(min.tinge, 1, length = k))
        colvec <- rgb(r, g, rep(0, 2 * k))
    }
    else {
        g <- c(rev(seq(min.tinge, 1, length = k)), rep(0, k + 
            1))
        r <- c(rep(0, k + 1), seq(min.tinge, 1, length = k))
        colvec <- rgb(r, g, rep(0, 2 * k + 1))
    }
    colvec
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
