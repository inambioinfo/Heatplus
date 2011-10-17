\name{modifyExistingList}
\alias{modifyExistingList}
\alias{extractArg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Override existing list entries}
\description{
Override existing list entries and extract arguments that are specified as named lists
}
\usage{
modifyExistingList(x, val)

extractArg(arglist, deflist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{val}{
%%     ~~Describe \code{val} here~~
}
  \item{arglist}{
%%     ~~Describe \code{arglist} here~~
}
  \item{deflist}{
%%     ~~Describe \code{deflist} here~~
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
function(x, val)
{
    if (is.null(x)) x = list()
    if (is.null(val)) val = list()
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    for (v in intersect(xnames, vnames)) {
        x[[v]] <- if (is.list(x[[v]]) && is.list(val[[v]])) 
            modifyExistingList(x[[v]], val[[v]])
        else val[[v]]
    }
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
