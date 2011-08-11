\name{convAnnData}
\alias{convAnnData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converting data frames for display as annotation}
\description{
Converting data frames for display as annotation
}
\usage{
convAnnData(x, nval.fac = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{nval.fac}{
%%     ~~Describe \code{nval.fac} here~~
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
function(x, nval.fac=3)
{
    if (is.null(x)) return(NULL)
    x = as.data.frame(x)
    if (!is.null(nval.fac) & nval.fac>0) doConv = TRUE
    vv = colnames(x)
    for (v in vv) {
        if (is.logical(x[,v])) x[,v] = factor(as.numerical(x[,v]))
        if (doConv & length(unique(x[is.finite(x[,v]),v])) <= nval.fac) x[,v] = factor(x[,v])    
    }
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
