\name{cutree.dendrogram}
\alias{cutree.dendrogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cut a dendrogram at a specified height}
\description{
Cut a dendrogram at a specified height
}
\usage{
cutree.dendrogram(x, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
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
function(x, h)
{
    cutx = cut(x, h)
    cutl = lapply(cutx$lower, getLeaves)
    nclus = sapply(cutl, length)
    ret   = rep(1:length(nclus), nclus)
    ret[order.dendrogram(x)] = ret
    names(ret)[order.dendrogram(x)] = unlist(cutl)
    ret = as.integer(factor(ret, levels=unique(ret))) # recode for order of clus
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
