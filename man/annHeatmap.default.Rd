\name{annHeatmap.default}
\alias{annHeatmap.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default method for annotated heatmaps}
\description{
Default method for annotated heatmaps
}
\usage{
    \method{annHeatmap}{default}(x, annotation, dendrogram = list(clustfun = hclust, distfun = dist, Col = list(status = "yes"), Row = list(status = "hidden")), cluster = NULL, labels = NULL, legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{annotation}{
%%     ~~Describe \code{annotation} here~~
}
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{legend}{
%%     ~~Describe \code{legend} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, annotation, dendrogram=list(clustfun=hclust, distfun=dist, Col=list(status="yes"), Row=list(status="hidden")), cluster=NULL, labels=NULL, legend=TRUE, ...)
{
    ret = annHeatmap2(x, dendrogram=dendrogram, annotation=list(Col=list(data=annotation, fun=picketPlot)), cluster=cluster,  labels=labels, legend=TRUE, ...)
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
