\name{print.annHeatmap}
\alias{print.annHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Printing method for annotated heatmaos}
\description{
Printing method for annotated heatmaos
}
\usage{
    \method{print}{annHeatmap}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function(x, ...)
{
    cat("annotated Heatmap\n\n")
    cat("Rows: "); show(x$dendrogram$Row$dendro)
    cat("\t", if (is.null(x$annotation$Row$data)) 0 else ncol(x$annotation$Row$data), " annotation variable(s)\n")
    cat("Cols: "); show(x$dendrogram$Col$dendro)
    cat("\t", if (is.null(x$annotation$Col$data)) 0 else ncol(x$annotation$Col$data), " annotation variable(s)\n")
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
