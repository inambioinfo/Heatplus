\name{print.annHeatmap}
\alias{print.annHeatmap}
\title{Printing information about annotated heatmaps}
\description{
Printing method for annotated heatmaps
}
\usage{
    \method{print}{annHeatmap}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{an object of class \code{annHeatmap}
}
    \item{\dots}{extra arguments, currently ignored
}
}
\details{
A very simple printing method, displaying a minimum of information about dendrograms and annotation
}
\value{
    \code{x} is returned invisibly
}
\seealso{
    \code{\link{annHeatmap}}, \code{\link{annHeatmap2}}, \code{\link{plot.annHeatmap}}
}
\examples{
    set.seed(219)
    mat = matrix(rnorm(100), ncol=5)
    ann = data.frame(Class=c("A","A","B","A","B"))
    map1 = annHeatmap(mat, ann)
    map1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

