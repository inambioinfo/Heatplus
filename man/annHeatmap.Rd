\name{annHeatmap}
\alias{annHeatmap}
\alias{annHeatmap.default}
\alias{annHeatmap.ExpressionSet}
\title{Annotated heatmaps}
\description{
    Creating heatmaps with annotated columns
}
\usage{
    annHeatmap(x, ...)

    \method{annHeatmap}{default}(x, annotation, dendrogram = list(clustfun = hclust, distfun = dist, Col = list(status = "yes"), Row = list(status = "hidden")), cluster = NULL, labels = NULL, legend = TRUE, ...)

    \method{annHeatmap}{ExpressionSet}(x, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{either a numerical matrix with the data for the central heatmap (for the default method) or an object of class \code{ExpressionSet} 
}
    \item{annotation}{a data frame containing the annotation for the columns of \code{x}
}
    \item{dendrogram}{a list controlling the options for row- and column dendrogram, see \code{annHeatmap2}
}
\item{cluster}{a list controlling the options for clustering rows and columns of \code{x}, see \code{annHeatmap2}
}
  \item{labels}{a list controlling the row- and column labels as well as their location and size, see \code{annHeatmap2} 
}
\item{legend}{either a logical value, indicating whether to draw a legend at the default location determined by the function, or one of the sides of the plot (1-4), see \code{annHeatmap2}
}
  \item{\dots}{extra options passed to \code{annHeatmap2}
}
}
\details{
    These functions generate an object representing the heatmap; in order to produce graphical output, you have to invoke the \code{plot} method, see Examples.
}
\section{Warning}{These are currently simple convenience functions that allow quick plotting, but little control over the finer details. This may change in the future, but for now, if you want to do anything fancy, you should invoke \code{annHeatmap2} directly.}
\value{
    An object of class \code{annHeatmap}
}
\seealso{
    \code{\link{annHeatmap2}}, \code{\link{plot.annHeatmap}}
}
\examples{

    ## Default method
    set.seed(219)
    mat = matrix(rnorm(100), ncol=5)
    ann = data.frame(Class=c("A","A","B","A","B"))
    map1 = annHeatmap(mat, ann)
    plot(map1)
    
    ## Expression set
    require(Biobase)
    data(sample.ExpressionSet)
    map2 = annHeatmap(sample.ExpressionSet)
    plot(map2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

