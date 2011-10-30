\name{regHeatmap}
\alias{regHeatmap}
\alias{regHeatmap.default}
\title{Regular heatmaps with a legend}
\description{
    Creating regular heatmaps, without annotation, but allowing for a legend
}
\usage{
    regHeatmap(x, ...)

    \method{regHeatmap}{default}(x, dendrogram = list(clustfun = hclust, distfun = dist, status = "yes"), labels = NULL, legend = TRUE, ...)
}
\arguments{
  \item{x}{a numerical matrix
}
  \item{dendrogram}{a list controlling the options for row- and column dendrogram, see \code{annHeatmap2}
}
  \item{labels}{a list controlling the row- and column labels as well as their location and size, see \code{annHeatmap2}
}
  \item{legend}{either a logical value, indicating whether to draw a legend at the default location determined by the function, or one of the sides of the plot (1-4), see \code{annHeatmap2}
}
  \item{\dots}{extra options passed to \code{annHeatmap2}
}
}
\details{
    A gelded wrapper for \code{annHeatmap2} that allows for heatmaps without annotation or clustering on the dendrograms, but still offer some control over dendrograms, labels and legend.
    
    These functions generate an object representing the heatmap; in order to produce graphical output, you have to invoke the \code{plot} method, see Examples.
}
\value{
    An object of class \code{annHeatmap}
}
\seealso{
    \code{\link{annHeatmap}}, \code{\link{annHeatmap2}}, \code{\link{plot.annHeatmap}}
}
\examples{
    
    ## Default
    set.seed(219)
    mat = matrix(rnorm(100), ncol=5)
    map1 = regHeatmap(mat)
    plot(map1)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

