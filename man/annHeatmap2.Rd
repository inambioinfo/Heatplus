\name{annHeatmap2}
\alias{annHeatmap2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Annotated heatmaps}
\description{
This function plots a data matrix as intensity heatmap, with optional dendrograms, annotation panels and clustering for both rows and columns. This is the actual working function called by numerous wrappers.
}
\usage{
annHeatmap2(x, dendrogram, annotation, cluster, labels, scale = c("row", "col", "none"), breaks=256, col = g2r.colors, legend = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the numerical matrix to be shown as heatmap
}    
  \item{dendrogram}{
a list that controls how row- and column diagrams are determined and displayed
}
  \item{annotation}{
a list that controls the data and the way it is shown in row- and column annotation panels 
}
  \item{cluster}{
a list that controls how many clusters are chosen, and how these clusters are labelled and colored
}
  \item{labels}{
a list that controls the row- and column labels, as well as their size and placement
}
  \item{scale}{
      a character string indicating how the matrix \code{x} is standardized (by row, by column or not at all). This affects only display, not dendrograms or clustering
}
\item{breaks}{
    specifies the interval breaks for displaying the data in \code{x}; either a vector of explicit interval breaks, or just the desired number of intervals. See \code{niceBreaks} for details.
}
  \item{col}{
      specifies a palette of colors for the heatmap intensities; either a vector of explicit color definitions (one less than breaks) or a palette function. See \code{breakColors}.
}
  \item{legend}{
      whether and where to draw a legend for the colors/intervals in the heatmap. If \code{TRUE}, a legend is placed in a position determined by the function to be suitable; alternatively, integer values 1-4 indicate the side where the legend is to be drawn; and \code{FALSE} indicates that no legend should be drawn.
}
}
\details{

}
\value{
    An object of class \code{annHeatmap}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\seealso{
    \code{\link{heatmapLayout}}, \code{\link{niceBreaks}}, \code{\link{breakColors}}, \code{\link{g2r.colors}}
}
\examples{
require(Biobase)
data(sample.ExpressionSet)
ex1 = sample.ExpressionSet[51:85,]
map1 = annHeatmap2(exprs(ex1), ann=list(Col=list(data=pData(ex1))), cluster=list(Col=list(cuth=3000)))
plot(map1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

