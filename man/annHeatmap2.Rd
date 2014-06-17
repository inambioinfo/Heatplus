\name{annHeatmap2}
\alias{annHeatmap2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Annotated heatmaps}
\description{
This function plots a data matrix as intensity heatmap, with optional dendrograms, annotation panels and clustering for both rows and columns. This is the actual working function called by numerous wrappers.
}
\usage{
annHeatmap2(x, dendrogram, annotation, cluster, labels, 
            scale = c("row", "col", "none"), breaks=256, col = g2r.colors,
            legend = FALSE)
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
    Arguments \code{scale}, \code{breaks}, \code{col} and \code{legend} control different aspects of the whole plot directly as described. Arguments \code{dendrogram}, \code{annotation}, \code{cluster} and \code{labels} control aspects that may differ for the rows and columns of the central heatmap and have a special structure: each is a named list with different entries controling e.g. the look of a dendrogram, the data for annotation etc. Additionally, they can contain two extra entries called simply \code{Row} and \code{Col}; these are again named lists that can contain all the same entries as the parent list. Entries specified directly in the list apply to both rows and columns; entries specified as part of \code{Row} or \code{Col} override these defaults for the rows or columns only. 
    
    Recognized parameters for argument \code{dendrogram}:
    \describe{
        \item{clustfun}{the clustering function for generating the dendrogram; defaults to \code{hclust} for rows and columns}
        \item{distfun}{a function that returns the pairwise distances between samples/features as an object of class \code{dist}; defaults to \code{dist} for rows and columns}
        \item{status}{a string that controls the display of the dendrogram: \code{yes} means use the dendrogram to re-order the rows/columns and display the dendrogram; \code{hidden} means re-rorder, but do not display; \code{no} means do not use the dendrogram at all.}
        \item{lwd}{the line width of the branches of the dendrogram; defaults to 3.}
        \item{dendro}{an override argument that allows to pass in a dendrogram directly, bypassing the \code{clustfun} and \code{distfun} mechanism; defaults to \code{NULL} (i.e. is not used)}
    }
    
    Recognized entries for argument \code{annotation}: 
    \describe{
        \item{data}{a data frame containing the annotation data; defaults to \code{NULL}, i.e. no annotation is displayed}
        \item{control}{a list of fine-tuning parameters that is passed directly to \code{picketPlot}; defaults to an empty list, i.e. the default settings in \code{picketPlot}}
        \item{asIs}{logical value indicating whether the annotation \code{data} needs to be pre-processed via \code{convAnnData} or not; defaults to \code{TRUE}}
        \item{inclRef}{logical value indicating whether to include all levels of factor variables in \code{data}, or whether to drop the reference level (i.e. the first level). Defaults to \code{TRUE}}
    }
    
    Recognized entries for argument \code{cluster}:
    \describe{
        \item{cuth}{the height at which to cut through the dendrogram to define groups of similar features/samples; defaults to \code{NULL}, i.e. no    cutting}
        \item{label}{labels for the clusters; defaults to \code{NULL}, i.e. no labels}
        \item{col}{colors for the different clusters; the colors are used for coloring both the sub-trees of the dendrogram and the corresponding area in the annotation plot (if there is one). This is either a vector of colors, or a palette function that takes a number and returns a vector of colors of the specified length; defaults to \code{RainbowPastel}}
        \item{grp}{an override argument that directly specifies group memberships for the features/samples, completely bypassing the whole \code{dendrogram} and \code{cuth} mechanism. This probably only works for \code{dendrogram$status="no"}.}
    }

    Recognized entries for argument \code{labels}:
    \describe{
        \item{cex}{size of the text for the labels; defaults to \code{NULL}, i.e. use a hard-coded default guess}
        \item{nrow}{amount of space available for the labels between the central heatmap and the dendrogram, expressed as lines of text; defaults to 3.}
        \item{side}{side at which to draw the labels, coded as integer between 1 and 4 in the usual way (1 = below the plot, continuing clockwise). A common default for rows and columns does not make sense: rows only work with 2 and 4, columns only with 1 and 3. Defaults try to make use of empty space, depending on the presence of a dendrogram.}
        \item{labels}{labels for rows and columns; defaults to \code{NULL}, i.e. using the row- and column names of \code{x}.}
    }
}
\value{
    An object of class \code{annHeatmap}. Use \code{plot} to display it graphically.
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
map1 = annHeatmap2(exprs(ex1), ann=list(Col=list(data=pData(ex1))),
                   cluster=list(Col=list(cuth=3000)))
plot(map1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

