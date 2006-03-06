\name{heatmap_2}
\alias{heatmap_2}
\title{Display Data as Heatmap}
\description{
This function displays an expression data matrix as a heatmap. It is based on
an old version of \code{heatmap} in the \code{stats} package, but offers more
flexibility (e.g. skipping dendrograms, skipping row/column labelling, adding a
legend).
}
\usage{
heatmap_2(x, Rowv, Colv, distfun = dist, hclustfun = hclust, add.expr, 
          scale = c("row", "column", "none"), na.rm = TRUE, 
		  do.dendro = c(TRUE, TRUE), legend = 0, legfrac = 8, 
		  col = heat.colors(12), trim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the numerical data matrix to be displayed.}
  \item{Rowv}{either a dendrogram or a vector of reordering indexes for the rows.}
  \item{Colv}{either a dendrogram or a vector of reordering indexes for the columns.}
  \item{distfun}{function to compute the distances between rows and columns. 
                 Defaults to \code{dist}.}
  \item{hclustfun}{function used to cluster rows and columns. Defaults to \code{hclust}.}
  \item{add.expr}{Expression to be evaluated after the call to \code{image}. See Details.}
  \item{scale}{indicates whether values should be scaled by either by row, column, or not at all. Defaults to \code{row}.}
  \item{na.rm}{logical indicating whther to remove NAs.}
  \item{do.dendro}{logical vector of length two, indicating (in this order) whether to draw the row and column dendrograms.}
  \item{legend}{integer between 1 and 4, indicating on which side of the plot the legend should be drawn, as in \code{mtext}.}
  \item{legfrac}{fraction of the plot that is taken up by the legend; larger values correspond to smaller legends.}
  \item{col}{the color scheme for \code{image}. The default sucks.}
  \item{trim}{Percentage of values to be trimmed. This helps to keep an informative color scale, see Details.}
  \item{\dots}{extra arguments to \code{image}.}
}
\details{
With all parameters at their default, this gives the same result as a very old  version of \code{heatmap} that was the base for the modifications. All parameters of the same name have the same function as in \code{heatmap}, though \code{add.expr}, which can be used for adding graphical elements after the call to \code{image}, will probably not produce useful results. Note also that row- and column labels are optional, i.e. if the corresponding \code{dimname} of \code{x} is \code{NULL}, no labels are displayed.

Setting \code{trim} to a number between 0 and 1 uses equidistant classes between the (\code{trim})- and (1-\code{trim})-quantile, and lumps the values below and above this range into separate open-ended classes. If the data comes from a 
heavy-tailed distribution, this can save the display from putting too many values into to few classes.
}
\value{
Same as \code{heatmap} with \code{keep.dendro=FALSE}: an invisible list giving the reordered indices of the row- and column-elements as elements \code{rowInd} and \code{colInd}.}
\author{Original by Andy Liaw, with revisions by Robert Gentleman and Martin Maechler.

Alexander Ploner for this version.
}
\seealso{\code{\link{heatmap}}, \code{\link{hclust}}, \code{\link{heatmap_plus}}}
\examples{
# create data
mm = matrix(rnorm(1000, m=1), 100,10)
mm = cbind(mm, matrix(rnorm(2000), 100, 20))
mm = cbind(mm, matrix(rnorm(1500, m=-1), 100, 15))
mm2 = matrix(rnorm(450), 30, 15)
mm2 = cbind(mm2, matrix(rnorm(900,m=1.5), 30,30))
mm=rbind(mm, mm2)
colnames(mm) = paste("Sample", 1:45)
rownames(mm) = paste("Gene", 1:130)

# similar to base heatmap
heatmap_2(mm)

# remove column dendrogram
heatmap_2(mm, do.dendro=c(TRUE, FALSE))

# add a legend under the plot
heatmap_2(mm, legend=1)
# make it smaller
heatmap_2(mm, legend=1, legfrac=10)
# ... on the left side
heatmap_2(mm, legend=2, legfrac=10)

# remove the column labels by removing the column names
colnames(mm)=NULL
heatmap_2(mm, legend=1, legfrac=10)

# truncate the data drastically
heatmap_2(mm, legend=1, legfrac=10, trim=0.1)
}
\keyword{hplot}% at least one, from doc/KEYWORDS

