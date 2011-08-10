\name{oldCutplot.dendrogram}
\alias{oldCutplot.dendrogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Subtrees of a Dendrogram in Different Colors}
\description{
Plot a dendrogram, cut the tree at a given height, and draw the resulting
subtrees in different colors (OLD version, to be deprecated)
}
\usage{
oldCutplot.dendrogram(x, h, cluscol, leaflab = "none", horiz = FALSE, lwd = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a dendrogram.}
  \item{h}{the height at which the dendrogram is cut.}
  \item{cluscol}{the colors used for the subtrees; defaults to \code{rainbow}.}
  \item{leaflab}{indicates how leaf labels are to be drawn< defaults to 'perpendicular'.}
  \item{horiz}{logical indicating whether to plot the dendrogram horizontally or vertically.}
  \item{lwd}{the line width used for the color subtrees.}
  \item{\dots}{arguments to \code{plot.dendrogram}.}
}
\details{
This routine makes use of the functions \code{plot.dendrogram} and \code{plotNode} in package \code{stats}.
}
\author{Alexander Ploner <Alexander.Ploner@ki.se>}
\seealso{\code{\link{as.dendrogram}}}
\examples{
data(swiss)
cc = as.dendrogram(hclust(dist(swiss)))
oldCutplot.dendrogram(cc, h=80)
}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{aplot}% at least one, from doc/KEYWORDS
\keyword{cluster}% at least one, from doc/KEYWORDS
