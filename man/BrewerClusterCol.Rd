\name{BrewerClusterCol}
\alias{BrewerClusterCol}
\title{Color scheme for clusters}
\description{
This function returns a color vector based on one of the qualitative paletters supported by \code{RColorBrewer}. This allows visually distinct coloring of clusters and ensures sure that adjacent clusters have different colors. 
}
\usage{
BrewerClusterCol(n, name="Pastel1")
}
\arguments{
  \item{n}{desired number of colors}
  \item{name}{name of the qualitative palette from which colors are taken, see \code{\link{brewer.pal.info}} }
}
\details{
This is just a wrapper for \code{\link{brewer.pal}} that checks that the specified palette is qualitative, and allows for an arbitrary number of colors: for less than three colors, it just returns the first and second colors of the palette; for more than \code{maxcolors} colors, it recycles the basic palette as often as required. This is ok, because the main point is to have different colors for neighboring clusters. 
}
\value{
A character vector of length \code{n} of hexadecimal color codes.
}
\seealso{\code{\link{brewer.pal}}}
\examples{
## A Color Wheel: default palette with maximum number of colors
pie(rep(1,9), col=BrewerClusterCol(9))

## Double the number of colors 
pie(rep(1,18), col=BrewerClusterCol(18))

## Only two clusters/colors
pie(rep(1,2), col=BrewerClusterCol(2))

## Different qualitative palette: stronger colors
pie(rep(1,12), col=BrewerClusterCol(12, "Paired"))

}
\keyword{color}% at least one, from doc/KEYWORDS

