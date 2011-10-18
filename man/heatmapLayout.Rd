\name{heatmapLayout}
\alias{heatmapLayout}
\title{Generate a layout for an (annotated) heatmap}
\description{
    Generate a layout for an (annotated) heatmap. This function will generally not be called directly, but only via \code{annHeatmap2}.
}
\usage{
heatmapLayout(dendrogram, annotation, leg.side = NULL, show = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendrogram}{
      A list with named entries \code{Row} and \code{Col}. Each of these is a list with a named entry \code{status}. If the value of \code{status} is the string \code{"yes"}, space will be set aside for drawing a row- and/or column dendrogram. 
}
  \item{annotation}{
      A list with named entries \code{Row} and \code{Col}. Each of these is a list with a named entry \code{data}. If the value of \code{data} is not \code{NULL}, space will be set aside for a picket plot showing the row- and/or column annotation.
}
  \item{leg.side}{
      An integer indicating on where to reserve space for the legend: values 1-4 correspond to below, to the left, above and to the right, as in e.g. \code{axis}. For a value of \code{NULL}, the function provides a reasonable default where there is space left in the layout. For any other value, no space for a legend is put aside. 
}
  \item{show}{
      A logical value; if \code{TRUE}, the layout defined by the arguments is displayed graphically.
}
}
\details{
    Space for plots is reserved via the \code{layout} mechanism. The function starts with an empty maximum layout, fills in the plot, dendrograms, annotation plots and legend as required, and compresses the resulting layout by removing empty slots.
}
\value{
    A list with the following entries:
    \item{plot }{A matrix describing the plot layout; see \code{layout}}
    \item{width }{relative widths of plots (i.e. columns)}
    \item{height }{relative heights of plots (i.e. rows)}
    \item{legend.side }{side where to draw the legend}
}
\seealso{
    \code{\link{annHeatmap2}}, \code{\link{picketPlot}}, \code{\link{layout}}
}
\examples{
    ## Heatmap with column dendrogram, column annotation, default legend
    dnd = list(Row=list(status="no"), Col=list(status="yes"))
    ann = list(Row=list(data=NULL), Col=list(data=1))
    ## 1 = heatmap, 2=dendrogram, 3=annotation, 4=legend
    ll = heatmapLayout(dendrogram=dnd, annotation=ann, leg.side=NULL, show=TRUE)
    ll
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

