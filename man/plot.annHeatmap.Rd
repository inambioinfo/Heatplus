\name{plot.annHeatmap}
\alias{plot.annHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting method for annotated heatmaps}
\description{
Plotting method for annotated heatmaps
}
\usage{
plot.annHeatmap(x, ...)
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
    ## Set up the layout
    with(x$layout, layout(plot, width, height, respect=TRUE))
    
    ## Plot the central image, making space for labels, if required
    nc = ncol(x$data$x2); nr = nrow(x$data$x2)
    doRlab = !is.null(x$labels$Row$labels) 
    doClab = !is.null(x$labels$Col$labels)
    mmar = c(1, 0, 0, 2)
    if (doRlab) mmar[x$labels$Row$side] = x$labels$Row$nrow
    if (doClab) mmar[x$labels$Col$side] = x$labels$Col$nrow
    with(x$data, {
        par(mar=mmar)
        image(1:nc, 1:nr, t(x2), axes = FALSE, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5), xlab = "", ylab = "", col=col, ...)    
    })
    with (x$labels, {
        if (doRlab) axis(Row$side, 1:nr, las = 2, line = -0.5, tick = 0, labels = Row$labels, cex.axis = Row$cex)
        if (doClab) axis(Col$side, 1:nc, las = 2, line = -0.5, tick = 0, labels = Col$labels, cex.axis = Col$cex)
    })

    ## Plot the column/row dendrograms, as required
    with(x$dendrogram$Col,
        if (status=="yes") {
            par(mar=c(0, mmar[2], 3, mmar[4]))
            cutplot.dendrogram(dendro, h=x$cluster$Col$cuth, cluscol=x$cluster$Col$col, horiz=FALSE, axes = FALSE, xaxs = "i", leaflab = "none")
        })
    with(x$dendrogram$Row,
        if (status=="yes") {
            par(mar=c(mmar[1], 3, mmar[3], 0))
            cutplot.dendrogram(dendro, h=x$cluster$Row$cuth, cluscol=x$cluster$Row$col, horiz=TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
        })

    ## Plot the column/row annotation data, as required
    if (!is.null(x$annotation$Col$data)) {
        par(mar=c(1, mmar[2], 0, mmar[4]), xaxs="i", yaxs="i")
        picketPlot(x$annotation$Col$data[x$data$colInd,], grp=x$cluster$Col$grp, add=TRUE)
    }
    if (!is.null(x$annotation$Row$data)) {
        par(mar=c(mmar[1], 0, mmar[3], 1), xaxs="i", yaxs="i")
        picketPlot(x$annotation$Row$data[x$data$rowInd,], grp=x$cluster$Row$grp, add=TRUE, horizontal=FALSE)
    }

    ## Plot a legend, as required
    if (x$legend) {
        if (x$layout$legend.side \%in\% c(1,3)) {
            par(mar=c(2, mmar[2], 2, mmar[4]))
        } else {
            par(mar=c(mmar[1], 2, mmar[3], 2))            
        }        
        doLegend(range(x$data$x, na.rm=TRUE), x$layout$legend.side, col=x$data$col)
    }    
    
    invisible(x)    
    
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
