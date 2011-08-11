\name{annHeatmap2}
\alias{annHeatmap2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Working function for generating (annotated) heatmaps}
\description{
Working function for generating (annotated) heatmaps
}
\usage{
annHeatmap2(x, dendrogram, annotation, cluster, labels, scale = c("row", "col", "none"), col = heat.colors(12), equalize = FALSE, trim, legend = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
}
  \item{annotation}{
%%     ~~Describe \code{annotation} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{equalize}{
%%     ~~Describe \code{equalize} here~~
}
  \item{trim}{
%%     ~~Describe \code{trim} here~~
}
  \item{legend}{
%%     ~~Describe \code{legend} here~~
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
function(x, dendrogram, annotation, cluster, labels, scale=c("row", "col", "none"), col=heat.colors(12), equalize=FALSE, trim, legend=FALSE)
#
# Name: annHeatmap2
# Desc: a (possibly doubly) annotated heatmap
# Auth: Alexander.Ploner@ki.se 2010-07-12
#
# Chng: 
#
{
    ## Process arguments
    if (!is.matrix(x) | !is.numeric(x)) stop("x must be a numeric matrix")
    nc = ncol(x); nr = nrow(x)
    if (nc < 2 | nr < 2) stop("x must have at least two rows/columns")
    
    ## Process the different lists: dendrogram, cluster, annotation
    ## See lattice:::xyplot.formula, modifyLists, lattice:::construct.scales
    def = list(clustfun=hclust, distfun=dist, status="yes", dendro=NULL)
    dendrogram = extractArg(dendrogram, def)
    def = list(data=NULL, fun=picketPlot)
    annotation = extractArg(annotation, def)
    def = list(cuth=NULL, grp=NULL, label=NULL, col=RainbowPastel)
    cluster = extractArg(cluster, def)
    def = list(cex=NULL, nrow=3, side=NULL, labels=NULL)
    labels = extractArg(labels, def)
    
    ## Check values for the different lists

    
    ## Generate the layout
    layout = heatmapLayout(dendrogram, annotation)
    
    ## Copy the data for display, scale as required
    x2 = x
    scale = match.arg(scale)
    if (scale == "row") {
        x2 = sweep(x2, 1, rowMeans(x, na.rm = TRUE))
        sd = apply(x2, 1, sd, na.rm = TRUE)
        x2 = sweep(x2, 1, sd, "/")
    }
    else if (scale == "column") {
        x2 = sweep(x2, 2, colMeans(x, na.rm = TRUE))
        sd = apply(x2, 2, sd, na.rm = TRUE)
        x2 = sweep(x2, 2, sd, "/")
    }
    # More display magic - scae or trim, but generally not both
    if (!missing(trim)) {
        trim = min(trim[1], 1-trim[1])
        lo = quantile(x2, trim, na.rm=TRUE) 
        hi = quantile(x2, 1-trim, na.rm=TRUE)
        x2[x2<lo] = lo
        x2[x2>hi] = hi
    }
    # slightly dirty too - should not be used together with trim
    if (equalize) {
        att = attributes(x2)
        x2 = rank(x2, na.last=if (na.rm) TRUE else NA)
        attributes(x2) = att
    }
    
    ## Generate the dendrograms, if required; re-indexes in any cases
    ## We could put some sanity checks on the dendrograms in the else-branches
    dendrogram$Row = within(dendrogram$Row, 
        if (!inherits(dendro, "dendrogram")) {
            dendro = clustfun(distfun(x))
            dendro = reorder(as.dendrogram(dendro), rowMeans(x, na.rm=TRUE))
        }
    )
    dendrogram$Col = within(dendrogram$Col, 
        if (!inherits(dendro, "dendrogram")) {
            dendro = clustfun(distfun(t(x)))
            dendro = reorder(as.dendrogram(dendro), colMeans(x, na.rm=TRUE))
        }
    )
    ## Reorder the display data to agree with the dendrograms, if required
    rowInd = with(dendrogram$Row, if (status!="no") order.dendrogram(dendro) else 1:nr)
    colInd = with(dendrogram$Col, if (status!="no") order.dendrogram(dendro) else 1:nc)
    x2 = x2[rowInd, colInd]
    
    ## Set the defaults for the sample/variable labels
    labels$Row = within(labels$Row, {
        if (is.null(cex)) cex = 0.2 + 1/log10(nr)
        if (is.null(side)) side = if (is.null(annotation$Row$data)) 4 else 2
        if (is.null(labels)) labels = rownames(x2)
    })
    labels$Col = within(labels$Col, {
        if (is.null(cex)) cex = 0.2 + 1/log10(nc)
        if (is.null(side)) side = if (is.null(annotation$Col$data)) 1 else 3
        if (is.null(labels)) labels = colnames(x2)
    })
    
    ## Generate the clustering, if required (cut, or resort the cluster var)
    ## FIXME: does not deal with pre-defined grp form outside
    cluster$Row = within(cluster$Row, 
        if (!is.null(cuth) && (cuth > 0)) {
            grp = cutree.dendrogram(dendrogram$Row$dendro, cuth)[rowInd]
        })
    cluster$Col = within(cluster$Col, 
        if (!is.null(cuth) && (cuth > 0)) {
            grp = cutree.dendrogram(dendrogram$Col$dendro, cuth)[colInd]
        })

    ## Process the annotation data frames (factor/numeric, re-sort?)
    annotation$Row$data = convAnnData(annotation$Row$data)
    annotation$Col$data = convAnnData(annotation$Col$data)
        
    ## Generate the new object
    
    ## print, return invisibly
    ret = list(data=list(x=x, x2=x2, rowInd=rowInd, colInd=colInd, col=col), dendrogram=dendrogram, cluster=cluster, annotation=annotation, labels=labels, layout=layout, legend=legend)
    class(ret) = "annHeatmap"
    ret

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
