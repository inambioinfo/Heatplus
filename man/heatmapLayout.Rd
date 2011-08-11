\name{heatmapLayout}
\alias{heatmapLayout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a layout for an (annotated) heatmap}
\description{
Generate a layout for an (annotated) heatmap
}
\usage{
heatmapLayout(dendrogram, annotation, leg.side = NULL, show = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
}
  \item{annotation}{
%%     ~~Describe \code{annotation} here~~
}
  \item{leg.side}{
%%     ~~Describe \code{leg.side} here~~
}
  \item{show}{
%%     ~~Describe \code{show} here~~
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
function(dendrogram, annotation, leg.side=NULL, show=FALSE)
{
    ## Start: maximum matrix, 5 x 5, all zero
    ## Names for nice display post ante
    ll = matrix(0, nrow=5, ncol=5)
    ll.width = ll.height = rep(0, 5)
    cnt = 1
    rownames(ll) = c("leg3", "colDendro","image", "colAnn", "leg1")
    colnames(ll) = c("leg2", "rowDendro","image", "rowAnn", "leg4")
 
    ## The main plot
    ll[3,3] = cnt
    ll.width[3] = ll.height[3] = 5
    cnt = cnt+1
    ## The column dendrogram
    if (dendrogram$Col$status=="yes") {
        ll[2, 3] = 2
        ll.width[3]  = 5
        ll.height[2] = 2
        cnt = cnt+1
    }
    ## The row dendrogram
    if (dendrogram$Row$status=="yes") {
        ll[3, 2] = cnt
        ll.width[2]  = 2
        ll.height[3] = 5
        cnt = cnt+1
    }
    # Column annotation
    if (!is.null(annotation$Col$data)) {
        ll[4, 3] = cnt
        ll.width[3] = 5
        ll.height[4] = 2
        cnt = cnt+1
    }
    ## Row annotation
    if (!is.null(annotation$Row$data)) {
        ll[3, 4] = cnt
        ll.width[4]  = 2
        ll.height[3] = 5
        cnt = cnt+1
    }
    ## Legend: if no pref specified, go for empty, if possible
    if (is.null(leg.side)) {
        if (dendrogram$Row$status != "yes") {
            leg.side = 2
        } else if (is.null(annotation$Row$data)) {
            leg.side = 4
        } else if (is.null(annotation$Col$data)) {
            leg.side = 1
        } else if (dendrogram$Col$status != "yes") {
            leg.side = 3
        } else {
            leg.side = 4
        }
    }
    ## Add the legend space
    if (leg.side==1) {
        ll[5,3] = cnt
        ll.width[3] = 5
        ll.height[5] = 1
    } else if (leg.side==2) {
        ll[3,1] = cnt
        ll.width[1] = 1
        ll.height[3] = 5
    } else if (leg.side==3) {
        ll[1,3] = cnt
        ll.width[3]  = 5
        ll.height[1] = 1
    } else if (leg.side==4) {
        ll[3,5] = cnt
        ll.width[5]  = 1
        ll.height[3] = 5
    }
    
    ## Compress
    ndx = rowSums(ll)!=0
    ll  = ll[ndx, ]
    ll.height = ll.height[ndx]
    ndx = colSums(ll)!=0
    ll  = ll[, ndx]
    ll.width = ll.width[ndx]
    ## Do it - show it
    if (show) {
        layout(ll, width=ll.width, height=ll.height, respect=TRUE)            
        layout.show(max(ll))
    }
    return(list(plot=ll, width=ll.width, height=ll.height, legend.side=leg.side))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
