\name{cutplot.dendrogram}
\alias{cutplot.dendrogram}
\title{Plot a dendrogram in different colors}
\description{
Plot a dendrogram, optionally with different colors for all sub-dendrograms below a certain cutting height.
}
\usage{
cutplot.dendrogram(x, h, cluscol, leaflab = "none", horiz = FALSE, lwd = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{cluscol}{
%%     ~~Describe \code{cluscol} here~~
}
  \item{leaflab}{
%%     ~~Describe \code{leaflab} here~~
}
  \item{horiz}{
%%     ~~Describe \code{horiz} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
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
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, h, cluscol, leaflab= "none", horiz=FALSE, lwd=3, ...)
#
# Name: cutplot.dendrogram
# Desc: takes a dendrogram as described in library(mva), cuts it at level h,
#       and plots the dendrogram with the resulting subtrees in different 
#       colors
# Auth: obviously based on plot.dendrogram in library(mva)
#       modifications by Alexander.Ploner@meb.ki.se  211203
#
# Chng: 050204 AP 
#       changed environment(plot.hclust) to environment(as.dendrogram) to
#       make it work with R 1.8.1
#       250304 AP added RainbowPastel() to make it consistent with picketplot
#       030306 AP slightly more elegant access of plotNode
#       220710 AP also for horizontal plots
#
{
    ## If there is no cutting, we plot and leave
    if (missing(h) | is.null(h)) {
        return(plot(x, leaflab=leaflab, horiz=horiz, ...))
    }
    
    ## Some param processing
    if (missing(cluscol) | is.null(cluscol)) cluscol = RainbowPastel
    
    # Not nice, but necessary
    pn  = stats:::plotNode
    
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
   
    x = cut(x, h)
    plot(x[[1]], leaflab="none", horiz=horiz, ...)
    
    x = x[[2]]
    K = length(x)
    if (is.function(cluscol)) {
       cluscol = cluscol(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=horiz)
        x1 = x2 + 1
   }
   
    
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
