\name{picketPlot}
\alias{picketPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Display a data frame of annotation information}
\description{
Display a data frame of annotation information
}
\usage{
picketPlot(x, grp = NULL, grpcol, grplabel = NULL, add = FALSE, horizontal = TRUE, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{grp}{
%%     ~~Describe \code{grp} here~~
}
  \item{grpcol}{
%%     ~~Describe \code{grpcol} here~~
}
  \item{grplabel}{
%%     ~~Describe \code{grplabel} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
}
  \item{horizontal}{
%%     ~~Describe \code{horizontal} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
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
function (x, grp=NULL, grpcol, grplabel=NULL, add=FALSE, horizontal=TRUE, control=list()) 
#
# Name: picketPlot (looks like a picket fence with holes, and sounds like the
#                   pocketplot in geostatistics)
# Desc: visualizes a pattern of 0/1/NAs by using bars, great for annotating a 
#       heatmap
# Auth: Alexander.Ploner@meb.ki.se  181203
#
# Chng: 221203 AP loess() with degree < 2
#       260104 AP 
#       - made loess() optional
#       - better use of space if no covariate
#       030304 AP
#       - added RainbowPastel() as default colors
#       - check grplabel before passing it to axis
#       2010-07-08 AP
#       - complete re-write
#       2010-08-28
#       - re-arranged code for vertical/horizontal drawing
{    
    # deal with the setup
    cc = list(boxw=1, boxh=4, hbuff=0.1, vbuff=0.1, span=1/3, nacol=gray(0.85), 
              degree=1, cex.label=1.5, numfac=2)
    cc[names(control)] = control
    
    # Count variables, panels, types
    nsamp  = nrow(x)
    nvar   = ncol(x)
    npanel = sapply(x, function(x) if (is.factor(x)) nlevels(x) else 1)
    bpanel = sapply(x, is.factor)
    
    # Compute panel heights, widths
    panelw = nsamp*(cc$boxw+2*cc$hbuff)
    panelh = cc$boxh+2*cc$vbuff
    totalh = sum(npanel * panelh * ifelse(bpanel, 1, cc$numfac))
    LL = cbind(0, 0)
    UR = cbind(panelw, totalh)
    
    # Set up the x-values for a single panel
    xbase = seq(cc$hbuff, by=cc$boxw+2*cc$hbuff, length=nsamp)
    xcent = xbase + cc$boxw/2
    
    # if we get a cluster variable, we have to set differently colored 
    # backgrounds; this assumes that the grp variable is sorted in the 
    # way it appears on the plot
    if (!is.null(grp)) {
        grp = as.integer(factor(grp, levels=unique(grp)))
        tt = table(grp)
        gg = length(tt)
        grpcoord = c(0,cumsum(tt/sum(tt))*panelw)
        grp0 = cbind(grpcoord[1:gg], rep(0, gg))
        grp1 = cbind(grpcoord[2:(gg+1)], rep(totalh, gg))
        if (missing(grpcol)) {
            grpcol=RainbowPastel(gg)
        }
    }

    # Loop over vars and fill in the panels
    panels = list()
    voff = 0
    for (i in 1:nvar) {
        if (bpanel[i]) {
            ## Coordinates
            x0 = rep(xbase, npanel[i])
            x1 = x0+cc$boxw
            y0 = voff + rep(seq(cc$vbuff, by=cc$boxh+2*cc$vbuff, length=npanel[i]), rep(nsamp, npanel[i]))
            y1 = y0 + cc$boxh
            ## Convert factor to matrix & set fill
            naAction = attr(na.exclude(x[,i, drop=FALSE]), "na.action")
            bindata = naresid(naAction, model.matrix(~x[,i]-1))
            fill = ifelse(bindata==1, "black", "transparent")
            fill[is.na(fill)] = cc$nacol
            label = paste(colnames(x)[i], "=", levels(x[,i]))
            labcc = if (!is.null(label)) sort(unique((y0+y1)/2)) else NULL 
            panels[[i]] = list(ll=cbind(x0, y0), ur=cbind(x1, y1), fill=fill, label=label, labcc=labcc)
            voff = voff + panelh*npanel[i]
        } else {
            xv = x[,i]
            rr = range(xv, na.rm=TRUE)
            yval = voff + cc$vbuff*cc$numfac + ((xv - rr[1])/(rr[2] - rr[1]))*cc$boxh*cc$numfac
            if ((cc$degree>0) & (cc$span>0)){
                yy = predict(loess(yval~xcent, span=cc$span, degree=cc$degree))
            }
            label = colnames(x)[i]
            labcc = if (!is.null(label)) mean(range(yval, na.rm=TRUE)) else NULL
            axlab = pretty(range(xv, na.rm=TRUE))
            axcc  = voff + cc$vbuff*cc$numfac + ((axlab - rr[1])/(rr[2] - rr[1]))*cc$boxh*cc$numfac
            panels[[i]] = list(raw=cbind(xcent, yval), smo=cbind(xcent, yy), label=label, labcc=labcc, axlab=axlab, axcc=axcc)
            voff = voff + panelh*cc$numfac
        }
    }
    
    # if grplabels are given, we add another horizontal axis to the 
    # last plot (independent of whether it is binvar or contvar)
    if (!is.null(grp) & !is.null(grplabel)) {
        mids = (grpcoord[1:gg] + grpcoord[2:(gg+1)])/2
        # Is the grplabel ok?
        labelnum = length(grplabel)
        if (labelnum < gg) {
            warning("more groups than labels (filling up with blanks)")
            grplabel = c(grplabel, rep(" ", gg-labelnum))
        } else if (gg < labelnum) {
            warning("more labels than groups (ignoring the extras)")
            grplabel = grplabel[1:gg]
        }
    }
            
    ## Switch coordinates, if you have to
    h2v = function(cc) cbind(cc[,2]-totalh, cc[,1])
    if (horizontal) {
        grpaxis = 1
        labaxis = 2
        covaxis = 4
        las = 1        
    } else {
        grpaxis = 4
        labaxis = 3
        covaxis = 1
        las = 3
        ## Rotate
        LL = h2v(LL)
        UR = h2v(UR)
        if (!is.null(grp)) {
            grp0 = h2v(grp0)
            grp1 = h2v(grp1)
        }
        for (i in 1:nvar) {
            panels[[i]][[1]] = h2v(panels[[i]][[1]])
            panels[[i]][[2]] = h2v(panels[[i]][[2]])
            panels[[i]]$labcc = panels[[i]]$labcc - totalh 
            panels[[i]]$axcc  = panels[[i]]$axcc - totalh
        }
    }
    
    # Set up the plot
    plot(rbind(LL, UR), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
    # Add the colored rectangles, if required
    if (!is.null(grp)) {
        rect(grp0[,1], grp0[,2], grp1[,1], grp1[,2], col=grpcol, border="transparent")
    }
    # Loop over vars and fill in the panels
    for (i in 1:nvar) {
        if (bpanel[i]) {
            ## Do the rectangles
            with(panels[[i]], rect(ll[,1], ll[,2], ur[,1], ur[,2], col=fill, border="transparent") )
        } else {
            with(panels[[i]], points(raw[,1], raw[,2]))
            if ((cc$degree>0) & (cc$span>0)){
                with(panels[[i]], lines(smo[,1], smo[,2]))
            }
            with(panels[[i]], axis(covaxis, at=axcc, label=axlab))
        }
        ## Name panel (regardless of type)
        if (!is.null(panels[[i]]$label)) {
            axis(labaxis, at=panels[[i]]$labcc, label=panels[[i]]$label, las=las, tick=FALSE, font=2, col=par("bg"), col.axis=par("fg"))
        }
    }
    # if grplabels are given, we add another horizontal axis to the 
    # last plot (independent of whether it is binvar or contvar)
    if (!is.null(grp) & !is.null(grplabel)) {
        axis(grpaxis, grpcoord, label=FALSE, tcl=-1.5)
        axis(grpaxis, mids, label=grplabel, font=2, cex.axis=cc$cex.label, tick=FALSE)
    }                                

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
