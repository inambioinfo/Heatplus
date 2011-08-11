###################################################
### chunk number 1: heatmapLayout_Def
###################################################
#line 10 "annHeatmap.Rnw"
heatmapLayout = function(dendrogram, annotation, leg.side=NULL, show=FALSE)
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


###################################################
### chunk number 2: modifyExistingList_Def
###################################################
#line 109 "annHeatmap.Rnw"
modifyExistingList = function(x, val)
{
    if (is.null(x)) x = list()
    if (is.null(val)) val = list()
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    for (v in intersect(xnames, vnames)) {
        x[[v]] <- if (is.list(x[[v]]) && is.list(val[[v]])) 
            modifyExistingList(x[[v]], val[[v]])
        else val[[v]]
    }
    x
}


###################################################
### chunk number 3: extractArg_Def
###################################################
#line 128 "annHeatmap.Rnw"
extractArg = function(arglist, deflist)
{
    if (missing(arglist)) arglist = NULL
    al2 = modifyExistingList(deflist, arglist)
    row = col = al2
    row = modifyExistingList(row, arglist[["Row"]])    
    col = modifyExistingList(col, arglist[["Col"]])    
    list(Row=row, Col=col)  
}


###################################################
### chunk number 4: picketPlot_Def
###################################################
#line 150 "annHeatmap.Rnw"
picketPlot = function (x, grp=NULL, grpcol, grplabel=NULL, add=FALSE, horizontal=TRUE, control=list()) 
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


###################################################
### chunk number 5: findBreaks_Def
###################################################
#line 331 "annHeatmap.Rnw"
findBreaks = function(x, n, type)
{
    if (missing(type)) type = "length"
    if (missing(n))    n = nclass.Sturges(x)
    if (type=="length") {
        rx = range(x, na.rm=TRUE)    
        dx = diff(rx)
        if (dx==0) dx = abs(rx[1])
        brks = seq.int(from=rx[1]-dx/1000, to=rx[2]+x/1000, length.out=n+1)
    } else if (type=="count") {
        brks = quantile(x, (0:n)/n, na.rm=TRUE)
    } else stop("Invalid type argument")
    attr(brks, "type") = type
    brks
}


###################################################
### chunk number 6: doLegend_Def
###################################################
#line 360 "annHeatmap.Rnw"
doLegend = function(r, col, side)
{
    nc = length(col)
    zval = seq(r[1], r[2], length=nc)
    z  = matrix(zval, ncol=1)
    if (side %in% c(1,3)) {
        image(x=zval, y=1, z=z, xaxt="n", yaxt="n", col=col)
    } else {
        image(x=1, y=zval, z=t(z), xaxt="n", yaxt="n", col=col)
    }        
    axis(side, las=1)
}


###################################################
### chunk number 7: convAnnData_Def
###################################################
#line 382 "annHeatmap.Rnw"
convAnnData = function(x, nval.fac=3)
{
    if (is.null(x)) return(NULL)
    x = as.data.frame(x)
    if (!is.null(nval.fac) & nval.fac>0) doConv = TRUE
    vv = colnames(x)
    for (v in vv) {
        if (is.logical(x[,v])) x[,v] = factor(as.numerical(x[,v]))
        if (doConv & length(unique(x[is.finite(x[,v]),v])) <= nval.fac) x[,v] = factor(x[,v])    
    }
    x
}


###################################################
### chunk number 8: cut.dendrogram_Def
###################################################
#line 404 "annHeatmap.Rnw"
cutree.dendrogram = function(x, h)
{
    cutx = cut(x, h)
    cutl = lapply(cutx$lower, getLeaves)
    nclus = sapply(cutl, length)
    ret   = rep(1:length(nclus), nclus)
    ret[order.dendrogram(x)] = ret
    names(ret)[order.dendrogram(x)] = unlist(cutl)
    ret = as.integer(factor(ret, levels=unique(ret))) # recode for order of clus
    ret
}


###################################################
### chunk number 9: getLeaves_Def
###################################################
#line 419 "annHeatmap.Rnw"
getLeaves = function(x)
{
    if (is.leaf(x)) attr(x, "label")
    else c(getLeaves(x[[1]]), getLeaves(x[[2]]))
}


###################################################
### chunk number 10: print.annHeatmap_Def
###################################################
#line 434 "annHeatmap.Rnw"
print.annHeatmap = function(x, ...)
{
    cat("annotated Heatmap\n\n")
    cat("Rows: "); show(x$dendrogram$Row$dendro)
    cat("\t", if (is.null(x$annotation$Row$data)) 0 else ncol(x$annotation$Row$data), " annotation variable(s)\n")
    cat("Cols: "); show(x$dendrogram$Col$dendro)
    cat("\t", if (is.null(x$annotation$Col$data)) 0 else ncol(x$annotation$Col$data), " annotation variable(s)\n")
    invisible(x)
}


###################################################
### chunk number 12: cutplot_dendrogam_Def
###################################################
#line 472 "annHeatmap.Rnw"
cutplot.dendrogram = function(x, h, cluscol, leaflab= "none", horiz=FALSE, lwd=3, ...)
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


###################################################
### chunk number 13: annHeatmap2_Def
###################################################
#line 538 "annHeatmap.Rnw"
annHeatmap2 = function(x, dendrogram, annotation, cluster, labels, scale=c("row", "col", "none"), col=heat.colors(12), equalize=FALSE, trim, legend=FALSE)
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


###################################################
### chunk number 14: plot.annHeatmap_Def
###################################################
#line 660 "annHeatmap.Rnw"
plot.annHeatmap = function(x, ...)
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
        if (x$layout$legend.side %in% c(1,3)) {
            par(mar=c(2, mmar[2], 2, mmar[4]))
        } else {
            par(mar=c(mmar[1], 2, mmar[3], 2))            
        }        
        doLegend(range(x$data$x, na.rm=TRUE), x$layout$legend.side, col=x$data$col)
    }    
    
    invisible(x)    
    
}




###################################################
### chunk number 15: Generics_Def
###################################################
#line 731 "annHeatmap.Rnw"
regHeatmap = function(x, ...) UseMethod("regHeatmap")
annHeatmap = function(x, ...) UseMethod("annHeatmap")


###################################################
### chunk number 16: regHeatmap_Def
###################################################
#line 741 "annHeatmap.Rnw"
regHeatmap.default = function(x, dendrogram=list(clustfun=hclust, distfun=dist, status="yes"), labels=NULL, legend=TRUE, ...)
{
    ret = annHeatmap2(x, dendrogram=dendrogram, annotation=NULL, cluster=NULL,  labels=labels, legend=legend, ...)
    ret
}


###################################################
### chunk number 17: annHeatmap_Def
###################################################
#line 754 "annHeatmap.Rnw"
annHeatmap.default = function(x, annotation, dendrogram=list(clustfun=hclust, distfun=dist, Col=list(status="yes"), Row=list(status="hidden")), cluster=NULL, labels=NULL, legend=TRUE, ...)
{
    ret = annHeatmap2(x, dendrogram=dendrogram, annotation=list(Col=list(data=annotation, fun=picketPlot)), cluster=cluster,  labels=labels, legend=TRUE, ...)
    ret
}


###################################################
### chunk number 18: annHeatmapExpressionSet_Def
###################################################
#line 766 "annHeatmap.Rnw"
annHeatmap.ExpressionSet =function(x, ...)
{
    expmat = exprs(x)
    anndat = pData(x)
    annHeatmap(expmat, anndat, ...)
}




