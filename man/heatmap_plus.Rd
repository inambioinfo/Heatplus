\name{heatmap_plus}
\alias{heatmap_plus}
\title{Display an Annotated Heatmap}
\description{
This function displays an expression data matrix as a heatmap with a column dendrogram. A given clustering will be shown in color.  Additionally, a number of binary and interval scaled covariates can be added to characterize these clusters.
}
\usage{
heatmap_plus(x, addvar, covariate = NULL, picket.control = list(), h, clus,
             cluscol, cluslabel = NULL, Rowv, Colv, reorder = c(TRUE, TRUE),
			 distfun = dist, hclustfun = hclust, scale = c("row", "column",
			 "none"), na.rm = TRUE, do.dendro = TRUE, col = heat.colors(12),
			 trim, equalize = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the numerical data matrix to be displayed.}
  \item{addvar}{data frame with (mostly binary) covariates.}
  \item{covariate}{integer indicating the one column in \code{addvar} that is interval scaled.}
  \item{picket.control}{list of option for drawing the covariates, passed to \code{picketplot}.}
  \item{h}{height at which to cut the dendrogram, as in \code{cutree}; overrides \code{clus}.}
  \item{clus}{an explicit vector of cluster memberships for the columns of \code{x}, if no dendrogram is used; ignored if \code{do.dendro=TRUE} and \code{h} is specified.}
  \item{cluscol}{a vector of colors used to indicate clusters.}
  \item{cluslabel}{labels to designate cluster names.}
  \item{Rowv}{either a dendrogram or a vector of reordering indexes for the rows.}
  \item{Colv}{either a dendrogram or a vector of reordering indexes for the columns.}
  \item{reorder}{logical vector of length two, indicating whether the rows and columns (in this order) should be reordered using \code{order.dendrogram}.}
  \item{distfun}{function to compute the distances between rows and columns. Defaults to \code{dist}.}
  \item{hclustfun}{function used to cluster rows and columns. Defaults to \code{hclust}.}
  \item{scale}{indicates whether values should be scaled by either by row, column, or not at all. Defaults to \code{row}.}
  \item{na.rm}{logical indicating whther to remove NAs.}
  \item{do.dendro}{logical indicating whether to draw the column dendrogram.}
  \item{col}{the color scheme for \code{image}. The default sucks.}
  \item{trim}{Percentage of values to be trimmed. This helps to keep an informative color scale, see Details.}
  \item{equalize}{logical indicating whther to use the ranks of the data for setting the color scheme; alternative to \code{trim}, see Details.}
  \item{\dots}{extra arguments to \code{image}.}
}
\details{
This is a heavily modified version of \code{heatmap_2}, which is a heavily modfied version of an old version of \code{heatmap} in package \code{stats}, so some of the arguments are described in more detail there. The main distinguishing feature of this routine is the possibility to color a cluster solution, and to add a covariate display.

Covariates are assumed to be binary, coded as 0 and 1 (or \code{FALSE} and \code{TRUE} respectively). One of the covariates can be interval scaled, the column index of this variable is supplied via argument \code{covariate}. The details of the added display are handled by the function \code{picketplot}.
  
Setting \code{trim} to a number between 0 and 1 uses equidistant classes between the (\code{trim})- and (1-\code{trim})-quantile, and lumps the values below and above this range into separate open-ended classes. If the data comes from a heavy-tailed distribution, this can save the display from putting too many values into to few classes. Alternatively, you can set \code{equal=TRUE}, which uses an equidistant color scheme for the ranks of the values.
}
\value{
A list with components
  \item{rowInd}{indices of the rows of the display in terms of the rows of \code{x}.}
  \item{colInd}{ditto for the columns of the display.}
  \item{clus}{the cluster indices of the columns of the display.}
}
\author{Original by Andy Liaw, with revisions by Robert Gentleman and Martin Maechler.

Alexander Ploner for the modifications documented here.
}
\seealso{\code{\link{heatmap_2}}, \code{\link{heatmap}}, \code{\link{picketplot}}, \code{\link{cutplot.dendrogram}},
         \code{\link{RGBColVec}}}
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
addvar = data.frame(Var1=rep(c(0,1,0),c(10,20,15)),
                    Var2=rep(c(1,0,0),c(10,20,15)),
                    Var3=rep(c(1,0), c(15,30)),
                    Var4=rep(seq(0,1,length=4), c(10,5,15,15))+rnorm(45, sd=0.5))
addvar[3,3] = addvar[17,2] = addvar[34,1] =NA
colnames(addvar) = c("Variable X","Variable Y", "ZZ","Interval")


# the lame default, without clustering
# Labels do not look too hot that way
heatmap_plus(mm)

# without labels, but with cluster
dimnames(mm)=NULL
heatmap_plus(mm, h=40)

# add some covariates, with nice names
heatmap_plus(mm, addvar=addvar, cov=4)

# covariates and clustering
heatmap_plus(mm, addvar=addvar, cov=4, h=20, col=RGBColVec(64), equal=TRUE)

# Clustering without the dendrogram
cc = cutree(hclust(dist(t(mm))), k=5)
heatmap_plus(mm, addvar=addvar, cov=4, clus=cc, do.dendro=FALSE)
}
\keyword{hplot}% at least one, from doc/KEYWORDS

