\name{picketplot}
\alias{picketplot}
\title{Barplots for Several Binary Variables}
\description{
Display one or more binary variables by using black bars for presence/validity
of a condition, empty space for absence/invalidity, and an extra color for missing values. Additionally, an index plot for one interval scaled variable can be added, possibly with a smoothing function.

This routine is primarily intended for augmenting heatmaps. It might be useful in other contexts, but misses most frills for using it comfortably.
}
\usage{
picketplot(x, covariate = NULL, grp = NULL, grpcol, grplabel = NULL, 
           add = FALSE, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix or data frame containing the data.}
  \item{covariate}{the index of the column in \code{x} that contains the interval scaled variable, if any.}
  \item{grp}{cluster indices for the rows of \code{x}, used for assigning background color.}
  \item{grpcol}{colors corresponding to the clusters.}
  \item{grplabel}{cluster names.}
  \item{add}{logical indicating whether to start a new plot, or whether to add the plot to the existing one.}
  \item{control}{a list of parameters controlling the appearance of the plot, see Details.}
}
\details{
The following named list elements can be set to change the appearance of the
plot: 
  \describe{
    \item{boxw}{the relative width of a marking box.}
    \item{boxh}{the relative height of a marking box.}
    \item{hbuff}{the horizontal separation around marking boxes; equals half the 
                 horizontal distance between two marking boxes.}
    \item{vbuff}{ditto for vertical separation.}
    \item{span}{passed on to \code{loess} used for the smoothing curve.}
    \item{nacol}{color for missing values of binary variables.}
    \item{degree}{if 0, no smoothing line is drawn; otherwise passed on to
    \code{loess} used for the smoothing curve.}
    \item{cex.label}{the character size for \code{grplabel}.}
  }
}
\author{Alexander Ploner <Alexander.Ploner@ki.se>}
\note{The plot looks like a more or less derelict picket fence, and 'picketplot' sounds somewhat like the 'pocketplot' used in geostatistics.}
\seealso{\code{\link{heatmap_plus}}}
\examples{
# without covariate
mm = cbind(sample(0:1, 42, rep=TRUE), sample(0:1, 42, rep=TRUE))
mm[sample(42, 5), 1] = NA
picketplot(mm)

# with clustering
cl = rep(1:3, c(10,22,10))
cn = c("Cluster I","Cluster II","Cluster III")
cc = c("lightblue","lightgreen","lightpink") # windows palette
picketplot(mm, grp=cl, grplabel=cn, grpcol=cc)

# add a covariate; setting the colnames makes the variable labels
mm = cbind(mm, rnorm(42) + cl/2)
colnames(mm) = c("State A","State B", "X")
picketplot(mm, covariate=3,grp=cl, grplabel=cn, grpcol=cc)

# using extra controls
picketplot(mm, covariate=3,grp=cl, grplabel=cn, grpcol=cc, control=list(nacol="white", degree=0))

}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{aplot}% __ONLY ONE__ keyword per line
