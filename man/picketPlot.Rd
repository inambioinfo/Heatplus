\name{picketPlot}
\alias{picketPlot}
\title{Display a data frame of annotation information}
\description{
Displays a data frame of both factor and numerical variables in parallel panels. Factors levels are indicated by black rectangles, using dummy variables for more than two levels. Numerical variables are shown as simple index plots with an optional loess smoother. Panels can be arranged horizontally or vertically, and different groups of subjects can be indicated through different background colors. 
}
\usage{
picketPlot(x, grp = NULL, grpcol, grplabel = NULL, horizontal = TRUE, asIs = FALSE, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{usually a data frame, which is passed to \code{convAnnData} to be converted to a numerical matrix with dummy coding for the factor levels. Alternatively, such a numerical matrix can be constructed manually and passed in as \code{x}, see argument \code{asIs}}
    \item{grp}{an optional vector of cluster memberships, in the same order as the rows of \code{x}}
    \item{grpcol}{an optional vector of background colors for the clusters specified in \code{grp}}
    \item{grplabel}{an optional vector of names for the clusters specified in \code{grp}}
    \item{horizontal}{logical value whether to plot variables horizontally (default) or vertically}
    \item{asIs}{a logical value indicating whether \code{x} should be passed to \code{convAnnData} for pre-processing or not. Defaults to \code{FALSE}.}
  \item{control}{a list of control parameters that determines the appearance of the plot
      \itemize{
          \item{\code{boxw} is the relative length of the short side of a box marking (width for a horizontal plot). Default is 1.}
          \item{\code{boxh} is the relative length of the long side of a box marking (default: 4)}
          \item{\code{hbuff} is the relative distance between two box markings for the same variable (horizontal buffer for a horizontal plot). Default is 0.1}
          \item{\code{vbuff} is the relative distance between two box markings for the same subject, but different variables (default: 0.1)}
          \item{\code{cex.label} is the expansion factor for plotting cluster labels}
          \item{\code{numfac} is the expansion factor indicating how much higher (for a horizontal plot) or wider (for a vertical plot) panels with numerical variables are than panels for factor variables.}
          \item{\code{nacol} is the color for box markings indicating missing values (default: \code{gray(0.85)})}
          \item{\code{span} is the span argument for the loess smoother. Default is 1/3; setting this to zero switches off smoothing.}
          \item{\code{degree} is the degree of loess smoothing. Default is 1; setting this to zero switches off smoothing}
          \item{\code{pch} is the plotting character for numerical variables}
          \item{\code{cex.pch} is the size of the plotting character for numerical variables}
          \item{\code{col.pch} is the color of the plotting character for numerical variables}          
      }
  }
}
\details{
    Missing values are indicated by a box marking in \code{nacol} for factor values.
}
\value{Depending on data and arguments, the value of the last \code{axis} function call. Irrelevant, as this function is called for its side effect of producing a plot.}
\seealso{
\code{\link{annHeatmap2}}, \code{\link{convAnnData}}, \code{\link{par}} 
}
\examples{
    ## Standard call
    data(mtcars)
    picketPlot(mtcars)
    
    ## Pre-process the data for display
    mm = convAnnData(mtcars, inclRef=FALSE)
    picketPlot(mm, asIs=TRUE)
    
    ## Higher panels for continous traits
    picketPlot(mm, asIs=TRUE, control=list(numfac=3))
 
    ## With clusters
    picketPlot(mtcars, grp = rep(1:2, c(16, 16)), grpcol = c("pink","lightblue"), grplabel=c("Cluster 1", "Cluster 2"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

