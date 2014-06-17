\name{convAnnData}
\alias{convAnnData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converting data frames for display as annotation}
\description{
Converts a data frames for display as annotation in a heatmap. This is mostly intended as an internal function, but might be useful for finetuning an annotation data frame manually. 
}
\usage{
convAnnData(x, nval.fac = 3, inclRef = TRUE, asIs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the data frame to be converted}
  \item{nval.fac}{lower limit for unique values in numerical variables}
  \item{inclRef}{logical value indicating whether to include the reference level among the dummy variables for factors}
  \item{asIs}{logical value indicating whether to perform a conversion; if \code{TRUE}, the input \code{x} is  simply returned, provided it is a numerical matrix (otherwise, the function stops with an error message)}
}
\details{Logical variables are converted to factors. So are numerical variables with less than \code{nval.fac} unique values.
}
\value{\code{convAnnData} returns the converted data frame, which is a numerical matrix}
\seealso{\code{\link{annHeatmap2}}}
\examples{
data(mtcars)
summary(mtcars)
summary(convAnnData(mtcars))
summary(convAnnData(mtcars, nval.fac=2))    
summary(convAnnData(mtcars, nval.fac=2, inclRef=FALSE))    

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

