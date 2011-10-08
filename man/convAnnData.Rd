\name{convAnnData}
\alias{convAnnData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converting data frames for display as annotation}
\description{
Converts a data frames for display as annotation in a heatmap. This is mostly intended as an internal function, but might be useful for finetuning an annotation data frame manually. 
}
\usage{
convAnnData(x, nval.fac = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the data frame to be converted}
  \item{nval.fac}{lower limit for unique values in numerical variables}
}
\details{Logical variables are converted to factors. So are numerical variables with less than \code{nval.fac} unique values.
}
\value{\code{convAnnData} returns the converted data frame.}
\seealso{\code{\link{annHeatmap2}}}
\examples{
data(mtcars)
summary(mtcars)
summary(convAnnData(mtcars))
summary(convAnnData(mtcars, nval.fac=2))    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

