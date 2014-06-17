\name{Setting defaults}
\alias{modifyExistingList}
\alias{extractArg}
\title{Override existing list entries}
\description{
Override existing list entries and extract arguments that are specified as named lists
}
\usage{
modifyExistingList(x, val)

extractArg(arglist, deflist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{a named list, the target for replacing with entries with the same name from \code{val}
}
    \item{val}{a named list that serves as template for filling in values in \code{x}
}
    \item{arglist}{a named list; these are the specified arguments that override the defaults.
}
  \item{deflist}{a named list whose entries are all possible slots (with default values) that can be filled. 
}
}
\details{
    \code{modifyExistingList} is a general function that recursively overwrites named items in \code{x} with the value of items of \code{val} with the same name. Items in \code{val} that have no name, or do not correspond to an item in \code{x} with the same name, are ignored.
    
    \code{extractArg} is a specific helper function for setting default values for the \code{annHeatmap2}-family of functions, where arguments are given as a list with two named items, \code{Row} and \code{Col}. Each of these items is again a named list of actual parameters. At the same time, all items with other names than \code{Row} and \code{Col} at the top level are assumed to be shared items with the same value for both sub-lists. \code{extractArg} uses \code{modifyExistingList} to overwrite the default values specified in \code{deflist} with the actual values specified in \code{arglist}, see Examples.
}
\value{
    \code{modifyExistingList} returns \code{x}, with values replaced from \code{val} where names match. \code{extractArg} returns a list with items \code{Row} and \code{Col} fully specified according to both \code{deflist} and \code{arglist}.
}
\seealso{
    \code{\link{annHeatmap2}}
}
\examples{
    ## Replace items with matching names recursively
    x   = list(a=1, b=2, c=list(a=31, b=32), 135)
    val = list(a=2, c=list(b=1114), d=92)
    modifyExistingList(x, val)
    
    ## Same defaults for rows/columns, no arguments specified
    defs = list(a="A", b="B", c="C")
    extractArg(NULL, defs)

    ## Shared and non-shared defaults
    defs = list(common.1=134, common.2=72, Row=list(row.only=14), Col=list(col.only=134))
    args = list(common.1 = -1, Row=list(row.only=94, common.2=-15))
    extractArg(args, defs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
