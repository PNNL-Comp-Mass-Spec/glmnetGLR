\name{sort.data.frame}
\alias{sort.data.frame}
\title{Sorts a dataframe by one or more variables}
\usage{
sort.data.frame(form, dat)
sort.data.frame(dat, form)
}
\arguments{
  \item{form}{A formula with the variable names to use for
  sorting}

  \item{dat}{The dataframe to sort}
}
\value{
The sorted dataframe
}
\description{
Sorts a dataframe by one or more variables, ascending or
descending The order of the arguments \code{form} and
\code{dat} is interchangeable.  In the formula, use + for
ascending, - for decending. Sorting is left to right in the
formula.  For example, to sort the dataframe \code{Oats}
with sorting variables \code{Block} and \code{Variety},
useage is either of the following:
}
\details{
\code{sort.data.frame(~ Block - Variety, Oats)} or

\code{sort.data.frame(Oats, ~ - Variety + Block)}
}
\examples{
d <- data.frame(b=factor(c("Hi","Med","Hi","Low"),levels=c("Low","Med","Hi"),
                ordered=TRUE),
                x=c("A","D","A","C"),y=c(8,3,9,9),z=c(1,1,1,2))

# Sort by descending z, descending b
sort.data.frame(~-z-b,d)

# Sort by ascending x, ascending y, and ascending z
sort.data.frame(~x+y+z,d)

# Sort by descending x, ascending y, ascending z
sort.data.frame(~-x+y+z,d)

# Sort by ascending x, descending y, ascending z
sort.data.frame(d,~x-y+z)
}
\author{
Kevin Wright with ideas from Any Liaw and small edits by
Landon Sego
}

