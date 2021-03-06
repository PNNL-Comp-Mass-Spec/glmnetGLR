\name{parseJob}
\alias{parseJob}
\title{Parses a collection  of elements into (almost) equal sized groups}
\usage{
parseJob(n, njobs, collate = FALSE, random.seed = NULL,
  text.to.eval = FALSE)
}
\arguments{
  \item{n}{The number elements to be parsed}

  \item{njobs}{The number of groups}

  \item{collate}{\code{=TRUE} alternative ordering of the
  grouping. See example below.}

  \item{random.seed}{An integer setting the random seed,
  which will result in randomizing the elements among the
  jobs.  If \code{NULL}, no randomization is performed.
  Randomization cannot be performed if \code{collate =
  TRUE} or if \code{text.to.eval = TRUE}. Randomization is
  useful when the computing time for each element varies
  significantly because it helps to even out the run times
  of parallel jobs.}

  \item{text.to.eval}{If \code{= TRUE}, a text expression
  is returned, that when evaluated, will produce the
  sequence of elements for that group.  This is especially
  useful when \code{n} is very large. (See \code{Value}
  section below).}
}
\value{
When \code{text.to.eval = FALSE}, a list with \code{njobs}
elements is returned, each element containing a numeric
vector of the element numbers which correspond to that
group.  When \code{text.to.eval = TRUE}, a list with
\code{njobs} elements is returned, each element containing
text (let's call it \code{val}), that when evaluated using
\code{eval(parse(text = val))}, will produce the sequence
of numbers corresponding to the group.
}
\description{
Parses a collection of elements into (almost) equal sized
groups. Useful for splitting up an R job that operates over
a large dataframe or list into multiple jobs.
}
\examples{
x <- parseJob(29, 6)
print(x)

# To see the range of each group
lapply(x, range)

# To see the length of each group
lapply(x, length)

# Randomize the outcome
parseJob(32, 5, random.seed = 231)

# Example of text.to.eval
parseJob(11, 3, text.to.eval = TRUE)
lapply(parseJob(11, 3, text.to.eval = TRUE), function(x) eval(parse(text = x)))

# Example of collated and text.to.eval
parseJob(11, 3, collate = TRUE)
parseJob(11, 3, collate = TRUE, text.to.eval = TRUE)
}
\author{
Landon Sego
}

