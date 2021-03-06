##' Silent wrapper for coercing a vector to numeric
##'
##' @param x vector of any type
##' 
##' @return If \code{as.numeric(x)} produces an error or warning,
##' \code{x} is returned unchanged.  Otherwise, \code{as.numeric(x)} is returned.
##' 
##' @seealso \code{\link{as.numeric}}
##' @examples
##' as.numeric.silent(c("this","that"))
##' as.numeric.silent(c("2893.9","9423.48"))

# Function to convert a vector to numeric (without producing errors or warnings)
as.numeric.silent <- function(x) {

  # Set warnings to errors
  op <- options(warn = 2)

  # If an error results, the data aren't numeric
  if (class(x.num <- try(as.numeric(x), silent = TRUE)) == "try-error")
    x.num <- x

  # Restore the previous setting
  options(op)

  return(x.num)
    
} # as.numeric.silent()
