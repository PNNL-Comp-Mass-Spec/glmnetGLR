##' Converts all factor variables in a dataframe to character variables
##' 
##' @param dframe A dataframe
##' 
##' @return The same dataframe with the factor variables converted to character variables.
##' 
##' @author Landon Sego
##' 
##' @examples
##' x <- data.frame(a=factor(c(rep(1,4),rep(2,4),rep(3,4))), y=rnorm(12))
##' str(x)
##' x <- factor2character(x)
##' str(x)


# Function for converting all factors to characters in a data frame

factor2character <- function(dframe) {

  if (!is.data.frame(dframe))
    stop("'", deparse(substitute(dframe)), "' must be a dataframe.\n")

  for (cname in colnames(dframe)) {

    if (is.factor(dframe[,cname]))
      dframe[,cname] <- as.character(dframe[,cname])
    
  }

  return(dframe)

} # factor2character
