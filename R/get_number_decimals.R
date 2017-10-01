#' Get number of decimals for any vector of numbers (i.e. return total number of digits after decimal point)
#' 
#' @param x Numerical vector
#' 
#' @return Vector giving the number of decimals for each element in x
#'
#' @examples
#' get_number_decimals(c(1,2))
#' get_number_decimals(c(0.0006, 0.0750, 0.0420, 0.0031, 0.0001, 0.0100))
#' get_number_decimals(c(6.5*10^-5, 0.0100))
#' get_number_decimals(c(10.5,6.57*10^-5))
#' get_number_decimals(c(-1000.05,-6.57*10^-5))
#'
#' @export
get_number_decimals <- function(x)
{
  n_dec <- rep(NA, length=length(x))
  
  if((class(x) != "numeric") & (class(x) != "integer"))
  {
    stop("Input must be numeric")
  }
  if(class(x) == "integer")
  {
    n_dec[] <- 0
  } else {
    ##convert to absolute value
    x <- abs(x)
    ##get rid of whole part 
    x <- x-trunc(x)
    
    ##first convert into scientific notation
    sci_note_x <- as.character(format(x,scientific = TRUE))
    ##split on "e-"
    split_on_e <- strsplit(sci_note_x, "e-")
    ##get the negative exponents
    neg_exp <- sapply(split_on_e, function(x){as.numeric(x[2])})
    ##if NA, change to 0
    neg_exp[is.na(neg_exp)] <- 0
    
    ##get the number of digits after "."
    after_dot <- sapply(split_on_e, function(x){y <- strsplit(x[1], "\\.")[[1]][2]})
    ##remove the "0" at the end
    after_dot <- gsub("(0)+$","",after_dot)
    after_dot <- nchar(after_dot)
    after_dot[is.na(after_dot)] <- 0
    
    ##if NA, change to 0
    n_dec <- after_dot + neg_exp
  }
  n_dec
}

