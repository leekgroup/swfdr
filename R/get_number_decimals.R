#' Get number of decimals (i.e. return total number of digits after decimal point) for any vector of numbers in [0,1) if number of decimals <= 6
#' 
#' @param x Numerical vector where all elements are in [0,1)
#' 
#' @return Vector giving the number of decimals for each element in x if the number is <= 6; otherwise return 7 with a warning
#'
#' @examples
#' get_number_decimals(c(0.0006, 0.0750, 0.0420, 0.0031, 0.0001, 0.0100))
#' get_number_decimals(c(6*10^-4, 7.5*10^-2, 4.2*10^-2, 3.1*10^-3, 10^-4, 10^-2))
#' get_number_decimals(c(6.5*10^-4, 0.0100)) 
#' get_number_decimals(c(6.5e-4, 0.0100))
#' get_number_decimals(c(0.00065, 0.0100))
#' get_number_decimals(c(10^-7, 10e-7, 10e-3))
#'
#' @export
get_number_decimals <- function(x)
{
  if((any(x<0))|(any(x>=1)))
  {
    stop("All elements of x should be in [0,1)")
  }
  
  ##get the maximum number of digits
  max_digits <- 6
  
  ##get all numbers spaced 10^-k apart from 0 to 1
  list_grid <- lapply(1:max_digits, function(k){(1:10^k)/(10^k)})
  
  ##function for a single number
  n_dec_single <- function(x_single){
    ##round to get rid of funny numerical issues
    x_single <- round(x_single, 12)
    
    if(x_single < 10^-max_digits)
    {
      n_dec <- max_digits + 1
      warning(paste(max_digits + 1, " is a place-holder. Number of decimals seems to be >", 
                    max_digits, ". This case is not implemented. Beware floating point arithmetic!",
                    sep=""))
    } else {
      ##get which vector the query number is in, which corresponds to the "number of digits"
      grid_x_is_in <- sapply(list_grid, function(l,a){min(abs(l-a)) < 10^-12}, x_single)
      
      if(x_single==0)
      {
        n_dec <- 0
      } else {
        if(sum(grid_x_is_in) >= 1)
        {
          n_dec <- min((1:max_digits)[grid_x_is_in])
        } else {
          n_dec <- max_digits + 1
          warning(paste(max_digits + 1, " is a place-holder. Number of decimals seems to be >", 
                        max_digits, ". This case is not implemented. Beware floating point arithmetic!",
                        sep=""))        
        }
      }
    }
    n_dec
  }
  n_d <- sapply(x, n_dec_single) 
  n_d
}
