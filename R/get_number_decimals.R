#' Get number of decimals (i.e. return total number of digits after decimal point) for any vector of numbers in [0,1) if number of decimals <= 6
#' 
#' @param x Numerical vector where all elements are in [0,1)
#' 
#' @return Vector giving the number of decimals for each element in x if the number is <= 6; otherwise return 7 with a warning
#'
#' @examples
#' get_number_decimals(c(0.0006, 0.0750, 0.0420, 0.0031, 0.0001, 0.0100))
#' get_number_decimals(c(6.5*10^-4, 0.0100)) ##This does not work correctly!
#' get_number_decimals(c(6.5e-4, 0.0100))
#' get_number_decimals(c(0.00065, 0.0100))
#'
#' @export
get_number_decimals <- function(x)
{
  if((any(x<0))|(any(x>=1)))
  {
    stop("All elements of x should be in [0,1)")
  }
  
  ##get all numbers spaced 10^-k apart from 0 to 1
  list_grid <- lapply(1:6, function(k){(1:10^k)/(10^k)})
  
  ##function for a single number
  n_dec_single <- function(x_single){
    ##get which vector the query number is in, which corresponds to the "number of digits"
    grid_x_is_in <- sapply(list_grid, function(l,a){a %in% l}, x_single)
    
    if(x_single==0)
    {
      n_dec <- 0
    } else {
      if(sum(grid_x_is_in) >= 1)
      {
        n_dec <- min((1:6)[grid_x_is_in])
      } else {
        n_dec <- 7
        warning("Number of decimals seems to be > 6. This case is not implemented. Beware floating point arithmetic!")
      }
    }
    n_dec
  }
  n_d <- sapply(x, n_dec_single) 
  n_d
}
