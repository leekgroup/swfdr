# Functions for validation of inputs for lm_pi0 and lm_qvalue
#
# All these functions are for internal use for the package only


#' validate pvalues. They must be finite, in range [0,1]
#'
#' @keywords internal
#' @noRd
#' @param p numeric vector of p-values
#'
#' @return numeric vector of length 2 with (min(p), max(p))
check_p <- function(p) {
  if (missing(p)) {
    stop("p is a required argument\n", call. = FALSE)
  }
  if (!is(p, "numeric")) {
    stop("p must be a numeric vector\n", call. = FALSE)
  }
  prange <- range(p, na.rm=TRUE)
  if (prange[1] < 0 | prange[2] > 1) {
    stop("p must be numeric in range [0, 1]\n", call. = FALSE)
  }
  if (any(is.na(p))) {
    stop("p must not have missing values\n", call. = FALSE)
  }
  prange
}


#' validate lambda. They must be numeric, finite, sorted with unique value
#' 
#' @keywords internal
#' @noRd
#' @param x vector of lambda values
#' @param pmax numeric, maximal pvalue
#'
#' @return numeric vector of sorted unique x,
#' or stop if x does not satisfy criteria
check_lambda <- function(x, pmax) {
  if (!is(x, "numeric")) {
    stop("lambda must be a numeric vector \n", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("lambda must not contain NAs, NULL, or non-finite elements\n",
         call. = FALSE)
  }
  x <- sort(unique(x))
  if (length(x)<4) {
    stop("lambda must be of length >=4\n", call. = FALSE)
  }
  if (min(x)<0 | max(x)>=1) {
    stop("lambda values must all be in range [0, 1)\n", call. = FALSE)
  }
  if (pmax < max(x)) {
    warning("maximal p is smaller than maximal lambda", call. = FALSE)
  }
  x
}


#' validate degrees of freedom. 
#'
#' @keywords internal
#' @noRd
#' @param x expect a single number
#' @param max.value numeric, maximal value allowed for x
#'
#' @return integer derived from x
check_df <- function(x, max.value) {
  if (!is(x, "numeric") & !is(x, "integer")) {
    stop("df must be a number")
  }
  if (length(x) != 1 | any(!is.finite(x))) {
    stop("df must be a single finite number", call. = FALSE)
  }
  x <- round(x)
  if (x <= 1 | x > max.value) {
    stop("df must be in range 1 < df < length(lambda)\n", call. = FALSE)
  }
  x
}


#' validate matrix of covariates.
#'
#' The matrix must be compatible with a vector of pvalues
#'
#' @keywords internal
#' @noRd
#' @param X vector or matrix of covariates
#' @param p vector of p-values
#'
#' @return matrix
check_X <- function(X, p) {
  # allow for null input (no covariates)
  if (missing(X)) {
    X <- NULL
  }
  if (is.null(X)) {
    warning("X is missing or NULL - modeling will have no effect",
            call.=FALSE)
    X <- cbind(rep(1, length(p)))
    rownames(X) <- names(p)
  }
  # allow for a single covariates specified as a vector
  if (is.null(dim(X))) {
    X <- cbind(X=X)
  }
  # ensure that X and pvalues are compatible
  if (!is(X, "matrix")) {
    if (is(X, "data.frame")) {
      X <- as.matrix(X)
    }
  }
  if (length(p) != nrow(X)) {
    stop("incompatible X and p - different lengths\n", call. = FALSE)
  }
  if (!is.null(names(p)) & !identical(rownames(X), names(p))) {
    stop("X and p have different names", call. = FALSE)
  }
  # ensure that all columns in X are numeric
  is.number.class = function(z) {
    is(z, "numeric") | is(z, "integer") | is(z, "factor")
  }
  if (!all(apply(X, 2, is.number.class))) {
    stop("X must be a numeric vector or numeric matrix\n", call. = FALSE)
  }
  # ensure that all values are set
  num.bad <- sum(!is.finite(X))
  if (num.bad>0) {
    stop("X must not contain missing or non-finite values\n", call. = FALSE)
  }
  X
}


#' check if an object is of a certain class
#'
#' @keywords internal
#' @noRd
#' @param x object
#' @param classname character
#' 
#' @return nothing, emit error if check not satisfied
check_class <- function(x, classname) {
  if (!is(x, classname)) {
    stop(paste0("object is not of class '", classname, "'\n"), call. = FALSE)
  }
}

