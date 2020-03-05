# Functions for display of package objects via print()


#' Display a summary of object lm_pi0
#'
#' @keywords internal
#' @noRd
#' @param x object of class lm_pi0
#' @param ... other arguments ignored
#'
#' @method print lm_pi0
#' @export
print.lm_pi0 <- function(x, ...) {
  check_class(x, "lm_pi0")
  components <- get.components(c("call", "lambda", "X", "pi0"),
                               list(...)[["components"]])
  compound.message(x, components)
  invisible(x)
}


#' Display a summary of an lm_qvalue object
#'
#' @keywords internal
#' @noRd
#' @param x lm_qvalue object
#' @param ... ignored
#'
#' @method print lm_qvalue
#' @export
print.lm_qvalue <- function(x, ...) {
  check_class(x, "lm_qvalue")
  components <- get.components(c("call", "lambda", "X", "pi0", "hits"),
                               list(...)[["components"]])
  compound.message(x, components)
  invisible(x)
}


#' Display a summary of object lm_pi0
#'
#' (This is the same as print.lm_pi0)
#'
#' @keywords internal
#' @noRd
#' @param object object of class lm_pi0
#' @param ... other arguments ignored
#'
#' @method summary lm_pi0
#' @export
summary.lm_pi0 <- function(object, ...) {
  invisible(print.lm_pi0(object))
}


#' Display a summary of object lm_qvalue
#'
#' (This is the same as print.lm_qvalue)
#'
#' @keywords internal
#' @noRd
#' @param object object of class lm_qvalue
#' @param ... other arguments ignored
#'
#' @method summary lm_qvalue
#' @export
summary.lm_qvalue <- function(object, ...) {
  invisible(print.lm_qvalue(object))
}




# #############################################################################
# Some helper functions to the print() and summary()


#' Helper to create a single string by concatenating items from a vector
#'
#' @keywords internal
#' @noRd
#' @param x vector of things
#' @param width vector of character widths for each item in x
#'
#' @return a single string
v2s <- function(x, width=8) {
  empty <- paste(rep(" ", max(width)), collapse="")
  xlen <- length(x)
  if (length(width)<xlen) {
    width <- rep(width, length=xlen)[1:xlen]
  }
  if (is(x, "numeric")) {
    x <- as.character(round(x, 4))
  }
  result <- as.character(x)
  for (i in seq_along(x)) {
    ichars <- nchar(x[i])
    if (ichars < width[i]) {
      result[i] <- paste0(substr(empty, 1, width[i]-ichars), result[i])
    }
  }
  paste(result, collapse=" ")
}


#' Get a set intersection, but when second set is null default to first set
#'
#' This is meant to identify a subset of supported features that are requested
#'
#' @keywords internal
#' @noRd
#' @param supported vector of supported feature names
#' @param requested vector of requested feature names
#'
#' @return character vector with an intersection, or all supported features
get.components <- function(supported, requested) {
  if (is.null(requested)) {
    return(supported)
  }
  intersect(supported, requested)
}


#' Compose a two line report about a numeric vector
#'
#' @keywords internal
#' @noRd
#' @param v numeric vector
#'
#' @return vector with two strings a header line and a data line
compose.stats <- function(v) {
  header <- c("(Length)", "Min", "Mean", "Median", "Max")
  data <- c(length(v), min(v), mean(v), median(v), max(v))
  c(v2s(header), v2s(data))
}


#' Compose and output a compound message and output
#'
#' @keywords internal
#' @noRd
#' @param x list object of type lm_qvalue or lm_pi0 (not checked)
#' @param components character vector, identifiers suggesting what to include
#' in output
#'
#' @return character vector
compound.message <- function(x, components=c("call", "lambda",
                                             "X", "pi0", "hits")) {
  comps = components
  output = setNames(vector("list", length=length(comps)), comps)
  if ("call" %in% comps) {
    output$call = c("Call:", x$call)
  }
  if ("lambda" %in% comps) {
    output$lambda =   c("lambda:", compose.stats(x$lambda))
  }
  if ("X" %in% comps) {
    output$X = c("covariates:", "(Length)", v2s(length(x$X.names)))
  }
  if ("pi0" %in% comps) {
    output$pi0 = c("pi0:", compose.stats(x$pi0))
  }
  if ("hits" %in% comps) {
    hits.header <- c(" ", "<1e-4", "<1e-3", "<0.01", "<0.05", "<0.1", "<1")
    hits.widths = c(9, rep(7, 6))
    thresholds = c(1e-4, 1e-3, 1e-2, 0.05, 0.1, 1)
    hits.p <- c("p-value",
                vapply(thresholds,
                       function(t) { sum(x$pvalues<t) },
                       integer(1)))
    hits.q <- c("q-value",
                vapply(thresholds,
                       function(t) { sum(x$qvalues<t) },
                       integer(1)
                       ))
    output$hits <- c("Cumulative number of significant calls:",
                     v2s(hits.header, hits.widths),
                     v2s(hits.p, hits.widths),
                     v2s(hits.q, hits.widths))
  }
  # add a separator line
  output = lapply(output, function(x) { c(x, "") })
  message(paste(c("", unlist(output)), collapse="\n")) 
}

