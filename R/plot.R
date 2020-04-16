# Plot functions for lm_pi0 and lm_qvalue
#


#' default plotting style
#'
#' @keywords internal
#' @noRd
#' @importFrom Rcssplot Rcss 
cssSwfdr = Rcss(system.file("css", "swfdr.css", package="swfdr"))


#' Display a histogram of pi0 estimates
#'
#' @param x object of class lm_pi0
#' @param breaks integer or vector, position of histogram breakpoints
#' @param xlab character, label for x-axis
#' @param ylab character, label for y-axis
#' @param main character, label for title
#' @param css style object of class Rcss
#' @param cssclass character, style class
#' @param ... other arguments, passed to hist
#'
#' @importFrom Rcssplot par plot lines axis mtext hist
#' @importFrom Rcssplot RcssGetCompulsoryClass RcssValue
#' 
#' @examples
#' # generate a lm_pi0 object
#' pValues <- rbeta(1000, 0.5, 1)
#' X <- pValues + rnorm(1000, 0, 0.5)
#' pi0 <- lm_pi0(pValues, X=X)
#' # draw histogram of pi0 estimates
#' hist(pi0)
#' 
#' @method hist lm_pi0
#' @export
hist.lm_pi0 <- function(x, breaks=31,
                        xlab="pi0 estimate", ylab="Counts", main="",
                        css="swfdr", cssclass=NULL, ...) {
  
  if (identical(css, "swfdr")) {
    RcssDefaultStyle <- cssSwfdr
  } else {
    RcssDefaultStyle <- css
  }
  RcssCompulsoryClass <- RcssGetCompulsoryClass(c("swfdr", "hist", cssclass))
  
  pi0 <- regularize.interval(x$pi0)
  pi0.raw <- regularize.interval(x$pi0.raw)
  if (length(breaks)==1) {
    breaks <- seq(0, 1, length=breaks)
  }
  hcol <- RcssValue("hist", "col", default="#000000")
  
  # draw the distribution of pi0 estimates
  par()
  result <- hist(pi0, breaks=breaks, xlim=c(0, 1),
                 xlab="", ylab="", main="", ...)
  axis(1, at=c(0, 1), label=rep("", 2), line=0, Rcssclass="x")
  axis(1, label=NA, line=0, Rcssclass="x")
  axis(1, lwd=0, Rcssclass="x")
  axis(2, label=NA, line=0, Rcssclass="y")
  axis(2, lwd=0, Rcssclass="y")
  mtext(xlab, side=1, Rcssclass="xlab")
  mtext(ylab, side=2, Rcssclass="ylab")
  mtext(main, side=3, Rcssclass="main")
  
  # show the average pi0, i.e. the no-covariates pi0
  lines(rep(pi0.raw, 2), c(0, 1.05*max(result$counts)), Rcssclass="center")
  invisible(result)
}


#' Display a scatter-plot compring p-values and q-values
#'
#' @keywords internal
#' @param x object of class lm_pi0
#' @param var.x characer, variable to show on the x-axis
#' @param log character, use 'xy' to display both axes on logarithmic scale
#' @param threshold numeric, threshold to call significance is displayed using
#' lines, set to NULL to avoid showing the threshold lines
#' @param xlab character, label for x-axis
#' @param ylab character, label for y-axis
#' @param main character, label for plot title
#' @param xlim numeric of length 2, limits for x-axis
#' @param ylim numeric of length 2, limits for y-axis
#' @param css object of class Rcss
#' @param cssclass character, style class
#' @param ... other arguments, passed to points()
#'
#' @importFrom Rcssplot par plot lines axis mtext points abline
#' @importFrom Rcssplot RcssGetCompulsoryClass RcssValue
#'
#' @examples
#' # generate a lm_pi0 object
#' pValues <- rbeta(1000, 0.5, 1)
#' X <- pValues + rnorm(1000, 0, 0.5)
#' qValues <- lm_qvalue(pValues, X=X)
#' # draw a comparison of p-values and q-values
#' plot(qValues)
#' plot(qValues, threshold=NA)
#' 
#' @method plot lm_qvalue
#' @export
plot.lm_qvalue <- function(x, var.x=c("pvalues", "qvalues.raw"),
                           log="xy", threshold=0.05,
                           xlab=NULL, ylab=NULL, main="",
                           xlim=NULL, ylim=NULL, 
                           css="swfdr", cssclass=NULL, ...) {
  
  # adjust settings to show qvalue vs pvalue, or qvalue with vs without covariates
  var.x = match.arg(var.x)
  if (is.null(ylab)) {
    ylab <- "conditioned q-value"
  }
  if (var.x=="pvalues") {
    .data <- cbind(x$pvalues, x$qvalues)
    if (is.null(xlab)) {
      xlab <- "p-value"
    }
  } else if (var.x=="qvalues.raw") {
    .data <- cbind(x$qvalues.raw, x$qvalues)
    if (is.null(xlab)) {
      xlab <- "q-value"
    }
  }
  .min.nonzero <- min(.data[.data>0])
  if (is.null(xlim)) {
    xlim = c(.min.nonzero, 1)
  }
  if (is.null(ylim)) {
    ylim = xlim
  }

  if (identical(css, "swfdr")) {
    RcssDefaultStyle <- cssSwfdr
  } else {
    RcssDefaultStyle <- css
  }
  RcssCompulsoryClass <- RcssGetCompulsoryClass(c("swfdr", "scatter", cssclass))
  
  # creat the plot area and add points
  par()
  plot(.data[,1], .data[,2], log=log, type="n",
       xlim=xlim, ylim=ylim, 
       xlab="", ylab="", main="")
  axis(1, at=xlim, label=c("", ""), tck=0, line=0, Rcssclass="x") 
  axis(1, label=NA, line=0, Rcssclass="x")
  axis(1, lwd=0, Rcssclass="x")
  axis(2, at=ylim, label=c("", ""), tck=0, line=0, Rcssclass="y")
  axis(2, label=NA, line=0, Rcssclass="y")
  axis(2, lwd=0, Rcssclass="y")
  mtext(xlab, side=1, Rcssclass="xlab")
  mtext(ylab, side=2, Rcssclass="ylab")
  mtext(main, side=3, Rcssclass="main")
  
  count.class <- NULL
  if (length(x$pvalues)>1000) {
    count.class <- "many"
  }
  points(.data[,1], .data[,2], Rcssclass=count.class, ...)
  
  if (is.finite(threshold)) {
    abline(h=threshold, v=threshold, Rcssclass="threshold")
  }
}

