#' @title create a plot of estimated coefficients
#' @description
#' Plot point estimates of regression coefficients and 95\% credible intervals
#' 
#' @name plot.psbcSpeedUp
#' 
#' @importFrom stats quantile
#' @importFrom GGally ggcoef 
#' @importFrom ggplot2 xlab ylab
#' 
#' @param x an object of class \code{psbcSpeedUp}
#' @param type type of point estimates of regression coefficients. One of 
#' \code{c("mean", "median", "mode")}. Default is \code{mean}
#' @param interval logical argument to show 95\% credible intervals. Default 
#' is \code{TRUE}
#' @param ... other arguments
#' 
#'
#' @examples
#' 
#' # Load the example dataset
#' data("exampleData", package = "psbcSpeedUp")
#' p <- exampleData$p
#' q <- exampleData$q
#' survObj <- exampleData[1:3]
#' 
#' # Set hyperparameters
#' mypriorPara <- list('kappa0'=1, 'c0'=2, 'r'=10/9, 'delta'=1e-05, 'groupInd'=c(1:p),
#'                     'beta.prop.var'=1, 'beta.clin.var'=1, 'beta.ini'= rep(0,p+q), 
#'                     'lambdaSq'=1, 'sigmaSq'= runif(1, 0.1, 10))
#' 
#' # run Bayesian Lasso Cox
#' library(psbcSpeedUp)
#' set.seed(123)
#' fitBayesCox =  psbcSpeedUp(survObj, p=p, q=q, hyperpar=mypriorPara, 
#' nIter=10, burnin=0, outFilePath=tempdir())
#' plot(fitBayesCox)
#'
#' @export
plot.psbcSpeedUp <- function(x, type = "mean", interval = TRUE, ...) {
  if (!inherits(x, "psbcSpeedUp")) {
    stop("Use only with \"psbcSpeedUp\" objects")
  }

  if (length(type)==1) {
    if (! type %in% c("mean", "median", "mode")) {
      stop("'type' should be one of c('mean', 'median', 'mode')!")
    }
  } else {
    stop("'type' should be one of c('mean', 'median', 'mode')!")
  }
  
  if (!is.logical(interval))
    stop("Argument 'interval' must be a logical value!")
  
  #pdf("psbcBeta.pdf", height = 5, width = 3.5)
  if (is.null(colnames(x$output$beta.p))) {
    x_names <- paste0("x", 1:ncol(x$output$beta.p))
  } else {
    x_names <- colnames(x$output$beta.p)
  }

  beta_p <- x$output$beta.p[-(1:(x$input$burnin + 1)), ]
  beta_est <- apply(beta_p, 2, type)
  beta_L <- apply(beta_p, 2, quantile, 0.025)
  beta_U <- apply(beta_p, 2, quantile, 0.975)
  tbl <- data.frame(term = x_names, estimate = beta_est,  conf.low = beta_L,  conf.high = beta_U)
  tbl$term <- factor(tbl$term, levels = tbl$term)
  
  ggcoef(tbl) + xlab(expression(Posterior~~beta)) + ylab("")
  #dev.off()
}
