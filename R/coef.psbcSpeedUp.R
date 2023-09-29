#' @title coef method for class \code{psbcSpeedUp}
#'
#' @description
#' Extract the point estimates of the regression coefficients
#'
#' @name coef.psbcSpeedUp
#'
#' @param object an object of class \code{psbcSpeedUp}
#' @param type type of point estimates of regressions. One of
#' \code{c("mean", "median")}. Default is \code{mean}
#' @param ... not used
#'
#' @return Estimated coefficients are from an object of class \code{psbcSpeedUp}.
#' If the \code{psbcSpeedUp} specified data standardization, the fitted values
#' are base based on standardized data.
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
#' mypriorPara <- list(
#'   "groupInd" = 1:p, "eta0" = 0.02, "kappa0" = 1, "c0" = 2, "r" = 10 / 9, 
#'   "delta" = 1e-05, "lambdaSq" = 1, "sigmaSq" = runif(1, 0.1, 10), 
#'   "beta.prop.var" = 1, "beta.clin.var" = 1)
#'
#' # run Bayesian Lasso Cox
#' library(psbcSpeedUp)
#' set.seed(123)
#' fitBayesCox <- psbcSpeedUp(survObj,
#'   p = p, q = q, hyperpar = mypriorPara,
#'   nIter = 10, burnin = 0, outFilePath = tempdir()
#' )
#' coef(fitBayesCox)
#'
#' @export
coef.psbcSpeedUp <- function(object, type = "mean", ...) {
  if (!inherits(object, "psbcSpeedUp")) {
    stop("Use only with \"psbcSpeedUp\" objects")
  }

  if (length(type) == 1) {
    if (!type %in% c("mean", "median")) {
      stop("'type' should be one of c('mean', 'median')!")
    }
  } else {
    stop("'type' should be one of c('mean', 'median')!")
  }

  beta_p <- object$output$beta.p[-(1:(object$input$burnin + 1)), ]
  beta_est <- apply(beta_p, 2, type)

  return(beta_est)
}
