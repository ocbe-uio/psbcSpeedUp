#' psbcSpeedUp
#'
#' @title Time-dependent Brier scores
#'
#' @description
#' Predict time-dependent Brier scores based on Cox regression models
#'
#' @name plotBrier
#' 
#' @importFrom ggplot2 ggplot aes geom_step theme element_blank
#' 
#' @param object fitted object obtained with \code{psbcSpeedUp}
#' @param survObj.new a list containing observed data from new subjects with
#' components \code{t}, \code{di}, \code{x}. If \code{NULL}, the prediction is
#' based on the training data
#' @param times maximum time point to evaluate the prediction
#' @param method option to use the posterior mean (\code{"mean"}) of coefficients
#' for prediction or Bayesian model averaging (\code{"BMA"}) for prediction
#' @param \dots not used
#'
#' @keywords survival
##' @examples
##'
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
#' \donttest{
#' # run Bayesian Lasso Cox
#' library("psbcSpeedUp")
#' set.seed(123)
#' fitBayesCox <- psbcSpeedUp(survObj,
#'   p = p, q = q, hyperpar = mypriorPara,
#'   nIter = 10, burnin = 0, outFilePath = tempdir()
#' )
#' # predict survival probabilities of the train data
#' plotBrier(fitBayesCox, times = 80)
#' }
#'
#' @export
plotBrier <- function(object, survObj.new = NULL, 
                      method = "mean", times = NULL, ...) {
  
  Brier_score <- predict.psbcSpeedUp(object, 
                                      survObj.new = survObj.new, 
                                      method = method, 
                                      times = times)
  Brier <- model <- NULL
  #Brier_score %>%
    ggplot2::ggplot(Brier_score, 
                    aes(times, Brier, group = model, color = model)) +
    xlab("Evaluation time points") +
    ylab("Brier score") +
    geom_step(direction = "vh") +
    theme(legend.title = element_blank())
}
