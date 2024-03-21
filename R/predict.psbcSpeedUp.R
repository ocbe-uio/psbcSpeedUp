#' psbcSpeedUp
#'
#' @title Predict survival risk
#'
#' @description
#' Predict survival probability, (cumulative) hazard or (integrated) Brier scores based on Cox regression models
#'
#' @name predict.psbcSpeedUp
#' 
#' @importFrom riskRegression Score predictCox
#' @importFrom utils read.table
#' @importFrom survival coxph Surv
#'
#' @param object fitted object obtained with \code{psbcSpeedUp}
#' @param survObj.new a list containing observed data from new subjects with
#' components \code{t}, \code{di}, \code{x}. If \code{NULL}, the prediction is
#' based on the training data. If \code{type} is among
#' \code{c("hazard", "cumhazard", "survival")}, only \code{survObj.new$x} is needed
#' @param times time points at which to evaluate the risks. If \code{NULL}
#' (default), the event/censoring times are used. If \code{type="brier"}, the
#' largest one of the \code{times} is used
#' @param type option to chose for predicting survival probabilities (one of
#' \code{c('hazard','cumhazard','survival')}) or brier scores (\code{type="brier"})
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
#' library("survival")
#' set.seed(123)
#' fitBayesCox <- psbcSpeedUp(survObj,
#'   p = p, q = q, hyperpar = mypriorPara,
#'   nIter = 10, burnin = 0, outFilePath = tempdir()
#' )
#' # predict survival probabilities of the train data
#' predict(fitBayesCox)
#' }
#'
#' @export
predict.psbcSpeedUp <- function(object, survObj.new = NULL, type = "brier",
                                method = "mean", times = NULL, ...) {
  if (!inherits(object, "psbcSpeedUp")) {
    stop("Use only with 'psbcSpeedUp' object!")
  }

  if (any(!type %in% c("brier", "hazard", "cumhazard", "survival"))) {
    stop("Argument 'type' has to be among 'brier', 'hazard', 'cumhazard'' or 'survival'!")
  }

  if (!method %in% c("mean", "BMA")) {
    stop("Argument 'method' has to be 'mean' or 'BMA'!")
  }

  survObj_train <- as.matrix(read.table(paste0(object$input$outFilePath, "data.txt")))
  # Check the survival input object
  if (!is.null(survObj.new)) {
    if (!is.list(survObj.new)) {
      stop("Argument 'survObj' must be an input")
    }

    if (length(type) == 1 && type == "brier") {
      if (any(!names(survObj.new) %in% c("t", "di", "x"))) {
        stop("List 'survObj' must have three compoents 't', 'di' and 'x'!")
      }
    } else {
      if ("x" %in% names(survObj.new)) {
        stop("List 'survObj' must have a compoent 'x'!")
      }
    }

    if (is.numeric(survObj.new$x) && (is.data.frame(survObj.new$x) || is.matrix(survObj.new$x))) {
      if (ncol(survObj.new$x) != ncol(object$output$beta.p)) {
        stop("The 'survObj.new$x' has wrong columns!")
      }
      if (is.data.frame(survObj.new$x)) {
        survObj.new$x <- data.matrix(survObj.new$x)
      }
    } else {
      stop("Data 'survObj.new$x' must be a numeric dataframe or matrix!")
    }
  } else {
    survObj.new <- list()
    survObj.new$t <- survObj_train[, 1]
    survObj.new$di <- survObj_train[, 2]
    survObj.new$x <- survObj_train[, -c(1:2)]
  }

  betas <- object$output$beta.p[-(1:(object$input$burnin / object$input$thin + 1)), ]
  ibs <- data.frame("Null model" = NA, "Bayesian Cox" = NA)

  if (is.null(times)) {
    times <- sort(unique(survObj.new$t))
  }
  if (length(type) == 1 && type == "brier") {
    times <- seq(0, max(times), length = 100)
  }

  if (!"brier" %in% type) {
    if (method == "mean") {
      beta_m <- colMeans(betas)
      lp <- as.vector(survObj_train[, -c(1:2)] %*% beta_m)
      data_train <- data.frame(
        time = survObj_train[, 1],
        status = survObj_train[, 2], lp = lp
      )
      model_train <- coxph(Surv(time, status) ~ lp,
        data = data_train, y = TRUE, x = TRUE
      )

      lp_test <- as.vector(survObj.new$x %*% beta_m)
      # data_test <- data.frame(time = survObj.new$t, status = survObj.new$di,
      #                         lp = lp_test)
      fit <- predictCox(model_train,
        times = times, type = type,
        newdata = data.frame(lp = lp_test)
      )
      fit[names(fit) %in% c(
        "lastEventTime", "nTimes",
        "var.lp", "var.strata"
      )] <- NULL
    } else {
      # obtain linear predictors corresponding to MCMC estimates
      lp_all_train <- survObj_train[, -c(1:2)] %*% t(betas)
      lp_all_test <- survObj.new$x %*% t(betas)
      data_train <- data.frame(
        time = survObj_train[, 1],
        status = survObj_train[, 2],
        lp = lp_all_train[, 1]
      )
      # data_test <- data.frame(time = survObj.new$t,
      #                         status = survObj.new$di,
      #                         lp = lp_all_test[, 1])
      model_train <- coxph(Surv(time, status) ~ lp,
        data = data_train,
        y = TRUE, x = TRUE
      )
      # calculate survival prediction based on the 1st MCMC estimates
      fit <- predictCox(model_train,
        times = times, type = type,
        newdata = data.frame(lp = lp_all_test[, 1])
      )
      fit0 <- fit[names(fit) %in% type]

      # calculate Brier scores based on other MCMC estimates
      for (i in 2:nrow(betas)) {
        data_train$lp <- lp_all_train[, i]
        # data_test$lp <- lp_all_test[, i]
        model_train <- coxph(Surv(time, status) ~ lp,
          data = data_train,
          y = TRUE, x = TRUE
        )
        fit <- predictCox(model_train,
          times = times, type = type,
          newdata = data.frame(lp = lp_all_test[, i])
        )
        fit0 <- Map("+", fit0, fit[names(fit) %in% type])
      }
      fit0 <- Map("/", fit0, nrow(betas))
      fit[names(fit) %in% type] <- fit0

      fit[names(fit) %in% c(
        "lastEventTime", "nTimes",
        "var.lp", "var.strata"
      )] <- NULL
    }
    fit
  } else {
    if (method == "mean") {
      beta_m <- colMeans(betas)
      lp <- as.vector(survObj_train[, -c(1:2)] %*% beta_m)
      data_train <- data.frame(time = survObj_train[, 1], status = survObj_train[, 2], lp = lp)
      model_train <- coxph(Surv(time, status) ~ lp, data = data_train, y = TRUE, x = TRUE)

      lp_test <- as.vector(survObj.new$x %*% beta_m)
      data_test <- data.frame(time = survObj.new$t, status = survObj.new$di, lp = lp_test)

      if (any(type %in% c("hazard", "cumhazard", "survival"))) {
        Brier <- predictCox(model_train,
          times = times, type = type,
          newdata = data.frame(survObj.new$x)
        )
        Brier <- data.frame(times = fit$times, fit[names(fit) %in% type])
      } else {
        Brier <- Score(list("Bayesian Cox" = model_train),
          formula = Surv(time, status) ~ 1,
          data = data_test, conf.int = FALSE,
          metrics = "brier", summary = "ibs",
          times = times
        )$Brier$score
        # Brier <- BrierScore[BrierScore$model != "Null model", ]
        # extract IBS for Null model and the Bayesian Cox model
        ibs[1, ] <- Brier$IBS[c(length(times), length(times) * 2)]
      }
    } else {
      # obtain linear predictors corresponding to MCMC estimates
      lp_all_train <- survObj_train[, -c(1:2)] %*% t(betas)
      lp_all_test <- survObj.new$x %*% t(betas)
      Brier <- matrix(0, nrow = length(times) * 2, ncol = 2)
      data_train <- data.frame(
        time = survObj_train[, 1],
        status = survObj_train[, 2],
        lp = lp_all_train[, 1]
      )
      data_test <- data.frame(
        time = survObj.new$t,
        status = survObj.new$di,
        lp = lp_all_test[, 1]
      )
      model_train <- coxph(Surv(time, status) ~ lp,
        data = data_train,
        y = TRUE, x = TRUE
      )
      # calculate Brier scores based on the 1st MCMC estimates
      BrierScore0 <- Score(list("Bayesian Cox" = model_train),
        formula = Surv(time, status) ~ 1,
        data = data_test, conf.int = FALSE,
        metrics = "brier", summary = "ibs",
        times = times
      )$Brier$score
      Brier <- as.matrix(BrierScore0[1:(nrow(BrierScore0) / 2), -c(1:2)])
      # calculate Brier scores based on other MCMC estimates
      for (i in 2:nrow(betas)) {
        data_train$lp <- lp_all_train[, i]
        data_test$lp <- lp_all_test[, i]
        model_train <- coxph(Surv(time, status) ~ lp,
          data = data_train,
          y = TRUE, x = TRUE
        )
        BrierScore <- Score(list("Bayesian Cox" = model_train),
          formula = Surv(time, status) ~ 1,
          data = data_test, conf.int = FALSE,
          metrics = "brier", summary = "ibs",
          null.model = FALSE, times = times
        )$Brier$score
        Brier <- Brier + as.matrix(BrierScore[, -c(1:2)])
      }
      BrierScore0 <- data.frame(BrierScore0)
      BrierScore0[-c(1:(nrow(BrierScore0) / 2)), 3:4] <- Brier / nrow(betas)
      Brier <- BrierScore0
      # extract IBS for Null model and the Bayesian Cox model
      ibs[1, ] <- Brier$IBS[c(length(times), length(times) * 2)]
    }
    rownames(ibs) <- "IBS"
    print(ibs)
    invisible(Brier)
  }

  # invisible(fit)
}
