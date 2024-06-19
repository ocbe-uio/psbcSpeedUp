#' psbcSpeedUp
#' @title Function to Fit the Bayesian Cox Lasso Model
#'
#' @description
#' This a speed-up and extended version of the function \code{psbcGL()} in the R package \code{psbcGrouup}
#'
#' @name psbcSpeedUp
#' @useDynLib psbcSpeedUp
#' @aliases psbcSpeedUp-package
#' @importFrom Rcpp evalCpp
#' @importFrom xml2 as_xml_document write_xml
#' @importFrom stats rexp rgamma runif
#' @importFrom utils write.table
#' @importFrom survival survreg Surv
#'
#' @param survObj a list containing observed data from \code{n} subjects;
#' \code{t}, \code{di}, \code{x}. See details for more information
#' @param p number of covariates for variable selection
#' @param q number of mandatory covariates
#' @param hyperpar a list containing prior parameter values; among
#' \code{c('groupInd', 'beta.ini', 'eta0', 'kappa0', 'c0', 'r', 'delta',
#' 'lambdaSq', 'sigmaSq', 'tauSq', 's', 'h', 'beta.prop.var',
#' 'beta.clin.var')}. See details for more information
#' @param nIter the number of iterations of the chain
#' @param burnin number of iterations to discard at the start of the chain.
#' Default is 0
#' @param thin thinning MCMC intermediate results to be stored
#' @param rw when setting to "TRUE", the conventional random walk Metropolis
#' Hastings algorithm is used. Otherwise, the mean and the variance of the
#' proposal density is updated using the jumping rule described in
#' Lee et al. (2011)
#' @param outFilePath path to where the output files are to be written
#' @param tmpFolder the path to a temporary folder where intermediate data
#' files are stored (will be erased at the end of the chain). It is specified
#' relative to \code{outFilePath}
#'
#' @details
#' \tabular{ll}{
#' \code{t} \tab a vector of \code{n} times to the event \cr
#' \code{di} \tab a vector of \code{n} censoring indicators for the event time (1=event occurred, 0=censored) \cr
#' \code{x} \tab covariate matrix, \code{n} observations by \code{p} variables\cr
#' \code{groupInd} \tab a vector of \code{p} group indicator for each variable\cr
#' \code{beta.ini} \tab the starting values for coefficients \eqn{\beta}\cr
#' \code{eta0} \tab scale parameter of gamma process prior for the cumulative baseline hazard, \eqn{\eta_0 > 0}\cr
#' \code{kappa0} \tab shape parameter of gamma process prior for the cumulative baseline hazard, \eqn{\kappa_0 > 0}\cr
#' \code{c0} \tab the confidence parameter of gamma process prior for the cumulative baseline hazard, \eqn{c_0 > 0}\cr
#' \code{r} \tab the shape parameter of the gamma prior for \eqn{\lambda^2}\cr
#' \code{delta} \tab the rate parameter of the gamma prior for \eqn{\lambda^2}\cr
#' \code{lambdaSq} \tab the starting value for \eqn{\lambda^2}\cr
#' \code{sigmaSq} \tab the starting value for \eqn{\sigma^2}\cr
#' \code{tauSq} \tab the starting values for \eqn{\tau^2}\cr
#' \code{s} \tab the set of time partitions for specification of the cumulative baseline hazard function\cr
#' \code{h} \tab the starting values for \eqn{h}\cr
#' \code{beta.prop.var} \tab the variance of the proposal density for \eqn{\beta} in a random walk M-H sampler\cr
#' \code{beta.clin.var} \tab the starting value for the variance of \eqn{\beta}\cr
#' }
#'
#' @return An object of class \code{psbcSpeedUp} is saved as
#' \code{obj_psbcSpeedUp.rda} in the output file, including the following components:
#' \itemize{
#' \item input - a list of all input parameters by the user
#' \item output - a list of the all output estimates:
#' \itemize{
#' \item "\code{beta.p}" - a matrix with MCMC intermediate estimates of the regression coefficients.
#' \item "\code{h.p}" - a matrix with MCMC intermediate estimates of the increments in the cumulative baseline hazard in each interval.
#' \item "\code{tauSq.p}" - a vector MCMC intermediate estimates of the hyperparameter "tauSq".
#' \item "\code{sigmaSq.p}" - a vector MCMC intermediate estimates of the hyperparameter "sigmaSq".
#' \item "\code{lambdaSq.p}" - a vector MCMC intermediate estimates of the hyperparameter "lambdaSq".
#' \item "\code{accept.rate}" - a vector acceptance rates of individual regression coefficients.
#' }
#' \item call - the matched call.
#' }
#'
#'
#' @references Lee KH, Chakraborty S, and Sun J (2011). Bayesian Variable
#' Selection in Semiparametric Proportional Hazards Model for High Dimensional
#' Survival Data. \emph{The International Journal of Biostatistics}, 7(1):1-32.
#' @references Zucknick M, Saadati M, and Benner A (2015). Nonidentical twins:
#' Comparison of frequentist and Bayesian lasso for Cox models.
#' \emph{Biometrical Journal}, 57(6):959–81.
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
#'   "beta.prop.var" = 1, "beta.clin.var" = 1
#' )
#'
#' \donttest{
#' # run Bayesian Lasso Cox
#' library("psbcSpeedUp")
#' set.seed(123)
#' fitBayesCox <- psbcSpeedUp(survObj,
#'   p = p, q = q, hyperpar = mypriorPara,
#'   nIter = 10, burnin = 0, outFilePath = tempdir()
#' )
#' plot(fitBayesCox, color = "blue")
#' }
#'
#' @export
psbcSpeedUp <- function(survObj = NULL, p = 0, q = 0, hyperpar = list(),
                        nIter = 1, burnin = 0, thin = 1, rw = FALSE,
                        outFilePath, tmpFolder = "tmp/") {
  # Check the survival input object
  if (!is.list(survObj)) {
    stop("Argument 'survObj' must be an input")
  }
  if (sum(names(survObj) %in% c("t", "di", "x")) != 3) {
    stop("List 'survObj' must have three compoents 't', 'di' and 'x'!")
  }
  if (is.data.frame(survObj$x) || is.matrix(survObj$x)) {
    if (is.data.frame(survObj$x)) {
      survObj$x <- data.matrix(survObj$x)
    }
  } else {
    stop("Data 'survObj$x' must be a dataframe or a matrix!")
  }

  if (p + q == 0) {
    p <- ncol(survObj$x)
  } else {
    if (p %% 1 != 0 || p < 1) {
      stop("Argument 'p' must be a positive integer!")
    }
    if (q %% 1 != 0 || q < 0) {
      stop("Argument 'q' must be a positive integer!")
    }
    if (p + q != ncol(survObj$x)) {
      stop("The sum of 'p' and 'q' must equal the number of columns of 'survObj$x'!")
    }
  }
  if (thin %% 1 != 0 || thin < 1) {
    stop("Argument 'thin' must be a positive integer!")
  }

  # Check the directory for the output files
  if (outFilePath == "") {
    stop("Please specify a directory to save all output files!")
  }

  outFilePathLength <- nchar(outFilePath)
  if (substr(outFilePath, outFilePathLength, outFilePathLength) != "/") {
    outFilePath <- paste0(outFilePath, "/")
  }
  if (!file.exists(outFilePath)) {
    dir.create(outFilePath)
  }

  # Create temporary directory
  tmpFolderLength <- nchar(tmpFolder)
  if (substr(tmpFolder, tmpFolderLength, tmpFolderLength) != "/") {
    tmpFolder <- paste0(tmpFolder, "/")
  }
  tmpFolder <- paste0(outFilePath, tmpFolder)
  if (!file.exists(tmpFolder)) {
    dir.create(tmpFolder)
  }

  # check the formula
  cl <- match.call()

  # Write down in a single data file
  write.table(cbind(survObj$t, survObj$di, survObj$x),
    paste0(tmpFolder, "data.txt"),
    row.names = FALSE, col.names = FALSE
  )
  data <- paste0(tmpFolder, "data.txt")

  # Check hyperparameters
  if (!is.list(hyperpar)) {
    stop("Argument 'hyperpar' must be a list!")
  }
  if (any(!names(hyperpar) %in%
    c(
      "beta.ini", "lambdaSq", "sigmaSq", "rate", "tauSq",
      "s", "h", "groupInd", "eta0", "kappa0", "c0", "r", "delta",
      "beta.prop.var", "beta.clin.var"
    ))) {
    stop("Hyperparameters must be among c('beta.ini', 'lambdaSq', 'sigmaSq',
         'rate', 'tauSq', 'h', 'groupInd', 'eta0', 'kappa0', 'c0', 'r',
         'delta', 'beta.prop.var', 'beta.clin.var')!")
  }
  if ("groupInd" %in% names(hyperpar)) {
    if (length(hyperpar$groupInd) != p) {
      stop("Please specify correct hyperpar$groupInd!")
    }
    groupInd <- hyperpar$groupInd
  } else {
    groupInd <- 1:p
  }
  if ("beta.ini" %in% names(hyperpar)) {
    ini_beta <- hyperpar$beta.ini
  } else {
    ini_beta <- rep(0, p + q)
  }
  if (!"lambdaSq" %in% names(hyperpar)) {
    hyperpar$lambdaSq <- 1
  }
  if (!"sigmaSq" %in% names(hyperpar)) {
    hyperpar$sigmaSq <- runif(1, 0.1, 10)
  }
  if (!"rate" %in% names(hyperpar)) {
    hyperpar$rate <- hyperpar$lambdaSq / 2
  }
  if ("tauSq" %in% names(hyperpar)) {
    ini_tauSq <- hyperpar$tauSq
  } else {
    ini_tauSq <- rexp(length(unique(groupInd))) # , rate = hyperpar$lambdaSq/2)
  }
  s0 <- survObj$t[survObj$di == 1]
  if ("s" %in% names(hyperpar)) {
    if (min(s) > min(s0) | max(s) < max(s0)) {
      stop("Parameter 's' does not cover all event times!")
    }
    if (any(diff(s) < 0)) {
      s <- sort(s)
    }
    if (any(duplicated(s))) {
      s <- unique(s)
    }
  } else {
    s <- sort(unique(s0)) # event times that are actually observed
    s <- c(s, 2 * max(survObj$t) - max(survObj$t[-which.max(survObj$t)]))
  }
  if ("h" %in% names(hyperpar)) {
    ini_h <- hyperpar$h
  } else {
    ini_h <- rgamma(length(s), 1, 1)
  }
  fit <- survreg(Surv(survObj$t, survObj$di, type = c('right')) ~ 1, 
                 dist = "weibull", x = TRUE, y = TRUE)
  if (!"kappa0" %in% names(hyperpar)) {
    hyperpar$kappa0  <- 1 / exp(fit$icoef["Log(scale)"])
  }
  if (!"eta0" %in% names(hyperpar)) {
    # hyperpar$eta0 <- round(log(2) / 36, 2)
    hyperpar$eta0 <- exp(fit$coefficients)^(-hyperpar$kappa0) 
  }
  if (!"c0" %in% names(hyperpar)) {
    hyperpar$c0 <- 2
  }
  if (!"r" %in% names(hyperpar)) {
    hyperpar$r <- 10 / 9
  }
  if (!"delta" %in% names(hyperpar)) {
    hyperpar$delta <- 0.00001
  }
  if ("beta.prop.var" %in% names(hyperpar)) {
    hyperpar$beta_prop_var <- hyperpar$beta.prop.var
  } else {
    hyperpar$beta_prop_var <- 1
  }
  if ("beta.clin.var" %in% names(hyperpar)) {
    hyperpar$beta_clin_var <- hyperpar$beta.clin.var
  } else {
    hyperpar$beta_clin_var <- 1
  }

  hyperParFile <- paste0(tmpFolder, "hyperpar.xml")

  ## Create the return object
  ret <- list(input = list(), output = list(), call = cl)
  class(ret) <- "psbcSpeedUp"

  # Copy the inputs
  ret$input["p"] <- p
  ret$input["q"] <- q
  ret$input["nIter"] <- nIter
  ret$input["burnin"] <- burnin
  ret$input["thin"] <- thin
  ret$input["rw"] <- rw
  ret$input["outFilePath"] <- outFilePath
  ret$input["tmpFolder"] <- tmpFolder
  ret$input["hyperpar"] <- list(hyperpar)

  hyperpar$s <- hyperpar$beta.ini <- hyperpar$tauSq <- hyperpar$h <- hyperpar$groupInd <-
    hyperpar$beta.prop.var <- hyperpar$beta.clin.var <- NULL

  ## Set up the XML file for hyperparameters
  xml <- xml2::as_xml_document(
    list(hyperparameters = list(
      lapply(hyperpar, function(x) list(format(x, scientific = FALSE))) # every element in the list should be a list
    ))
  )
  xml2::write_xml(xml, file = hyperParFile)

  # Run Bayesian Cox model
  nChains <- 1
  ret$output <- psbcSpeedUp_internal(
    data, p, q, hyperParFile, outFilePath,
    ini_beta, ini_tauSq, ini_h, groupInd, # hyperparameters which are vectors
    nIter, nChains, thin, rw
  )

  ret$output$accept.rate <- as.vector(ret$output$accept.rate) / nIter

  if (is.null(colnames(survObj$x))) {
    colnames(ret$output$beta.p) <- paste0("x", seq_len(ncol(survObj$x)))
  } else {
    colnames(ret$output$beta.p) <- colnames(survObj$x)
  }

  ## Save fitted object
  obj_psbc <- list(input = ret$input, output = ret$output)
  save(obj_psbc, file = paste0(outFilePath, "obj_psbc.rda"))

  if (outFilePath != tmpFolder) {
    unlink(tmpFolder, recursive = TRUE)

    write.table(cbind(survObj$t, survObj$di, survObj$x),
      paste0(outFilePath, "data.txt"),
      row.names = FALSE, col.names = FALSE
    )
  }

  return(ret)
}
