##' Reduced-rank mixture models with optimal tuning parameter(s)
##'
##' Reduced-rank mixture models with optimal tuning parameter(s)
##'
#' @usage
#' tune.rrmix(K = NULL, K.max = NULL, X, Y, est = c("FR", "RP", "ANNP"),
#'            lambda = NULL, n.lambda = 20, gamma = 2,
#'            ind0 = NULL, para0 = NULL, seed = NULL, kmscale = FALSE, km.nstart = 20,
#'            n.init = 100, commonvar = FALSE, maxiter = 1000, maxiter.int = 100,
#'            thres = 1e-05, thres.int = 1e-05, 
#'            para.true = NULL, ind.true = NULL)
#'
#' @param K number of mixture components. Required when K.max is `NULL'.
#' @param K.max maximum of mixture components. Default is `NULL'. When provided,
#'   the argument K is ignored.
#' @param X n by p design matrix where n is the number of observations and
#'   p is the number of predictors.
#' @param Y n by q response matrix where n is the number of observations and
#'   q is the number of responses.
#' @param est character, specifying the estimation method. `FR', `RP', and `ANNP'
#'   refers to as the full-ranked, rank-penalized, and adaptive nuclear norm
#'   penalized method, respectively.
#' @param lambda vector consisting of lambda candidates. Only used in the
#'   estimation method of `RP' and `ANNP'. If 0, all estimation methods (`FR',
#'   `RP', and `ANNP') provide the same estimation results. Default is 'NULL'.
#'   If 'NULL', data-adaptive range of lambda will be provided internally.
#' @param n.lambda number of lambda candidates to explore. Only used when
#'   'lambda' is 'NULL'. Default is 20.
#' @param gamma numerical value, specifying additional tuning parameter, only
#'   used in the estimation method of `ANNP'. It must be nonnegative.
#' @param ind0 vector of length n, specifying the initial assignment of the
#'   mixture membership of n observations when there is prior information on
#'   the membership. If `NULL', K-means clustering technique is used to assign
#'   the membership for n observations. Default is `NULL'.
#' @param para0 array of length K. It consists of K lists, each of which contains
#'   initial values of membership probability, coefficient matrix, and variance-
#'   covariance matrix.
#' @param seed seed number for the reproducibility of results. Default of `NULL'.
#' @param km.nstart number of random sets considered to perform K-means
#'   clustering for initialization. Default is 20.
#' @param kmscale logical value, indicating whether Y is scaled prior to K-means
#'   clustering for initialization. Default is `FALSE'.
#' @param n.init number of initializations to try. Two methods for initial
#'   clustering are used: K-means and random clustering.
#' @param commonvar logical value, indicating the homogeneity assumption of
#'   variance-covariance matrices across K mixture components. Default is `FALSE'.
#' @param maxiter maximum number of iterations for external iterative
#'   algorithm, used in all estimation methods.
#' @param maxiter.int maximum number of iterations for internal iterative
#'   algorithm, only used in the estimation method of `ANNP'.
#' @param thres threshold value for external EM algorithm, used in all
#'  estimation methods. It controls the termination of the EM algorithm.
#' @param thres.int threshold value for internal iterative algorithm, only used
#'   in the estimation method of `ANNP'. It controls the termination of the
#'  internal algorithm. 
#' @param ind.true vector of length n, specifying the true mixture membership
#'   for n observations. Only used when true models are known, e.g., in a
#'   simulation study.
#' @param para.true array of length K. It consists of K lists, each of which
#'   contains a coefficient matrix and its true rank. Only used when true models
#'   are known, e.g., in a simulation study.
#'   
#' @return
#'
#'   \item{lambda.cand}{lambda values used as input.}
#'   \item{penloglik}{penalized log-likelihood values corresponding to the set
#'     of lambda values.}
#'   \item{bic}{BIC values corresponding to the set of lambda values.}
#'   \item{est.rank}{estimated ranks corresponding to the set of lambda values.}
#'
#' @author
#'   Suyeon Kang, University of California, Riverside, \email{skang062@@ucr.edu};
#'   Weixin Yao, University of California, Riverside, \email{weixin.yao@@ucr.edu};
#'   Kun Chen, University of Connecticut, \email{kun.chen@@uconn.edu}.
#'
#' @seealso \code{\link{rrmix}}
#'
#' @references
#'   Kang, S., Chen, K., and Yao, W. (2022+). "Reduced rank estimation in mixtures
#'   of multivariate linear regression".
#'   
#' @export
#' @examples
#' #-----------------------------------------------------------#
#' # Real Data Example: Tuna Data
#' #-----------------------------------------------------------#
#' require(bayesm)
#' data(tuna)
#' tunaY <- log(tuna[, c("MOVE1", "MOVE2", "MOVE3", "MOVE4",
#'                   "MOVE5", "MOVE6", "MOVE7")])
#' tunaX <- tuna[, c("NSALE1", "NSALE2", "NSALE3", "NSALE4",
#'               "NSALE5", "NSALE6", "NSALE7",
#'               "LPRICE1", "LPRICE2", "LPRICE3", "LPRICE4",
#'               "LPRICE5", "LPRICE6", "LPRICE7")]
#' tunaX <- cbind(intercept = 1, tunaX)
#' 
#' \donttest{
#' tuna.tune <- tune.rrmix(K.max = 3, X = tunaX, Y = tunaY, est = "RP",
#'              lambda = exp(seq(0, log(100), length = 20)),
#'              seed = 100, n.init = 100)
#' summary(tuna.tune)
#' plot(tuna.tune, transform.y = log, ylab = "log(lambda)")} 
tune.rrmix <- function (K = NULL,
                        K.max = NULL,
                        X,
                        Y,
                        est = c("FR", "RP", "ANNP"),
                        lambda = NULL,
                        n.lambda = 20,
                        gamma = 2,
                        ind0 = NULL,
                        para0 = NULL,
                        seed = NULL,
                        kmscale = FALSE,
                        km.nstart = 20,
                        n.init = 100,
                        commonvar = FALSE,
                        maxiter = 1000,
                        maxiter.int = 100,
                        thres = 1e-05,
                        thres.int = 1e-05,
                        para.true = NULL,
                        ind.true = NULL)
{
  
  # Match arguments
  call <- match.call()
  est <- match.arg(est)
  
  # Warning if inappropriate setting
  if(!is.null(K.max)){
    if(!is.null(para.true) || !is.null(ind.true)){
      stop("'K.max' must be NULL when true parameters are known")
    } else {
      if(!is.null(K)){
        message("'K' is ignored when 'K.max' is provided")
      }
      K <- 1:K.max
      #para0 <- ind0 <- call$para0 <- call$ind0 <- NULL
    }
    if(!is.null(para0) || !is.null(ind0)){
      para0 <- ind0 <- call$para0 <- call$ind0 <- NULL
      message("'para0' and/or 'ind0' is ignored when 'K.max' is provided")
    }
  }
  if(nrow(X) != nrow(Y)){
    stop("number of rows of X and Y must be match.")
  }
  if(nrow(X) < ncol(X) | nrow(Y) < ncol(Y)){
    stop("currently the case when n < min(p, q) is not supported.")
  }
  
  if(sum(X[, 1] == 1) == nrow(X)){  # if intercept column exists in X
    colnames(X)[1] <- "intercept"
  } else {
    message("Intercept column added to X.")
    X <- as.matrix(cbind(intercept = 1, X))
  }
  
  if(sum(gamma < 0) > 0){
    stop("gamma must be nonnegative.")
  }
  
  # When lambda = NULL, provide data-adaptive lambda candidates
  if(is.null(lambda)){
    if(est == "FR"){
      lambda <- 0
    } else if (est == "RP" | est == "ANNP"){
      X <- as.matrix(X); Y <- as.matrix(Y)
      PY <- X %*% MASS::ginv(t(X) %*% X) %*% (t(X) %*% Y)
      n <- nrow(Y); q <- ncol(Y)
      S2 <- sum(svd(Y - PY)$d[1:min(n, q)]^2)/(n * q - as.numeric(Matrix::rankMatrix(X)) * q)
      muadapt <- 2 * S2 * (q + as.numeric(Matrix::rankMatrix(X)))
      if(est == "RP"){
        lambmax <- sqrt(muadapt * 2) * 4
      } else if (est == "ANNP"){
        lambmax <- muadapt * 4
      }
      lambda <- sort(seq(1, lambmax, length = n.lambda))
    }
  }
  
  if(est == "FR" || est == "RP"){
    ranges <- list(K = K, lambda = lambda, gamma = NA)
  } else if (est == "ANNP"){
    ranges <- list(K = K, lambda = lambda, gamma = gamma)
  }
  
  
  ## find best model
  parameters <- expand.grid(ranges)
  parameters <- parameters[order(parameters[, 1]), ]
  np <- nrow(parameters)
  plls <- bics <- scls <- hcls <- ees <- pes <- rep(NA, np)
  lls <- ps <- npars <- rep(NA, np)
  ranks <- matrix(NA, np, max(K))
  
  ## - loop over all models
  #if(!is.null(seed)){
  #  set.seed(seed)
  #}
  
  cat("\nNumber of parameter combinations to examine: ", np, "\n")
  
  if(is.null(K.max)){
    para0 <- initialize.para(K = K, X = X, Y = Y, ind0 = ind0,
                             seed = seed, km.nstart = km.nstart,
                             kmscale = kmscale, n.init = n.init, commonvar = commonvar)
  } else {
    para0all <- array(list(), K.max)
    for(i in 1:K.max){
      para0all[[i]] <- initialize.para(K = i, X = X, Y = Y, ind0 = ind0,
                                       seed = seed, km.nstart = km.nstart,
                                       kmscale = kmscale, n.init = n.init, commonvar = commonvar)
    }
  }
  
  for(i in 1:np){
    if(is.na(parameters[i, ]$gamma)){
      gamma <- NULL
    } else {
      gamma <- parameters[i, ]$gamma
    }
    
    if(!is.null(K.max)){
      if((parameters[i, ]$K - 1) * np/K.max < i & i <= parameters[i, ]$K * np/K.max){
        para0 <- para0all[[parameters[i, ]$K]]
      }
    }
    
    mod <- suppressMessages(rrmix(K = parameters[i, ]$K, X = X, Y = Y, est = est,
                 lambda = parameters[i, ]$lambda, gamma = gamma,
                 ind0 = ind0, para0 = para0, seed = seed,
                 kmscale = kmscale, km.nstart = km.nstart, n.init = n.init,
                 commonvar = commonvar, maxiter = maxiter,
                 maxiter.int = maxiter.int,
                 thres = thres, thres.int = thres.int, 
                 para.true = para.true, ind.true = ind.true))
    lls[i] <- mod$loglik
    ps[i] <- mod$penalty
    npars[i] <- mod$npar
    plls[i] <- mod$penloglik
    bics[i] <- mod$bic
    if(!is.null(para.true) || !is.null(ind.true)){
      scls[i] <- mod$class.err[1]
      hcls[i] <- mod$class.err[2]
      ees[i] <- mean(mod$est.err)
      pes[i] <- mean(mod$pred.err)
    }
    ranks[i, 1:parameters[i, ]$K] <- mod$est.rank
    cat("\nexamined: ", i)
  }
  
  if(is.null(para.true)){
    if(est == "FR" || est == "RP"){
      tab <- cbind(parameters$K, parameters$lambda, lls, ps, plls, npars, bics, ranks)
      colnames(tab) <- c("K", "lambda", "loglik", "penalty", "penloglik", "npar", "bic",
                         paste("est.rank", 1:max(parameters$K), sep = ""))
    } else if (est == "ANNP"){
      tab <- cbind(parameters$K, parameters$lambda, parameters$gamma, lls, ps, plls, npars, bics, ranks)
      colnames(tab) <- c("K", "lambda", "gamma", "loglik", "penalty", "penloglik", "npar", "bic",
                         paste("est.rank", 1:max(parameters$K), sep = ""))
    }
  } else {
    if(est == "FR" || est == "RP"){
      tab <- cbind(parameters$K, parameters$lambda, lls, ps, plls, npars, bics,
                   scls, hcls, ees, pes, ranks)
      colnames(tab) <- c("K", "lambda", "loglik", "penalty", "penloglik", "npar", "bic",
                         "soft.class.err", "hard.class.err", "est.err", "pred.err",
                         paste("est.rank", 1:max(parameters$K), sep = ""))
    } else if (est == "ANNP"){
      tab <- cbind(parameters$K, parameters$lambda, parameters$gamma, lls, ps, plls, npars, bics,
                   scls, hcls, ees, pes, ranks)
      colnames(tab) <- c("K", "lambda", "gamma", "loglik", "penalty", "penloglik", "npar", "bic",
                         "soft.class.err", "hard.class.err", "est.err", "pred.err",
                         paste("est.rank", 1:max(parameters$K), sep = ""))
    }
  }
  rownames(tab) <- 1:nrow(tab)
  #tab <- tab[, -c("penalty", "npar")]
  
  if(is.null(para.true)){
    if(est == "FR" || est == "RP"){
      output <- list(call = call, est = est, K = parameters$K,
                     lambda = parameters$lambda,
                     penloglik = plls, bic = bics,
                     est.rank = ranks, details = tab)
    } else if (est == "ANNP"){
      output <- list(call = call, est = est, K = parameters$K,
                     lambda = parameters$lambda, gamma = parameters$gamma,
                     penloglik = plls, bic = bics,
                     est.rank = ranks, details = tab)
    }
  } else {
    if(est == "FR" || est == "RP"){
      output <- list(call = call, est = est, K = parameters$K,
                     lambda = parameters$lambda,
                     penloglik = plls, bic = bics, soft.class.err = scls,
                     hard.class.err = hcls, est.err = ees, pred.err = pes,
                     est.rank = ranks, details = tab)
    } else if (est == "ANNP"){
      output <- list(call = call, est = est, K = parameters$K,
                     lambda = parameters$lambda, gamma = parameters$gamma,
                     penloglik = plls, bic = bics, soft.class.err = scls,
                     hard.class.err = hcls, est.err = ees, pred.err = pes,
                     est.rank = ranks, details = tab)
    }
  }
  
  class(output) <- "tune.rrmix"
  return(output)
}

summary.tune.rrmix <- function(object, metric = c("bic", "soft.class.err",
                                                  "hard.class.err", "est.err",
                                                  "pred.err"), ...){
  
  # Match arguments
  metric <- match.arg(metric)
  
  # Warning for inappropriate settings
  if(is.null(object$call$para.true) && metric != "bic"){
    stop("'metric' must be 'bic' when true parameters are unknown")
  }
  
  best.ind <- as.numeric(which.min(object$details[, metric]))
  if(object$est == "FR" || object$est == "RP"){
    parameters <- c("K", "lambda")
  } else if (object$est == "ANNP"){
    parameters <- c("K", "lambda", "gamma")
  }
  tmp <- data.frame(matrix(object$details[best.ind, parameters], nrow = 1))
  colnames(tmp) <- parameters
  rownames(tmp) <- ""
  
  object <- c(object, metric = metric, best = tmp)
  structure(object, class = "summary.tune.rrmix")
  
  
}

plot.tune.rrmix <- function(x,
                            metric = c("bic", "soft.class.err",
                                       "hard.class.err", "est.err", "pred.err"),
                            col = "blue",
                            main = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            swapxy = FALSE,
                            transform.x = NULL,
                            transform.y = NULL,
                            transform.z = NULL,
                            color.palette = hsv_palette(),
                            nlevels = 20,
                            ...){
  
  # Match arguments
  metric <- match.arg(metric)
  
  if(x$est == "FR" || x$est == "RP"){
    parameters <- c("K", "lambda")
  } else if (x$est == "ANNP"){
    parameters <- c("K", "lambda", "gamma")
  }
  
  tab <- data.frame(x$details[, parameters])
  tab <- tab[vapply(tab, function(x) length(unique(x)) > 1, logical(1L))]
  k <- ncol(tab)
  if (k > 2) stop("Cannot visualize more than 2 parameters")
  if (is.null(main)){
    main <- paste("Performance of `",
                  switch(x$est, "FR" = "full ranked mixture",
                         "RP" = "rank penalized mixture",
                         "ANNP" = "adaptive nuclear norm penalized mixture"),
                  "'", sep = "")
  }
  tab <- cbind(tab, metric = x$details[, metric])
  
  if (k == 1){
    if (!is.null(transform.x))
      tab[, 1] <- transform.x(tab[, 1])
    if (!is.null(transform.y))
      tab[, 2] <- transform.y(tab[, 2])
    if (is.null(xlab)) xlab <- colnames(tab)[1]
    if (is.null(ylab)) ylab <- metric
    plot(tab[, 1], tab[, 2], type = "b", main = main, xlab = xlab, ylab = ylab)
    graphics::abline(v = tab[which.min(tab[, 2]), 1], col = "blue", lty = 2)
  } else if (k == 2){
    if (!is.null(transform.x))
      tab[, 1] <- transform.x(tab[, 1])
    if (!is.null(transform.y))
      tab[, 2] <- transform.y(tab[, 2])
    if (swapxy)
      tab[, 1:2] <- tab[, 2:1]
    xx <- stats::xtabs(metric ~ ., data = tab)
    if (is.null(xlab)) xlab <- names(dimnames(xx))[1 + swapxy]
    if (is.null(ylab)) ylab <- names(dimnames(xx))[2 - swapxy]
    if (!is.null(transform.z))
      xx <- transform.z(xx)
    
    ## plotting
    graphics::filled.contour(x = as.double(rownames(xx)),
                             y = as.double(colnames(xx)),
                             z = xx,
                             xlab = xlab,
                             ylab = ylab,
                             nlevels = nlevels,
                             color.palette = color.palette,
                             main = main)
  }
  
}