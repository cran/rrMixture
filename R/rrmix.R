#' Reduced-Rank Mixture Models in Multivariate Regression
#'
#' `rrmix' is used to estimate parameters of reduced-rank mixture models in
#' multivariate linear regression using the full-ranked, rank-penalized, and
#' adaptive nuclear norm penalized estimators proposed by Kang et. al. (2022+).
#'
#' @usage
#' rrmix(K = 2, X, Y, est = c("FR", "RP", "ANNP"),
#'       lambda = 0, gamma = 2, ind0 = NULL, para0 = NULL, seed = NULL,
#'       kmscale = FALSE, km.nstart = 20, n.init = 100, commonvar = FALSE,
#'       maxiter = 1000, maxiter.int = 100, thres = 1e-05, thres.int = 1e-05,
#'       visible = FALSE, para.true = NULL, ind.true = NULL)       
#'
#' @param K number of mixture components.  
#' @param X n by p design matrix where n is the number of observations and
#'   p is the number of predictors.
#' @param Y n by q response matrix where n is the number of observations and
#'   q is the number of responses.
#' @param est character, specifying the estimation method. `FR', `RP', and `ANNP'
#'   refers to as the full-ranked, rank-penalized, and adaptive nuclear norm
#'   penalized method, respectively.
#' @param lambda numerical value, specifying tuning parameter. Only used in the
#'   estimation method of `RP' and `ANNP'. If 0, all estimation methods (`FR',
#'   `RP', and `ANNP') provide the same estimation results.
#' @param gamma numerical value, specifying additional tuning parameter, only
#'   used in the estimation method of `ANNP'. It must be nonnegative.
#' @param ind0 vector of length n, specifying the initial assignment of the
#'   mixture membership of n observations when there is prior information on
#'   the membership. If `NULL', K-means clustering technique is used to assign
#'   the membership for n observations. Default is `NULL'.
#' @param para0 array of length K. It consists of K lists, each of which contains
#'   initial values of membership probability, coefficient matrix, and variance-
#'   covariance matrix.
#' @param seed seed number for the reproducibility of initialization results in
#'   the EM algorithm. Default is `NULL'.
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
#' @param visible logical value, indicating whether the outputs from each iteration
#'   are printed. Useful when the whole algorithm takes long. Default is `FALSE'.
#' @param ind.true vector of length n, specifying the true mixture membership
#'   for n observations. Only used when true models are known, e.g., in a
#'   simulation study.
#' @param para.true array of length K. It consists of K lists, each of which
#'   contains a coefficient matrix and its true rank. Only used when true models
#'   are known, e.g., in a simulation study.
#'   
#' @return
#' 
#'   An object of class \code{rrmix} containing the fitted model, including: 
#'   \item{call}{original function call.}
#'   \item{seed}{seed number which is set for the initilization.}
#'   \item{n.est}{vector of length K, specifying the estimated number of
#'     observations in each mixture components.}
#'   \item{para}{array of length K. It consists of K lists, each of which contains
#'     final estimates of membership probability, coefficient matrix, and variance-
#'     covariance matrix.}
#'   \item{est.rank}{vector of length K, specifying the estimated ranks of
#'     coefficient matrices.}
#'   \item{npar}{number of parameters in the model, used to estimate the BIC.}
#'   \item{n.iter}{number of iterations (external EM algorithm).}
#'   \item{lambda}{tuning parameter for the estimation method of 'RP' or 'ANNP'.}
#'   \item{gamma}{tuning parameter for the estimation method of 'ANNP'.}
#'   \item{ind}{vector of length n, specifying the estimated mixture membership
#'     for n observations.}
#'   \item{ind.true}{vector of length n, specifying the true mixture membership
#'     for n observations. Only returned when the true models are known.}
#'   \item{loglik}{log-likelihood of the final model.}
#'   \item{penloglik}{penalized log-likelihood of the final model.}
#'   \item{penalty}{penalty in the penalized log-likelihood of the final model.}
#'   \item{bic}{BIC of the final model.}
#'   \item{avg.nn.iter}{average number of iterations for internal iterative
#'     algorithm, only returned for the estimation method of 'ANNP'.}
#'   \item{resmat}{matrix containing the information for each iteration of the
#'     EM algorithm, e.g., iteration number, log-likelihood, penalized log-
#'     likelihood, difference between penalized log-likelihood values from two
#'     consecutive iterations, and computing time.}   
#'   \item{class.err}{Soft and hard classification errors for mixture membership.
#'     Only returned when the true models are known.} 
#'   \item{est.err}{estimation error from the comparison between the estimated
#'     and true coefficient matrices. Only returned when the true models are
#'     known.} 
#'   \item{pred.err}{prediction error. Only returned when the true models are
#'   known.} 
#'
#' @author
#'   Suyeon Kang, University of California, Riverside, \email{skang062@@ucr.edu};
#'   Weixin Yao, University of California, Riverside, \email{weixin.yao@@ucr.edu};
#'   Kun Chen, University of Connecticut, \email{kun.chen@@uconn.edu}.
#'
#' @seealso \code{\link{rrmix.sim.norm}}, \code{\link{initialize.para}} 
#'
#' @references
#'   Kang, S., Chen, K., and Yao, W. (2022+). "Reduced rank estimation in mixtures
#'   of multivariate linear regression".
#'
#' @export
#' @examples
#' library(rrMixture)
#' 
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
#' # Rank-penalized estimation
#' \donttest{
#' tuna.rp <- rrmix(K = 2, X = tunaX, Y = tunaY, lambda = 3, est = "RP",
#'            seed = 100, n.init = 100)
#' summary(tuna.rp)
#' plot(tuna.rp)} 
#' 
#' # Adaptive nuclear norm penalized estimation
#' \donttest{
#' tuna.annp <- rrmix(K = 2, X = tunaX, Y = tunaY, lambda = 3, gamma = 2, est = "ANNP",
#'              seed = 100, n.init = 100)
#' summary(tuna.annp)
#' plot(tuna.annp)}       
#' 
#' #-----------------------------------------------------------#
#' # Simulation: Two Components Case
#' #-----------------------------------------------------------#
#' # Simulation Data
#' K2mod <- rrmix.sim.norm(K = 2, n = 100, p = 5, q = 5, rho = .5,
#'          b = 1, shift = 1, r.star = c(1, 3), sigma = c(1, 1),
#'          pr = c(.5, .5), seed = 1215)
#'          
#' # Rank-penalized estimation
#' \donttest{
#' K2.rp <- rrmix(K = 2, X = K2mod$X, Y = K2mod$Y, lambda = 1,
#'          seed = 17, est = "RP", ind.true = K2mod$ind.true,
#'          para.true = K2mod$para.true, n.init = 100)
#' summary(K2.rp)
#' plot(K2.rp)}
#'          
#' # Adaptive nuclear norm penalized estimation
#' \donttest{
#' K2.annp <- rrmix(K = 2, X = K2mod$X, Y = K2mod$Y, lambda = 1,
#'            seed = 17, est = "ANNP", ind.true = K2mod$ind.true,
#'            para.true = K2mod$para.true, n.init = 100)
#' summary(K2.annp)
#' plot(K2.annp)}
rrmix <- function(K = 2,
                  X,
                  Y,
                  est = c("FR", "RP", "ANNP"),
                  lambda = 0,
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
                  visible = FALSE,
                  para.true = NULL,
                  ind.true = NULL)
{
  # Match arguments
  call <- match.call()
  est <- match.arg(est)
  force <- T
  X <- as.matrix(X); Y <- as.matrix(Y)
  
  # Warning messages
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
  
  if(!is.null(gamma)){
    if(gamma < 0){
      stop("gamma must be nonnegative.")
    }
  }
    
  # To measure the computing time
  start_time = Sys.time()
  
  # If para.true = NULL, NULL
  r.star <- unlist(lapply(para.true, '[[', "r"))
  
  if(est == "FR"){
    lambda <- NULL
  }
  if(est == "RP"|est == "FR"){
    gamma <- NULL
  }
  
  n <- nrow(Y); p <- ncol(X); q <- ncol(Y)
  SVD <- svd(X %*% MASS::ginv(t(X)%*% X) %*% (t(X)%*% Y))
  if(est == "ANNP"){
    wi <- SVD$d^(- gamma)
  }
  
  # If there is prior information for initial parameter estimates, use it.
  # Otherwise, use 'initialize.para' function in the package, which provides the MLEs.
  if(!is.null(para0)){
    para <- para0
  } else {
    para <- initialize.para(K = K, X = X, Y = Y, ind0 = ind0,
                            seed = seed, km.nstart = km.nstart,
                            kmscale = kmscale, n.init = n.init, commonvar = commonvar)
  }
  
  if(K == 1){
    
    ## K = 1: Estimate without EM algorithm
    if(est == "FR"){
      para[[1]]$B <- MASS::ginv(t(X) %*% X) %*% (t(X) %*% Y)
      penalty <- 0
      r.hat <- min(p, q)
    } else if (est == "RP"){
      res1 <- CH.11(Y = Y, X = X, lambda = lambda)
      para[[1]][["B"]] <- res1$B
      penalty <- res1$penalty
      r.hat <- ifelse(lambda != 0, res1$penalty/lambda^2, min(p, q))
    } else if (est == "ANNP"){
      res2 <- CS.13(Y = Y, X = X, lambda = lambda, wi = wi, getpen = T)
      para[[1]][["B"]] <- res2$B
      SVD <- svd(X %*% MASS::ginv(t(X) %*% X) %*% (t(X) %*% Y))
      r.hat <- sum(SVD$d[1:min(n, q)] - lambda^(1/(gamma + 1)) > 0)
      penalty <- g(Bk = para[[1]][["B"]], X = X, Y = Y, lambda = lambda, wi = wi)
    }
    para[[1]][["Sigma"]] <- (1/n) * t(Y - X %*% para[[1]]$B) %*% (Y - X %*% para[[1]]$B)
    para[[1]][["Sigma"]] <- as.symmetric.matrix(para[[1]][["Sigma"]])
    
    # results
    n.est <- n
    est.rank <- r.hat
    numc <- 1
    const <- 1
    n.iter <- NA
    ind <- ind.true <- rep(1, n)
    loglh <- sum(log(sapply(1:nrow(Y), temp.func, Y, X, para, K, pd.ind = 1)))
    new.loglik <- loglh - penalty
    npar <- (numc - 1) + const * q * (q + 1)/2 + est.rank * (p + q - est.rank)
    bic <- -2 * loglh + log(n) * npar
    sum.nn.iter <- avg.nn.iter <- NA
    resmat <- NA
    class.err <- rep(NA, 2)
    if(!is.null(para.true)){
      est.err <- 100 * sum(svd(para.true[[1]][["B"]] - para[[1]][["B"]])$d[1:min(p, q)]^2)/(p * q)
      pred.err <- 100 * sum(svd(X %*% para.true[[1]][["B"]] - X %*% para[[1]][["B"]])$d[1:min(n, q)]^2)/(n * q)
    }
    
  } else if (K >= 2){
    
    ## K >= 2: Estimate via EM algorithm
    n.iter <- 0; err <- -10 # any negative number
    loglh <- loglik(Y, X, para, K)
    if(est == "FR"){
      penalty <- 0
    } else if(est == "RP"){
      penalty <- (1/2) * lambda^2 * sum(sapply(X = 1:K, FUN = function(kk){
        Matrix::rankMatrix(para[[kk]][["B"]])}[1]))
    } else if(est == "ANNP"){
      penalty <- sum(sapply(X = 1:K, FUN = function(kk){
        g(Bk = para[[kk]][["B"]], X = X, Y = Y, lambda = lambda, wi = wi)}))
    }
    old.loglik <- loglh - penalty
    
    #lambda <- lambda * nrow(Y)
    p.ik <- matrix(NA, n, K)
    penalty <- 0
    r.hat <- rep(NA, K)
    sum.nn.iter <- avg.nn.iter <- 0
    resmat <- NULL
    resmat <- rbind(resmat, c(0, loglh, old.loglik, NA, NA, NA))
    colnames(resmat) <- c("t", "L", "F", "old.F_new.F", "avg.r", "CPU.time")
    
    while(abs(err) > thres && n.iter < maxiter){
      
      old.para <- para
      old.penalty <- penalty
      
      # E-Step
      for(kk in 1:K){ 
        para[[kk]][["Sigma"]] <- as.symmetric.matrix(para[[kk]][["Sigma"]])
      }
      pd.ind <- which(sapply(1:K, FUN = function(kk){is.valid.pd.mat(para[[kk]][["Sigma"]])}))
      res.temp <- sapply(1:n, temp.func, Y, X, para, K, pd.ind)
      for(kk in 1:K){
        p.ik[, kk] <- res.temp[kk, ]/colSums(res.temp)
      }
      
      # M-step
      if(est == "FR"){
        for(kk in 1:K){
          # update for pi_k
          para[[kk]][["pi"]] <- colSums(p.ik)[kk]/n
          # update for B_k
          Wk <- diag(p.ik[, kk])
          Yk.star <- sqrt(Wk) %*% Y; Xk.star <- sqrt(Wk) %*% X
          para[[kk]][["B"]] <- MASS::ginv(t(Xk.star) %*% Xk.star) %*% (t(Xk.star) %*% Yk.star)
          # update for Sigma_k
          if(!commonvar){
            para[[kk]][["Sigma"]] <- temp.func2(kk, Y, X, para, K, p.ik)/colSums(p.ik)[kk]
            para[[kk]][["Sigma"]] <- as.symmetric.matrix(para[[kk]][["Sigma"]])
          }
        }
        # update for Sigma if commonvar
        if(commonvar){
          eachSigma <- Reduce('+', lapply(1:K, temp.func2, Y, X, para, K, p.ik))/n
          eachSigma <- as.symmetric.matrix(eachSigma)
          for(kk in 1:K){
            para[[kk]][["Sigma"]] <- eachSigma
          }
        }
        penalty <- 0
        r.hat <- rep(min(p, q), K)
        
      } else if(est == "RP"){
        penalty <- rep(NA, K)
        for(kk in 1:K){
          # update for pi_k
          para[[kk]][["pi"]] <- colSums(p.ik)[kk]/n
          # update for B_k
          Wk <- diag(p.ik[, kk])
          if(!is.nan(para[[kk]][["Sigma"]][1])){
            svdsigma <- svd(para[[kk]][["Sigma"]])
            Sigmanegahalf <- svdsigma$v %*% diag(svdsigma$d^(-1/2)) %*% t(svdsigma$u)
            Sigmahalf <- svdsigma$v %*% diag(svdsigma$d^(1/2)) %*% t(svdsigma$u)
            Yk.tilde.star <- sqrt(Wk) %*% Y %*% Sigmanegahalf; Xk.star <- sqrt(Wk) %*% X
            res1 <- CH.11(Yk.tilde.star, Xk.star, lambda = lambda)
            para[[kk]][["B"]] <- res1$B %*% Sigmahalf
            penalty[kk] <- res1$penalty
          } else {
            Yk.star <- sqrt(Wk) %*% Y; Xk.star <- sqrt(Wk) %*% X
            res1 <- CH.11(Yk.star, Xk.star, lambda = lambda)
            para[[kk]][["B"]] <- res1$B
            penalty[kk] <- res1$penalty
          }
          # update for Sigma_k
          if(!commonvar){
            para[[kk]][["Sigma"]] <- temp.func2(kk, Y, X, para, K, p.ik)/colSums(p.ik)[kk]
          }
        }
        # update for Sigma if commonvar
        if(commonvar){
          eachSigma <- Reduce('+', lapply(1:K, temp.func2, Y, X, para, K, p.ik))/n
          for(kk in 1:K){
            para[[kk]][["Sigma"]] <- eachSigma
          }
        }
        penalty <- sum(penalty, na.rm = T)/2
        
      } else if(est == "ANNP"){
        penalty <- 0
        old.sum.nn.iter <- sum.nn.iter
        for(kk in 1:K){
          # update for pi_k
          para[[kk]][["pi"]] <- colSums(p.ik)[kk]/n
          # update for B_k: interval proximal gradient algorithm for B_k
          Wk <- diag(p.ik[, kk])
          if(!is.nan(para[[kk]][["Sigma"]][1])){
            svdsigma <- svd(para[[kk]][["Sigma"]])
            Sigmainv <- svdsigma$v %*% diag(svdsigma$d^(-1)) %*% t(svdsigma$u)
            tau <- sqrt(sum(sort(diag(Wk), decreasing = T)^2)) * sqrt(sum(svd(Sigmainv)$d^2))
            nn.iter <- 0; eerr <- -10;
            old.C <- C <- X %*% para[[kk]][["B"]]; old.t <- 1; t <- 1
            while(abs(eerr) > thres.int && nn.iter < maxiter.int){
              tilde.C <- C + (old.t - 1)/t * (C - old.C)
              H <- tilde.C - (1/tau) * Wk %*% (tilde.C - Y) %*% Sigmainv
              old.C <- C; old.tilde.C <- tilde.C
              prox.res <- CS.13(Y = H, X = X, lambda = lambda/tau,# gamma = gamma,
                                wi = svd(X %*% MASS::ginv(t(X) %*% X) %*% (t(X) %*% H))$d^(- gamma),
                                getpen = F)
              C <- X %*% prox.res$B
              old.t <- t
              t <- (1 + sqrt(1 + 4*t^2))/2
              old.para[[kk]][["B"]] <- para[[kk]][["B"]]
              para[[kk]][["B"]] <- prox.res$B
              r.hat[kk] <- prox.res$r.hat
              eerr <- sqrt(sum((old.C - C)^2))
              if(force){if(eerr > 0){break}}
              nn.iter <- nn.iter + 1
            }
            penalty <- sum(sapply(X = 1:K, FUN = function(kk){
              g(Bk = para[[kk]][["B"]], X = X, Y = Y, lambda = lambda, wi = wi)}))
            sum.nn.iter <- sum.nn.iter + nn.iter
          }
          # update for Sigma_k
          if(!commonvar){
            para[[kk]][["Sigma"]] <- temp.func2(kk, Y, X, para, K, p.ik)/colSums(p.ik)[kk]
          }
        }
        # update for Sigma if commonvar
        if(commonvar){
          eachSigma <- Reduce('+', lapply(1:K, temp.func2, Y, X, para, K, p.ik))/n
          for(kk in 1:K){
            para[[kk]][["Sigma"]] <- eachSigma
          }
        }
        if(est == "ANNP"){avg.nn.iter <- (sum.nn.iter-old.sum.nn.iter)/K}
      }
      
      loglh <- loglik(Y, X, para, K)
      new.loglik <- loglh - penalty
      err <- old.loglik - new.loglik
      if(force){if(err > 0){break}}
      old.loglik <- new.loglik
      n.iter <- n.iter + 1
      resmat <- rbind(resmat, c(n.iter, loglh, new.loglik, err, avg.nn.iter,
                                difftime(Sys.time(), start_time, units = "secs")))
      if(visible){
        print(c(n.iter, loglh, new.loglik, err, avg.nn.iter,
                difftime(Sys.time(), start_time, units = "secs")))
      }
    }
    
    # Deal with label switching issue when true parameters are known
    if(!is.null(para.true)){
      all.comb <- gtools::permutations(K, K)
      all.class.err <- rep(NA, nrow(all.comb))
      p.ik.true <- matrix(0, n, K)
      for(kk in 1:K){
        p.ik.true[ind.true == kk, kk] <- 1
      }
      for(cc in 1:nrow(all.comb)){
        all.class.err[cc] <- sum((p.ik.true - p.ik[, all.comb[cc, ]])^2)
      }
      ind.fin <- all.comb[which.min(all.class.err), ]
      para <- para[ind.fin]
      p.ik <- p.ik[, ind.fin]
    }
    
    # Compute the estimated ranks
    ind <- apply(p.ik, 1, which.max)
    p.ik.hard <- matrix(0, n, K)
    n.est <- rep(NA, K)
    est.rank <- est.rank0 <- rep(NA, K)
    for(kk in 1:K){
      p.ik.hard[ind == kk, kk] <- 1
      n.est[kk] <- length(which(ind == kk))
      est.rank0[kk] <- ifelse(is.na(Matrix::rankMatrix(para[[kk]][["B"]])[1]), 0,
                              Matrix::rankMatrix(para[[kk]][["B"]])[1])
    }
    if(est == "FR"){
      est.rank <- est.rank0
    } else if(est == "RP"){
      for(kk in 1:K){
        Wk <- diag(p.ik[, kk])
        if(!is.nan(para[[kk]][["Sigma"]][1])){
          svdsigma <- svd(para[[kk]][["Sigma"]])
          Sigmanegahalf <- svdsigma$v %*% diag(svdsigma$d^(-1/2)) %*% t(svdsigma$u)
          Sigmahalf <- svdsigma$v %*% diag(svdsigma$d^(1/2)) %*% t(svdsigma$u)
          Yk.tilde.star <- sqrt(Wk) %*% Y %*% Sigmanegahalf
          Xk.star <- sqrt(Wk) %*% X
          est.rank[kk] <- sum(svd(Xk.star %*% MASS::ginv(t(Xk.star) %*% Xk.star)
                                  %*% (t(Xk.star) %*% Yk.tilde.star))$d > lambda)
          #lambda^(1/(gamma + 1)))
        } else {
          Yk.star <- sqrt(Wk) %*% Y; Xk.star <- sqrt(Wk) %*% X
          est.rank[kk] <- sum(svd(Xk.star %*% MASS::ginv(t(Xk.star) %*% Xk.star)
                                  %*% (t(Xk.star) %*% Yk.star))$d > lambda)
          #lambda^(1/(gamma + 1)))
        }
      }
    } else if(est == "ANNP"){
      est.rank <- r.hat
    }
    for(kk in 1:K){
      est.rank[kk] <- ifelse(is.na(est.rank[kk]), 0, est.rank[kk])
    }
    
    # Compute BIC
    numc <- sum(n.est != 0)
    const <- ifelse(commonvar, 1, numc)
    npar <- (numc - 1) + const * q * (q + 1)/2 + sum(est.rank * (p + q - est.rank))
    npar0 <- (numc - 1) + const * q * (q + 1)/2 + sum(est.rank0 * (p + q - est.rank0))
    bic <- -2 * loglh + log(n) * npar
    bic0 <- -2 * loglh + log(n) * npar0
    
    # Compute estimation/prediction/classification errors when the true parameters are known
    if(!is.null(para.true)){
      
      if(length(para.true) != length(para)){ #10/3/2021 #add more possibilities?
        empty.ind <- which(lapply(para, '[[', 1) == 0)
        para <- para[- empty.ind]
        p.ik <- p.ik[, - empty.ind]
        p.ik.hard <- p.ik.hard[, - empty.ind]
        n.est <- n.est[- empty.ind]
        est.rank <- est.rank[- empty.ind]
        est.rank0 <- est.rank0[- empty.ind]
      }
      if(K != length(para.true)){
        K <- length(para.true)
      }
      
      p.ik.true <- matrix(0, n, K)
      est.err <- pred.err <- rep(NA, K)
      for(kk in 1:K){
        p.ik.true[ind.true == kk, kk] <- 1
        est.err[kk] <- 100 * sum(svd(para.true[[kk]][["B"]] - para[[kk]][["B"]])$d[1:min(p, q)]^2)/(p * q)
        pred.err[kk] <- 100 * sum(svd(X %*% para.true[[kk]][["B"]] - X %*% para[[kk]][["B"]])$d[1:min(n, q)]^2)/(n * q)
      }
      class.err <- c(sum((p.ik.true - p.ik)^2), sum((p.ik.true - p.ik.hard)^2))
    }
    
  }
  
  for(kk in 1:K){
    colnames(para[[kk]][["B"]]) <- colnames(Y)
    rownames(para[[kk]][["B"]]) <- colnames(X)
  }
  
  if(is.null(para.true)){
    output <- list(call = call, seed = seed, n.est = n.est, para = para,
                   est.rank = est.rank, npar = npar, n.iter = n.iter,
                   lambda = lambda, gamma = gamma, ind = ind, 
                   loglik = loglh, penalty = penalty, penloglik = new.loglik, bic = bic,
                   avg.nn.iter = sum.nn.iter/(n.iter * K), resmat = resmat)
  } else {
    output <- list(call = call, seed = seed, n.est = n.est, para = para,
                   est.rank = est.rank, npar = npar, n.iter = n.iter,
                   lambda = lambda, gamma = gamma, ind = ind, ind.true = ind.true,
                   loglik = loglh, penalty = penalty, penloglik = new.loglik, bic = bic,
                   avg.nn.iter = sum.nn.iter/(n.iter * K), resmat = resmat,
                   class.err = class.err, est.err = est.err, pred.err = pred.err)
  }
  
  class(output) <- "rrmix"
  return(output)
}

summary.rrmix <- function(object, ...){
  structure(object, class = "summary.rrmix")
}

plot.rrmix <- function(x, pch.L = 1, pch.F = 2, col.L = "red", col.F = "blue",
                       lty.L = 1, lty.F = 1, type = "b", ...){
  
  #graphics::par(mfrow = c(1, 2))
  plot(x$resmat[, "t"], x$resmat[, "F"], type = type, col = col.F, pch = pch.F,
       xlab = "Iteration", ylab = "(Penalized) log-likelihood",
       main = "(Penalized) log-likelihood",
       ylim = c(min(x$resmat[, c("L", "F")]), max(x$resmat[, c("L", "F")])))
  graphics::lines(x$resmat[, "t"], x$resmat[, "L"], type = type, col = col.L, pch = pch.L)
  graphics::legend("bottomright", legend = c("F", "L"), col = c(col.F, col.L),
                   pch = c(pch.F, pch.L), lty = c(lty.F, lty.L), bty = "n")
  #plot(x$resmat[, "t"], abs(x$resmat[, "old.F_new.F"]), type = type, col = col.F, pch = pch.F,
  #     xlab = "Iteration", ylab = "Absolute difference in penalized log-likelihood",
  #     main = "Absolute difference in F")
  
  cat("\nF: Penalized log-likelihood")
  cat("\nL: Log-likelihood")
  cat("\n\n")
  
}


### Estimation of coefficient matrix ###
# Rank penalization (Bunea et al., 2011)
CH.11 <- function(Y, X, lambda){
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X) 
  C.ols <- MASS::ginv(t(X)%*% X) %*% (t(X)%*% Y)  # OLS
  SVD <- svd(X %*% C.ols)
  U.hat <- SVD$u; V.hat <- SVD$v; D.hat <- diag(SVD$d)
  H <- diag(SVD$d[1:min(n, q)] * ifelse(SVD$d[1:min(n, q)] > lambda, 1, 0))
  out1 <- C.ols %*% V.hat %*% MASS::ginv(D.hat) %*% H %*% t(V.hat)
  out2 <- lambda^2 * sum(svd(X %*% MASS::ginv(t(X) %*% X) %*% (t(X) %*% Y))$d > lambda)
  return(list(B = out1, penalty = out2))
}
# Adaptive nuclear norm penalization (Chen et al., 2013)
CS.13 <- function(Y, X, lambda, wi, getpen = T){
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X) 
  C.ols <- MASS::ginv(t(X)%*% X) %*% (t(X)%*% Y)  # OLS
  SVD <- svd(X %*% C.ols)
  U.hat <- SVD$u; V.hat <- SVD$v; D.hat <- diag(SVD$d)
  S <- diag(ifelse(SVD$d[1:min(n, q)] - lambda * wi[1:min(n, q)] > 0,
                   SVD$d[1:min(n, q)] - lambda * wi[1:min(n, q)], 0))
  out1 <- C.ols %*% V.hat %*% MASS::ginv(D.hat) %*% S %*% t(V.hat)
  ## to further save time for "Chen2"
  if(getpen){
    di <- svd(X %*% out1)$d
    out2 <- 2 * lambda * sum(wi[1:min(n, q)] * di[1:min(n, q)])
  } else {
    out2 <- NA
  }
  out3 <- sum(SVD$d[1:min(n, q)] - lambda * wi[1:min(n, q)] > 0)
  return(list(B = out1, penalty = out2, r.hat = out3))
}

### Log-likelihood and penalized log-likelihood function
temp.func <- function(ii, Y, X, para, K, pd.ind = NULL){
  if(is.null(pd.ind)){
    for(kk in 1:K){
      para[[kk]][["Sigma"]] <- as.symmetric.matrix(para[[kk]][["Sigma"]])
    }
    pd.ind <- which(sapply(1:K, FUN = function(kk){is.valid.pd.mat(para[[kk]][["Sigma"]])}))
  }
  temp <- rep(0, K)
  for(kk in pd.ind){
    dec <- dmvnormC(t(as.matrix(Y[ii, ], nrow = 1)),
                        t(para[[kk]][["B"]]) %*% X[ii, ],
                        para[[kk]][["Sigma"]], FALSE)
    temp[kk] <- para[[kk]][["pi"]] * dec
  }
  return(temp)
}
# A function for calculating non-penalized log-likelihood when given parameter estimates
loglik <- function(Y, X, para, K){
  for(kk in 1:K){ 
    para[[kk]][["Sigma"]] <- as.symmetric.matrix(para[[kk]][["Sigma"]])
  }
  pd.ind <- which(sapply(1:K, FUN = function(kk){is.valid.pd.mat(para[[kk]][["Sigma"]])}))
  out <- sum(log(colSums(sapply(1:nrow(Y), temp.func, Y, X, para, K, pd.ind))))
  return(out)
}
# A function for checking whether a var-cov matrix is valid and positive definite
is.valid.pd.mat <- function(mat){
  if(!is.nan(mat)[1]){
    return(matrixcalc::is.positive.definite(mat))
  } else {
    return(FALSE)
  }
}
# A function for making a var-cov matrix symmetrical when it is not due to rounding errors
as.symmetric.matrix <- function(mat){
  if(!is.nan(mat)[1]){
    rr <- 10
    while(!matrixcalc::is.symmetric.matrix(mat) & rr > 3){
      mat <- round(mat, rr)
      rr <- rr - 1
    }
  }
  return(mat)
}
# 'f' function in the proximal gradient algorithm for ANNPM
f <- function(C, X, Y, W, SIGMA, lambda, gamma = 2){
  if(!is.nan(SIGMA[1])){
    svdsigma <- svd(SIGMA)
    Sigmanegahalf <- svdsigma$v %*% diag(svdsigma$d^(-1/2)) %*% t(svdsigma$u)
    Y.tilde.star <- sqrt(W) %*% Y %*% Sigmanegahalf
    out <- sum(svd(Y.tilde.star - sqrt(W) %*% C %*% Sigmanegahalf)$d^2)
  } else {
    Y.star <- sqrt(W) %*% Y
    out <- sum(svd(Y.star - sqrt(W) %*% C)$d^2)
  }
  return(out/2)
}
# 'g' function in the proximal gradient algorithm for ANNPM
g <- function(Bk, X, Y, lambda, wi){
  n <- nrow(Y); q <- ncol(Y); p <- ncol(X)
  di <- svd(X %*% Bk)$d
  out <- lambda * sum(wi[1:min(n, q)] * di[1:min(n, q)])
  return(out)
}
# A function for calculating the update of Sigma
temp.func2 <- function(kk, Y, X, para, K, p.ik){
  t(Y - X %*% para[[kk]][["B"]]) %*% diag(p.ik[, kk]) %*% (Y - X %*% para[[kk]][["B"]])
}