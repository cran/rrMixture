#' Initialization of Parameter Estimates
#'
#' `initialize.para' is used to initialize parameter estimates.
#' 
#' @usage
#' initialize.para(K, X, Y, ind0 = NULL,
#'                 seed = NULL, km.nstart = 20, kmscale = FALSE, n.init = 100,
#'                 commonvar = FALSE)
#'
#' @param K number of mixture components.
#' @param X n by p design matrix where n is the number of observations and
#'   p is the number of predictors.
#' @param Y n by q response matrix where n is the number of observations and
#'   q is the number of responses.
#' @param ind0 vector of length n, specifying the initial assignment of the
#'   mixture membership of n observations when there is prior information on
#'   the membership. If `NULL', K-means clustering technique is used to assign
#'   the membership for n observations. Default is `NULL'.
#' @param seed seed number for the reproducibility of results. Default is `NULL'.
#' @param km.nstart number of random sets considered to perform K-means
#'   clustering. Only used for K-means clustering. Default is 20.
#' @param kmscale logical value, indicating whether Y is scaled prior to K-means
#'   clustering. Only used for K-means clustering. Default is `FALSE'.
#' @param n.init number of initializations to try. Two methods for initial
#'   clustering are used: K-means and random clustering.
#' @param commonvar logical value, indicating the homogeneity assumption of
#'   variance-covariance matrices across K mixture components. Default is `FALSE'.
#'
#' @return
#'
#'   \item{para}{array of length K. It consists of K lists, each of which contains
#'   initial estimates of membership probability, coefficient matrix, and variance-
#'   covariance matrix.}
#'
#' @author
#'   Suyeon Kang, University of California, Riverside, \email{skang062@@ucr.edu};
#'   Weixin Yao, University of California, Riverside, \email{weixin.yao@@ucr.edu};
#'   Kun Chen, University of Connecticut, \email{kun.chen@@uconn.edu}.
#'
#' @seealso \code{\link{rrmix.sim.norm}}
#'
#' @references
#'   Kang, S., Chen, K., and Yao, W. (2022+). "Reduced rank estimation in mixtures
#'   of multivariate linear regression".
#'   
#' @importFrom stats kmeans
#' @export
#' @examples
#' #-----------------------------------------------------------#
#' # Simulation 1: Two Components Case
#' #-----------------------------------------------------------#
#' K2mod <- rrmix.sim.norm(K = 2, n = 100, p = 5, q = 5, rho = .5,
#'          b = 1, shift = 1, r.star = c(1, 3), sigma = c(1, 1),
#'          pr = c(.5, .5), seed = 1215)
#' K2ini <- initialize.para(K = 2, X = K2mod$X, Y = K2mod$Y,
#'          seed = 100)
#' 
#' #-----------------------------------------------------------#
#' # Simulation 2: Four Components Case
#' #-----------------------------------------------------------#
#' \donttest{
#' K4mod <- rrmix.sim.norm(K = 4, n = 600, p = 15, q = 15,
#'          rho = .5, b = 1, shift = 1, r.star = c(1, 1, 3, 3),
#'          sigma = c(1, 1, 1, 1), pr = c(.25, .25, .25, .25),
#'          seed = 1215)
#' K4ini <- initialize.para(K = 4, X = K4mod$X, Y = K4mod$Y,
#'          seed = 100)}
initialize.para <- function(K,
                            X,
                            Y,
                            ind0 = NULL,
                            seed = NULL,
                            km.nstart = 20,
                            kmscale = FALSE,
                            n.init = 100,
                            commonvar = FALSE)
{
  # Match arguments
  X <- as.matrix(X); Y <- as.matrix(Y)
  n <- nrow(Y)
  
  # Initial mixture membership
  if(is.null(ind0)){
    
    init.grp <- matrix(NA, n, n.init)
    init.size <- matrix(NA, K, n.init)
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    if(kmscale){
      km <- kmeans(scale(Y), K, nstart = km.nstart)
    } else {
      km <- kmeans(Y, K, nstart = km.nstart)
    }
    init.grp[, 1] <- km$cluster
    init.size[, 1] <- sort(km$size)
    for(i in 2:n.init){
      init.grp[, i] <- sample(1:K, n, replace = TRUE)
    }
 
  } else {
    init.grp <- matrix(ind0, n, 1)
  }
  
  # Make an array for parameter estimates
  init.llh <- rep(NA, ifelse(is.null(ind0), n.init, 1))
  for(i in 1:n.init){
    
    para <- array(list(), K)
    for(kk in 1:K){
      para[[kk]] <- as.list(para[[kk]])
    }
    if(commonvar){
      H <- X %*% MASS::ginv(t(X) %*% X) %*% t(X)
    }
    for(kk in 1:K){
      Xk <- X[init.grp[, i] == kk, ]; Yk <- Y[init.grp[, i] == kk, ]
      Hk <- Xk %*% MASS::ginv(t(Xk) %*% Xk) %*% t(Xk)
      nk <- length(which(init.grp[, i] == kk))
      para[[kk]][["pi"]] <- nk/n
      para[[kk]][["B"]] <- MASS::ginv(t(Xk) %*% Xk) %*% (t(Xk) %*% Yk)
      if(commonvar){
        para[[kk]][["Sigma"]] <- (1/n) * t(Y) %*% (diag(n) - H) %*% Y
        #J <- matrix(1, n, 1) %*% t(matrix(1, n, 1))
        #para[[kk]][["Sigma"]] <- (1/n) * t(Y) %*% (diag(n) - J) %*% Y
      } else {
        para[[kk]][["Sigma"]] <- (1/nk) * t(Yk) %*% (diag(nk) - Hk) %*% Yk
        #Jk <- matrix(1, nk, 1) %*% t(matrix(1, nk, 1))
        #para[[kk]][["Sigma"]] <- (1/nk) * t(Yk) %*% (diag(nk) - Jk) %*% Yk
      }
      para[[kk]][["Sigma"]] <- as.symmetric.matrix(para[[kk]][["Sigma"]])
    }
    if(K == 1){
      init.llh[i] <- sum(log(sapply(1:nrow(Y), temp.func, Y, X, para, K, pd.ind = 1)))
    } else {
      init.llh[i] <- loglik(Y, X, para, K)
    }
    
    # finding parameter estimates which provide the smallest log-likelihood
    if(is.null(ind0)){
      if(i == 1){
        para.min <- para
      } else {
        if(init.llh[i] > min(init.llh[1:(i - 1)], na.rm = TRUE)){
          para.min <- para
        }
      }
    } else {
      para.min <- para
    }
    
  }
  
  para <- para.min
  return(para)
  
}
