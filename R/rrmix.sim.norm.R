#' Simulation Data Generator
#'
#' `rrmix.sim.norm' is used to create synthetic data from the multivariate
#' normal distribution, which is used in a numerical study of Kang et. al.
#' (2022+).
#'
#' @param K number of mixture components.
#' @param n number of observations.
#' @param p number of predictors including an intercept.
#' @param q number of responses.
#' @param rho correlation between predictors used to make a design matrix.
#' @param b signal strength which controls the magnitude of coefficient matrices.
#' @param shift mean shift which measures how separate the mixture components are.
#' @param r.star vector of length K, specifying the true ranks of K coefficient
#'   matrices.
#' @param sigma vector of length K, specifying the noise strength of K
#'   multivariate normal distributions.
#' @param pr vector of length K, specifying the multinomial probabilities for
#'   the K mixture components.
#' @param seed seed number for the reproducibility of results. Default is `NULL'.
#'
#' @return
#'
#'   \item{X}{n by p design matrix.}
#'   \item{Y}{n by q response matrix.}
#'   \item{E}{p by q error matrix.}
#'   \item{ind.true}{vector of length n, specifying the true mixture membership
#'   for n observations.}
#'   \item{para.true}{array of length K. It consists of K lists, each of which
#'   contains a coefficient matrix and its true rank.}
#'
#' @author
#'   Suyeon Kang, University of California, Riverside, \email{skang062@@ucr.edu};
#'   Weixin Yao, University of California, Riverside, \email{weixin.yao@@ucr.edu};
#'   Kun Chen, University of Connecticut, \email{kun.chen@@uconn.edu}.
#'
#' @references
#'   Kang, S., Chen, K., and Yao, W. (2022+). "Reduced rank estimation in mixtures
#'   of multivariate linear regression".
#'   
#' @importFrom stats rnorm rmultinom
#' @export
#' @examples
#' #-----------------------------------------------------------#
#' # Simulation 1: Two Components Case
#' #-----------------------------------------------------------#
#' K2mod <- rrmix.sim.norm(K = 2, n = 100, p = 5, q = 5, rho = .5,
#'          b = 1, shift = 1, r.star = c(1, 3), sigma = c(1, 1),
#'          pr = c(.5, .5), seed = 1215)
#' 
#' #-----------------------------------------------------------#
#' # Simulation 2: Four Components Case
#' #-----------------------------------------------------------#
#' K4mod <- rrmix.sim.norm(K = 4, n = 600, p = 15, q = 15,
#'          rho = .5, b = 1, shift = 1, r.star = c(1, 1, 3, 3),
#'          sigma = c(1, 1, 1, 1), pr = c(.25, .25, .25, .25),
#'          seed = 1215)
rrmix.sim.norm <- function(K = 2,
                           n = 100,
                           p = 5,
                           q = 5,
                           rho = .5,
                           b = 1,
                           shift = 1,
                           r.star = NULL,
                           sigma = NULL,
                           pr = NULL,
                           seed = NULL){
  
  # Warning for inappropriate settings
  if(rho <= 0 | rho >= 1){
    stop("'rho' must be between 0 and 1")
  }
  
  if(!is.null(r.star) && length(r.star) != K){
    stop("length of 'r.star' must be equal to 'K'")
  } else if (!is.null(r.star) && length(r.star) == K){
    if(!all(r.star <= rep(min(p, q), K))){
      stop("number in 'r.star' cannot be greater than the minimum of 'p' and 'q'")
    }
  } else if (is.null(r.star)){
    r.star <- rep(min(p, q), K)
  }
  
  if(!is.null(sigma) && length(sigma) != K){
    stop("length of 'sigma' must be equal to 'K'")
  } else if (is.null(sigma)){
    sigma <- rep(1, K)
  }
  
  if(!is.null(pr) && length(pr) != K){
    stop("length of 'pr' must be equal to 'K'")
  } else if (!is.null(pr) && sum(pr) != 1){
    pr <- pr/sum(pr)
  } else if (is.null(pr)){
    pr <- rep(1/K, K)
  } 

  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #-----------------------------------------------------#
  # Data Generation
  #-----------------------------------------------------#
  # Design matrix
  Gamma <- rho^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, rep(0, p), Gamma)
  X <- cbind(intercept = 1, X); p <- p + 1
  # Number of units from each component
  ns <- rmultinom(1, n, pr)
  # True coefficient matrices, Error matrix, Response matrix
  B <- array(NA, dim = c(p, q, K))
  E <- Y <- NULL
  for(k in 1:K){
    B[, , k] <- b * matrix(rnorm(p * r.star[k], mean = (k - 1) * shift), p, r.star[k]) %*%
      t(matrix(rnorm(q * r.star[k], mean = (k - 1) * shift), q, r.star[k]))
    E <- rbind(E, matrix(rnorm(ns[k] * q, 0, sigma[k]), ns[k], q))
    Y <- rbind(Y, X[ifelse(k == 1, 1, sum(ns[1:(k - 1)]) + 1):sum(ns[1:k]), ] %*% B[, , k] +
                 E[ifelse(k == 1, 1, sum(ns[1:(k - 1)]) + 1):sum(ns[1:k]), ])
  }
  
  #-----------------------------------------------------#
  # True Parameters
  #-----------------------------------------------------#
  # True parameters
  ind.true <- rep(c(1:K), ns)
  para.true <- array(list(), K)
  for(k in 1:K){
    para.true[[k]][["r"]] <- r.star[k]
    para.true[[k]][["B"]] <- B[, , k]
  }
  
  out <- list(X = X, Y = Y, E = E, ind.true = ind.true, para.true = para.true)
  return(out)
  
}
