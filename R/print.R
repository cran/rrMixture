### objects from rrmix.R ========================================================
##' @export
print.rrmix <- function(x, ...)
{
  cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n")
  cat("\nMethod: ", switch(x$call$est, "FR" = "Full ranked mixture",
                           "RP" = "Rank penalized mixture",
                           "ANNP" = "Adaptive nuclear norm penalized mixture"), "\n")
  if(is.null(x$call$para0) && is.null(x$call$ind0)){
    #init <- switch(x$call$init, "Random" = "Random clustering",
    #               "Kmeans" = "K-means clustering")
    init <- "K-means and random clustering"
  } else {
    init <- "Provided"
  }
  cat("\nInitialization: ", init, "\n")
  
  if (x$call$est == "RP" || x$call$est == "ANNP"){
    cat("\n")
    cat("Tuning Parameters:\n")
    cat(" lambda: ", x$lambda, "\n")
    if (x$call$est == "ANNP"){
      cat("  gamma: ", x$gamma, "\n")
    }
  }
  
  cat("\n")
  cat("Fitted Model:\n")
  cat("     Number of components: ", sum(x$n.est != 0), "\n")
  cat("           Log-likelihood: ", round(x$loglik, 4), "\n")
  cat(" Penalized log-likelihood: ", round(x$penloglik, 4), "\n")
  cat("                      BIC: ", round(x$bic, 4), "\n")
  
  cat("\n")
  
  invisible(x)
}

##' @export
print.summary.rrmix <- function(x, ...)
{
  print.rrmix(x)
  
  cat("\nNumber of Iterations: ", x$n.iter, "\n")
  
  cat("\n")
  if (!is.null(x$call$ind.true) || !is.null(x$call$para.true)){
    cat("Errors:\n")
    cat(" Soft classification error: ", round(x$class.err[1], 4), "\n")
    cat(" Hard classification error: ", round(x$class.err[2], 4), "\n")
    cat("          Estimation error: ", round(mean(x$est.err), 4), "\n")
    cat("          Prediction error: ", round(mean(x$pred.err), 4), "\n")
  }
  
  cat("\n")
  if (!is.null(x$call$ind.true) || !is.null(x$call$para.true)){
    cat("Confusion Matrix for Mixture Membership (Estimated vs True):\n")
    print(table(x$ind, x$ind.true))
  }
  
  cat("\n\n")
  
  invisible(x)
}


### for objects from tune.R =====================================================
##' @export
print.tune.rrmix <- function(x, ...){
  
  ## print results
  #cat("\n")
  #print(tab)
  #cat("\n")
  cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep = "\n")
  cat("Parameter tuning of ", sQuote("rrmix"), ":\n\n", sep = "")
  cat("- Method:", switch(x$est, "FR" = "Full ranked mixture",
                          "RP" = "Rank penalized mixture",
                          "ANNP" = "Adaptive nuclear norm penalized mixture"),"\n\n")
  
  if(is.null(x$call$para0) && is.null(x$call$ind0)){
    #init <- switch(x$call$init, "Random" = "Random clustering",
    #               "Kmeans" = "K-means clustering")
    init <- "K-means and random clustering"
  } else {
    init <- "Provided"
  }
  cat("- Initialization: ", init, "\n")
  cat("\n")
  
  invisible(x)
  
}

##' @export
print.summary.tune.rrmix <- function(x, ...) {
  
  best.ind <- as.numeric(which.min(x$details[, x$metric]))
  if(x$est == "FR" || x$est == "RP"){
    parameters <- c("K", "lambda")
  } else if (x$est == "ANNP"){
    parameters <- c("K", "lambda", "gamma")
  }
  
  ## print
  print.tune.rrmix(x)
  cat("- Performance metric:", x$metric, "\n")
  cat("\n- Best performance:", x$details[best.ind, x$metric], "\n")
  cat("\n- Best parameters:\n")
  tmp <- data.frame(matrix(x$details[best.ind, parameters], nrow = 1))
  colnames(tmp) <- parameters
  rownames(tmp) <- ""
  print(tmp)
  cat("\n- Detailed performance results:\n")
  print(x$details)
  cat("\n")
  
  invisible(x)
  
}
