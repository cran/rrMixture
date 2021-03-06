##' Summarize rrmix Objects
##'
##' S3 methods summarizing objects generated by \code{rrmix} and \code{tune.rrmix}.
##'
##' @name summary
##'
##' @param object Object generated from \code{rrmix} or \code{tune.rrmix}.
##' @param metric performance metric to use for finding best `rrmix' model.
##'   `soft.class.err', `hard.class.err', `est.err', and `pred.err' can only be
##'   used when true parameter values are known.
##' @param ... Other arguments for future usage.
NULL


##' @rdname summary
##' @export
summary.rrmix <- function(object, ...)
{
  class(object) <- "summary.rrmix"
  object
}


##' @rdname summary
##' @export
summary.tune.rrmix <- function(object, metric, ...)
{
  class(object) <- "summary.tune.rrmix"
  object
}