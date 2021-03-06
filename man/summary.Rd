% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary}
\alias{summary}
\alias{summary.rrmix}
\alias{summary.tune.rrmix}
\title{Summarize rrmix Objects}
\usage{
\method{summary}{rrmix}(object, ...)

\method{summary}{tune.rrmix}(
  object,
  metric = c("bic", "soft.class.err", "hard.class.err", "est.err", "pred.err"),
  ...
)
}
\arguments{
\item{object}{Object generated from \code{rrmix} or \code{tune.rrmix}.}

\item{...}{Other arguments for future usage.}

\item{metric}{performance metric to use for finding best `rrmix' model.
`soft.class.err', `hard.class.err', `est.err', and `pred.err' can only be
used when true parameter values are known.}
}
\description{
S3 methods summarizing objects generated by \code{rrmix} and \code{tune.rrmix}.
}
