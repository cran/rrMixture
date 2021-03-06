#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <math.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat const &x,
                      arma::vec const &center,
                      arma::mat const &cov) {
  arma::mat x_cen = x.t();
  x_cen.each_col() -= center;
  arma::solve(x_cen, arma::trimatl(chol(cov).t()), x_cen);
  x_cen.for_each( [](arma::mat::elem_type& val) {val = val * val;});
  return arma::sum(x_cen, 0).t();
}

// [[Rcpp::export]]
arma::vec dmvnormC(arma::mat const &x, 
                       arma::vec const &mean, 
                       arma::mat const &sigma, 
                       bool const logd = false) { 
  arma::vec const distval = Mahalanobis(x,  mean, sigma);
  double const logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec const logretval = 
    -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
    
    if (logd)
      return logretval;
    return exp(logretval);
}
