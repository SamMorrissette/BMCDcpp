// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "commonvars.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}


double CalcSSR(Rcpp::NumericMatrix mat1, Rcpp::NumericMatrix mat2) {
  double SSR = 0;
  for(int i=0; i<n; i++)
    SSR += sum((mat1.row(i)- mat2.row(i))*(mat1.row(i)- mat2.row(i)));
  SSR = SSR/2.0;
  return SSR;
}

// Credit for this function goes to Matthew Denny: https://www.mjdenny.com/blog.html

// [[Rcpp::export]]
arma::mat rdirichlet_cpp(int num_samples,
                         arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each row will be a draw from a Dirichlet
  arma::mat distribution = arma::zeros(num_samples, distribution_size);

  for (int i = 0; i < num_samples; ++i) {
    double sum_term = 0;
    // loop through the distribution and draw Gamma variables
    for (int j = 0; j < distribution_size; ++j) {
      double cur = R::rgamma(alpha_m[j],1.0);
      distribution(i,j) = cur;
      sum_term += cur;
    }
    // now normalize
    for (int j = 0; j < distribution_size; ++j) {
      distribution(i,j) = distribution(i,j)/sum_term;
    }
  }
  return(distribution);
}

// [[Rcpp::export]]
arma::rowvec Arma_colSums(const arma::mat& x) {
  return arma::sum(x, 0);
}

// arma::mat ScaleMatrix(arma::mat X) {
//   arma::mat mat_out = X;
//   arma::vec colMeans(dim);
//   arma::vec colSD(dim);
//
//   for (int i = 0; i < n; i++) {
//     for (int d = 0; d < dim; d++) {
//       colMeans(d) = mean(X.col(d));
//       colSD(d) = sqrt(var(X.col(d)));
//       mat_out(i, d) = (X(i, d) - colMeans(d))/colSD(d);
//     }
//   }
//   return mat_out;
// }

