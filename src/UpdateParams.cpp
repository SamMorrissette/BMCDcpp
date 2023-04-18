// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
#include "commonvars.h"
#include "helpers.h"

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat CalculateClassProbs(arma::rowvec p,
                              arma::mat mean,
                              arma::cube cov) {
  IntegerVector ints = seq_len(num_comps);
  arma::rowvec probs(num_comps);
  double numerator, denom;
  arma::rowvec out_vec(n);
  arma::mat prob_mat(n, num_comps);
  arma::vec result(1);

  for (int i = 0; i < n; i++) {
    denom = 0;
    for (int k = 0; k < num_comps; k++) {
      result = dmvnorm(X_mat.row(i), mean.col(k), cov.slice(k));
      denom = denom + p(k) * result(0);
    }

    for (int k = 0; k < num_comps; k++) {
      result = dmvnorm(X_mat.row(i), mean.col(k), cov.slice(k));
      numerator = p(k) * result(0);
      probs(k) = numerator/denom;
    }
    prob_mat.row(i) = probs;
  }
  return prob_mat;
}

arma::rowvec UpdateClasses(arma::mat prob_mat) {
  IntegerVector ints = seq_len(num_comps);
  NumericVector probs(num_comps);
  arma::rowvec out_vec(n);
  for (int i = 0; i < n; i++) {
    probs = as<NumericVector>(wrap(prob_mat.row(i)));
    out_vec(i) = as<int>(wrap(sample(ints, 1, true, probs)));
  }
  return out_vec;
}

arma::rowvec UpdateWeights(arma::rowvec z) {
  arma::rowvec out_vec(num_comps);
  arma::vec n_vec(num_comps);
  arma::uvec pos;
  arma::mat dir_sample;

  for (int k = 0; k < num_comps; k++) {
    pos = find(z == k+1);
    n_vec(k) = pos.n_elem + 1;
  }

  dir_sample = rdirichlet_cpp(1, n_vec);
  out_vec = dir_sample.row(0);
  return out_vec;
}

Params UpdateTheta(arma::rowvec z, int modelIndex) {
  Params params_out;
  switch (modelIndex) {
  case 1: // Unequal Unrestricted
  {
    params_out = DrawUnequalUnrestricted(z);
    break;
  }
  case 2: // Unequal Spherical
  {
    params_out = DrawUnequalSpherical(z);
    break;
  }
  case 3: // Unequal Diagonal;
    {
    params_out = DrawUnequalDiagonal(z);
    break;
    }
  case 4: // Equal Unrestricted;
    {
      params_out = DrawEqualUnrestricted(z);
      break;
    }
  case 5: // Equal Spherical;
    {
      params_out = DrawEqualSpherical(z);
      break;
    }
  case 6: // Equal Diagonal;
    {
      params_out = DrawEqualDiagonal(z);
      break;
    }
  }
  return params_out;
}
