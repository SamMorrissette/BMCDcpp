// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppDist.h>
#include "commonvars.h"
#include "helpers.h"

using namespace Rcpp;

double CalcBIC(arma::cube X, arma::mat p, arma::cube means, arma::field<arma::cube> covs,
               const int G,
               int modelIndex,
               int burn, int iters) {
  double BIC;
  int num_params;
  int df = (G - 1) + (G*dim);
  double lik_iter = 0;
  double log_likelihood = 0;
  arma::mat mean_X(n, dim, arma::fill::zeros);
  arma::rowvec mean_p(G, arma::fill::zeros);
  arma::mat mean_mu(dim, G, arma::fill::zeros);
  arma::cube mean_cov(dim, dim, G);
  arma::vec normd_result(1);
  for (int t = burn; t < iters; t++) {
    mean_X = mean_X + X.slice(t);
    mean_p = mean_p + p.row(t);
    mean_mu = mean_mu + means.slice(t);
    for (int k = 0; k < num_comps; k++) {
      mean_cov.slice(k) = mean_cov.slice(k) + (covs(t).slice(k) / (iters - burn));
    }
  }
  mean_X = mean_X / (iters - burn);
  mean_p = mean_p / (iters - burn);
  mean_mu = mean_mu / (iters - burn);
  switch (modelIndex) {
  case 1: // Unequal Unrestricted
  {
    num_params = df + G*(dim*(dim+1) / 2);
    break;
  }
  case 2: // Unequal Spherical
  {
    num_params = df + G;
    break;
  }
  case 3: // Unequal Diagonal;
    {
    num_params = df + (G * dim);
    break;
    }
  case 4: // Equal Unrestricted;
    {
      num_params = df + (dim*(dim+1) / 2);
      break;
    }
  case 5: // Equal Spherical;
    {
      num_params = df + 1;
      break;
    }
  case 6: // Equal Diagonal;
    {
      num_params = df + dim;
      break;
    }
  }

  for (int i = 0; i < n; i++) {
    lik_iter = 0;
    for (int j = 0; j < G; j++) {
      normd_result = mean_p(j) * dmvnorm(mean_X.row(i), mean_mu.col(j), mean_cov.slice(j));
      lik_iter = lik_iter + normd_result(0);
    }
    log_likelihood = log_likelihood + log(lik_iter);
  }
  BIC = (log(n) * num_params) - 2 * (log_likelihood);
  return BIC;
}
