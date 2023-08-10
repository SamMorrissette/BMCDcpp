// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <progress.hpp>
#include "helpers.h"
#include "commonvars.h"


arma::mat X_mat;
int n, dim, num_comps, num_iters;
double sigmasq_prior_shape;
double sigmasq_prior_scale;
const double prior_shape = 4, prior_rate = 2; //maybe not needed?
arma::vec mu0;
arma::mat S;
double nu0, kappa;
double prior_IG_alpha;
double prior_IG_beta;
Rcpp::Function RDist("dist");
Rcpp::Function as_matrix("as.matrix");

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::export]]
List BMCD_MCMC(arma::mat obs_dist,
             arma::mat init_X, const int G, double init_sigmasq,
             const int burn, const int iters,
             int modelIndex) {
  X_mat = init_X;
  n = init_X.n_rows;
  dim = init_X.n_cols;
  num_comps = G;
  num_iters = iters;

  Params newParams;
  double s_res = CalcSSR(distRcpp(wrap(init_X)), wrap(obs_dist));

  // Priors for NIW (UU model)
  mu0 = arma::vec(dim, arma::fill::zeros);
  nu0 = dim+2;
  S = arma::mat(dim, dim, arma::fill::eye);
  kappa = 3;

  // Priors for IG (US model)
  arma::mat cov_X_init = cov(X_mat);
  prior_IG_alpha = (dim+2) / 2;
  prior_IG_beta = (dim+2) / 2; //(sum(cov_X_init.diag()) / dim) / 2; //

  // Priors for sigma squared (Inverse-Gamma prior)
  sigmasq_prior_shape = 5;
  sigmasq_prior_scale = (sigmasq_prior_shape-1) * (s_res / (n*(n-1)/2));

  // Declare Chains
  arma::cube X(n, dim, iters);
  arma::vec sigmasq(iters);
  arma::mat p(iters, num_comps);
  arma::mat z(iters, n);
  arma::cube class_probs(n, num_comps, iters);
  arma::cube means(dim, num_comps, iters);
  arma::field<arma::cube> covs(iters);
  covs.fill(arma::zeros<arma::cube>(dim, dim, num_comps));

  // Label-switching variables
  Theta_calcs Theta;
  arma::mat cost_mat(num_comps, num_comps);
  arma::uvec perms;
  arma::mat perm_mat(iters, num_comps, arma::fill::zeros);
  // arma::cube old_probs(n, num_comps, iters, arma::fill::zeros);

  // Initialize the chains
  sigmasq(0) = init_sigmasq;
  InitChains(X, p, z, means, covs);

  double BIC;
  Progress prog(iters, true);


  for (int t = 1; t < iters; t++) {
    if (Progress::check_abort() )
      return -1.0;
    prog.increment();

    // Step 1 (Update X)
    X.slice(t) = UpdateX(obs_dist, X.slice(t-1), sigmasq(t-1), z.row(t-1), means.slice(t-1), covs(t-1));

    // Step 1.5 (Transform X)
    X.slice(t) = TransformX(X.slice(t), init_X);

    X_mat = X.slice(t); // Update global X matrix

    // Step 2 (Update sigma_sq measurement error)
    sigmasq(t) = UpdateSigmasq(obs_dist, X.slice(t), s_res, sigmasq(t-1));

    // Step 3 (drawing classes for each data point)
    class_probs.slice(t) = CalculateClassProbs(p.row(t-1), means.slice(t-1), covs(t-1));

    // Step 5 (update allocations)
    z.row(t) = UpdateClasses(class_probs.slice(t));

    // Step 6 (update component weights)
    p.row(t) = UpdateWeights(z.row(t));

    // Step 7 (Update model parameters)
    newParams = UpdateTheta(z.row(t), modelIndex);

    means.slice(t) = newParams.mean;
    covs(t) = newParams.covariance;

    //Step 8 (Label switching)
    if (t == 100) {
      Theta = InitializeCenters(p, means, covs);
    } else if (t > 100) {
      cost_mat = CostMatrix(Theta, p.row(t), means.slice(t), covs(t));
      perms = Hungarian(cost_mat);
      perm_mat.row(t) = arma::conv_to<arma::rowvec>::from(perms);
      // Might want to add a condition in here...
      ChangeLabels(t, perms, p, means, covs);
      // What about z and class_probs?
      Theta = UpdateCenters(Theta, t, p.row(t), means.slice(t), covs(t));
    }
  }

  //Step 9 (Calculate BIC)
  BIC = CalcBIC(X, p, means, covs, G, modelIndex, burn, iters);
  return List::create(
    Named("X") = X,
    Named("sigmasq") = sigmasq,
    Named("means") = means,
    Named("covs") = covs,
    Named("z") = z,
    Named("class_probs") = class_probs,
    Named("p") = p,
    Named("BIC") = BIC);
}
