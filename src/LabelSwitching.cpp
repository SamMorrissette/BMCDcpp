// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "commonvars.h"
#include "RcppHungarian.h"
#include "helpers.h"

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppHungarian)]]

Theta_calcs InitializeCenters(arma::mat p, arma::cube means, arma::field<arma::cube> covs) {
  int switch_iters = 100;
  Theta_calcs Theta;
  arma::rowvec Theta_bar_eps(num_comps), var_eps(num_comps);
  arma::mat Theta_bar_mu(dim, num_comps), var_mu(dim, num_comps);
  arma::cube Theta_bar_cov(dim, dim, num_comps), var_cov(dim, dim, num_comps);

  // Find theta bar for each parameter in the model

  for (int t = 0; t < switch_iters; t++) {
    // Epsilons
    Theta_bar_eps = Theta_bar_eps + p.row(t);

    // Means
    Theta_bar_mu = Theta_bar_mu + means.slice(t);

    //Covariances
    for (int k = 0; k < num_comps; k++) {
      Theta_bar_cov.slice(k) = Theta_bar_cov.slice(k) + (covs(t).slice(k) / switch_iters);
    }
  }

  Theta_bar_eps = Theta_bar_eps / switch_iters;
  Theta_bar_mu = Theta_bar_mu / switch_iters;

  // Find variance (s) for each parameter

  for (int t = 0; t < switch_iters; t++) {
    // Variances of epsilon parameters (vector)
    var_eps = var_eps + arma::square((p.row(t) - Theta_bar_eps));

    // Variances of mean parameters (matrix)
    var_mu = var_mu + arma::square((means.slice(t) - Theta_bar_mu));

    // Variances of covariances parameters (cube)
    for (int k = 0; k < num_comps; k++) {
      var_cov.slice(k) = var_cov.slice(k) + (arma::square((covs(t).slice(k) - Theta_bar_cov.slice(k))) / switch_iters);
    }
  }

  var_eps = var_eps / switch_iters;
  var_mu = var_mu / switch_iters;

  Theta.eps_bar = Theta_bar_eps;
  Theta.mu_bar = Theta_bar_mu;
  Theta.cov_bar = Theta_bar_cov;
  Theta.eps_var = var_eps;
  Theta.mu_var = var_mu;
  Theta.cov_var = var_cov;
  return Theta;
}


arma::mat CostMatrix(Theta_calcs Theta,
                     arma::rowvec p,
                     arma::mat mean,
                     arma::cube cov) {

  arma::mat out_mat (num_comps, num_comps);
  arma::mat eps_mat (num_comps, num_comps);
  arma::mat mu_mat (num_comps, num_comps), mu_dist(num_comps, num_comps), mu_norm_dist(num_comps, num_comps);
  arma::mat cov_mat (num_comps, num_comps);
  arma::cube cov_dist(dim, dim, num_comps), cov_norm_dist(dim, dim, num_comps);

  for (int k = 0; k < num_comps; k++) {
    // Epsilon
    eps_mat.row(k) = arma::square(p(k) - Theta.eps_bar) / Theta.eps_var(k);

    // Mean
    mu_dist = arma::square(mean.col(k) - Theta.mu_bar.each_col());
    mu_norm_dist = mu_dist.each_col() / Theta.mu_var.col(k);
    mu_mat.row(k) = Arma_colSums(mu_norm_dist);

    // Covariance
    cov_dist = cov.slice(k)  - Theta.cov_bar.each_slice();
    for (int j = 0; j < num_comps; j++) {
      cov_dist.slice(j) = arma::square(cov_dist.slice(j));
      cov_norm_dist.slice(j) = cov_dist.slice(j) / Theta.cov_var.slice(k);
      cov_mat(k, j) = accu(cov_norm_dist.slice(j));
    }
  }
  out_mat = eps_mat + mu_mat + cov_mat;
  return out_mat;
}


// This function was taken from the Rcpp hungarian package!
IntegerVector HungarianSolver(NumericMatrix costMatrix) {
  int nr = costMatrix.nrow();
  int nc = costMatrix.ncol();

  vector<double> c(nc);
  vector<vector<double>> cm(nr, c);
  for (int i=0; i < nr; i++){
    for (int j=0; j < nc; j++){
      c[j] = costMatrix(i,j);
    }
    cm[i] = c;
  }

  HungarianAlgorithm HungAlgo;
  vector<int> assignment;
  double cost = HungAlgo.Solve(cm, assignment);
  IntegerVector assign(nr);
  for (int i=0; i < nr; i++){
    assign(i) = assignment[i]+1;
  }
  return assign;
}

arma::uvec Hungarian(arma::mat cost_mat) {
  NumericMatrix hung_input;
  IntegerVector hung_output;
  arma::uvec perm_vec(num_comps);

  hung_input = as<NumericMatrix>(wrap(cost_mat));
  hung_output = HungarianSolver(hung_input);
  perm_vec = as<arma::uvec>(wrap(hung_output));
  return perm_vec;
}

Theta_calcs UpdateCenters(Theta_calcs Theta, int t,
                        arma::rowvec p,
                        arma::mat mean,
                        arma::cube cov) {
  Theta_calcs newTheta;
  double m = 100;
  double r = t - m;
  const double OLD_WEIGHT = (m + r - 1) / (m + r);
  const double NEW_WEIGHT = 1 / (m + r);
  arma::rowvec new_eps_bar(num_comps), new_eps_var(num_comps);
  arma::mat new_mu_bar(dim, num_comps), new_mu_var(dim, num_comps);
  arma::cube new_cov_bar(dim, dim, num_comps), new_cov_var(dim, dim, num_comps);


  // Calculating new means

  new_eps_bar = OLD_WEIGHT * Theta.eps_bar + NEW_WEIGHT * p;
  new_mu_bar = OLD_WEIGHT * Theta.mu_bar + NEW_WEIGHT * mean;
  for (int k = 0; k < num_comps; k++) {
    new_cov_bar.slice(k) = OLD_WEIGHT * Theta.cov_bar.slice(k) + NEW_WEIGHT * cov.slice(k);
  }
  // Calculating new variances

  new_eps_var =
    OLD_WEIGHT * Theta.eps_var +
    OLD_WEIGHT * arma::square(Theta.eps_bar - new_eps_bar) +
    NEW_WEIGHT * arma::square(p - new_eps_bar);
  new_mu_var =
    OLD_WEIGHT * Theta.mu_var +
    OLD_WEIGHT * arma::square(Theta.mu_bar - new_mu_bar) +
    NEW_WEIGHT * arma::square(mean - new_mu_bar);

  for (int k = 0; k < num_comps; k++) {
    new_cov_var.slice(k) =
      OLD_WEIGHT * Theta.cov_var.slice(k) +
      OLD_WEIGHT * arma::square(Theta.cov_bar.slice(k) - new_cov_bar.slice(k)) +
      NEW_WEIGHT * arma::square(cov.slice(k) -  new_cov_bar.slice(k));
  }

  newTheta.eps_bar= new_eps_bar;
  newTheta.eps_var= new_eps_var;
  newTheta.mu_bar= new_mu_bar;
  newTheta.mu_var= new_mu_var;
  newTheta.cov_bar= new_cov_bar;
  newTheta.cov_var= new_cov_var;

  return newTheta;
}

void ChangeLabels(int t,
                  arma::uvec perms,
                  arma::mat &p,
                  arma::cube &means,
                  arma::field<arma::cube> &covs) {
  arma::rowvec old_p = p.row(t);
  arma::mat old_mean = means.slice(t);
  arma::cube old_cov = covs(t);
  for (int k = 0; k < num_comps; k++) {
    p.row(t)(k) = old_p(perms(k)-1);
    means.slice(t).col(k) = old_mean.col(perms(k)-1);
    covs(t).slice(k) = old_cov.slice(perms(k)-1);
  }
}















