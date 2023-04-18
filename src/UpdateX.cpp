// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "commonvars.h"
#include "helpers.h"

using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat UpdateX(arma::mat obs_distances, arma::mat current_X, double sigmasq, arma::rowvec z,
                  arma::mat means, arma::cube covs) {

  arma::mat new_X = current_X;
  arma::mat I_p(dim, dim, arma::fill::eye);
  arma::mat cd_var(dim, dim), cd_sigma(dim, dim), T_inv(dim, dim);
  arma::vec x_old(dim), rand_vec(dim), x_new(dim), mu(dim), center_new_x(dim), center_old_x(dim);

  const int scale = 2.38*2.38;
  double Q1, Q1_st;
  double delta_ij, delta_st_ij, delta_sq_ij, delta_sq_st_ij;
  double norm_d, norm_d_st, norm_ratio, t1, t1_st, s_temp;
  double rl_fw;
  int comp;
  int accepted = 0;

  for (int i = 0; i < n; i++) {

    comp = int (z(i)-1);
    T_inv = inv(covs.slice(comp));
    mu = means.col(comp);
    cd_var = scale * (sigmasq / (n-1)) * I_p;
    cd_sigma = arma::sqrt(cd_var);
    x_old = current_X.row(i).t();

    rand_vec = as<arma::vec>(wrap(rnorm(dim)));

    for (int j = 0; j < dim; j++) {
      x_new(j) = rand_vec(j) * cd_sigma.diag()(j) + x_old(j);
      center_new_x(j) = x_new(j) - mu(j);
      center_old_x(j) = x_old(j) - mu(j);
    }

    Q1 = as<double>(wrap(center_old_x.t() * T_inv * center_old_x));
    Q1_st = as<double>(wrap(center_new_x.t() * T_inv * center_new_x));
    t1 = 0;
    t1_st = 0;
    s_temp = 0;

    arma::mat x_prop = current_X;
    x_prop.row(i) = x_new.as_row();
    arma::mat D_old = square(current_X.each_row() - current_X.row(i));
    arma::mat D_new = square(x_prop.each_row() - x_prop.row(i));
    arma::vec delta_i = sqrt(sum(D_old, 1));
    arma::vec delta_i_new = sqrt(sum(D_new, 1));

    arma::rowvec DISTc_row_i_arma = obs_distances.row(i);

    t1 = (-0.5/sigmasq)*sum(square(delta_i.t() - DISTc_row_i_arma));
    t1_st = (-0.5/sigmasq)*sum(square(delta_i_new.t() - DISTc_row_i_arma));

    delta_i.shed_row(i);
    delta_i_new.shed_row(i);
    arma::vec a = log(arma::normcdf(delta_i_new/ sqrt(sigmasq)) / arma::normcdf(delta_i/ sqrt(sigmasq)));
    s_temp = sum(a);

    // for (int j = 0; j < n; j++) {
    //   if (j != i) {
    //     delta_sq_ij = 0;
    //     delta_sq_st_ij = 0;
    //     for (int d = 0; d < dim; d++) {
    //       delta_sq_ij += (current_X(j,d) - x_old(d)) * (current_X(j,d) - x_old(d));
    //       delta_sq_st_ij += (current_X(j,d) - x_new(d)) * (current_X(j,d) - x_new(d));
    //     }
    //
    //     delta_ij = sqrt(delta_sq_ij);
    //     delta_st_ij = sqrt(delta_sq_st_ij);
    //     norm_d = R::pnorm(delta_ij / sqrt(sigmasq), 0, 1, true, false);
    //     norm_d_st = R::pnorm(delta_st_ij / sqrt(sigmasq), 0, 1, true, false);
    //     norm_ratio = log(norm_d / norm_d_st);
    //     t1 = t1 - (1/(2*sigmasq)) * (delta_ij - obs_distances(i, j)) * (delta_ij - obs_distances(i, j));
    //     t1_st = t1_st - (1/(2*sigmasq)) * (delta_st_ij - obs_distances(i, j)) * (delta_st_ij - obs_distances(i, j));
    //     s_temp = s_temp - norm_ratio;
    //
    //   }
    // }

    rl_fw = t1_st - t1 - s_temp - 0.5*(Q1_st - Q1);

    if (log(R::runif(0, 1)) < rl_fw) {
      new_X.row(i) = x_new.as_row();
      accepted += 1;
    }
  }
  return new_X;
}


double UpdateSigmasq(arma::mat obs_distances, arma::mat current_X, double s_res, double sigmasq) {
  double sigmasq_out, sig_old, sig_new, norm_sig, rl_f, rl_g, rl_w;
  double scale = 2.38*2.38, cdvar_sig;
  double m = n*(n-1)/2;
  Rcpp::Function RDist("dist");
  Rcpp::Function as_matrix("as.matrix");
  Rcpp::NumericVector T1(n), T2(n), T3(n), T4(n);
  Rcpp::NumericMatrix delta = as_matrix(RDist(current_X));
  double ssr = CalcSSR(delta, wrap(obs_distances));

  sig_old = sigmasq;
  cdvar_sig = scale * 2 * s_res * s_res / ((m-2)*(m-2)*(m-4));
  sig_new = -1;
  while(sig_new<0){
    norm_sig = R::rnorm(0, 1);
    sig_new = norm_sig*sqrt(cdvar_sig)+sig_old;
  }

  T1 = wrap(delta / sqrt(sig_new));
  T2 = wrap(delta / sqrt(sig_old));
  rl_f = -sum(log(pnorm(T1)/pnorm(T2)))/2-
    (0.5*ssr+sigmasq_prior_scale)*(1/sig_new-1/sig_old)-
    (0.5*m+sigmasq_prior_shape+1)*log(sig_new/sig_old);
  T3 = wrap(sig_new/sqrt(cdvar_sig));
  T4 = wrap(sig_old/sqrt(cdvar_sig));

  rl_g = sum(log(pnorm(T3)/pnorm(T4)));
  rl_w = rl_f-rl_g;

  if(log(R::runif(0, 1)) <= rl_w){
    sig_old = sig_new;
  }
  sigmasq_out = sig_old;
  return sigmasq_out;
}


arma::mat TransformX(arma::mat X, arma::mat X_star) {
  arma::mat new_X(n, dim);
  arma::mat I_mat(n, n, arma::fill::eye), J(n, n), C(n, n), XT;
  arma::mat U, V, T, little_t;
  arma::vec s;
  arma::vec vec_1(n, arma::fill::ones);
  arma::mat test = vec_1 * vec_1.t();
  double new_n = (double) n;
  J = I_mat - (1/new_n) * (vec_1 * vec_1.t());
  C = X_star.t() * J * X;
  svd(U, s, V, C);
  T = V * U.t();
  XT = X * T;
  little_t = (1/new_n) * (X_star - (XT)).t() * vec_1;
  new_X = XT + (vec_1 * little_t.t());
  return new_X;
}
