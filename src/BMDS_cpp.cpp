// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
//#define ARMA_64BIT_WORD 1
#include "RcppArmadillo.h"
#include "Rcpp.h"

using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


//////////////////////////////////////////////////////////////////////////
//' @name distRcpp
//' @description calculate Euclidean distances between rows of matrix X
//' @title calculate Euclidean distances
//' @usage distRcpp(X)
//' @param X data matrix
//' @return distance matrix
//' @examples
//' x <- matrix(rnorm(100), nrow = 5)
//' distRcpp(x)
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix distRcpp(const NumericMatrix X) {
  const int n = X.nrow();
  NumericMatrix D(n, n);

  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      double dist = 0.0;
      for (int k = 0; k < X.ncol(); ++k) {
        const double tmp = X(i, k) - X(j, k);
        dist += tmp * tmp;
      }
      D(i, j) = sqrt(dist);
      D(j, i) = D(i, j);
    }
  }
  return D;
}

arma::mat eigVec(arma::mat xx){
  int p = xx.n_rows;
  arma::mat eig_V;
  arma::vec eigT;
  eig_sym(eigT,eig_V,xx);
  arma::mat evac(p,p);
  for(int i=0;i<p;i++)
    for(int j=0;j<p;j++)
      evac(i,j) = eig_V(i,p-1-j);
  return evac;
}



arma::mat rg(const arma::mat& w,double scale,int num) {
  int p = w.n_rows;
  arma::mat lambda(num,p,arma::fill::zeros);


  for (int j=0;j<num;j++){
    for(int i = 0; i<p;i++){
      lambda(j,i)=arma::conv_to<double>::from(arma::randg<arma::mat>(1,1,arma::distr_param(w(i),scale)));
    }
  }
  return mean(lambda.rows(0,num-1),0);
}


//////////////////////////////////////////////////////////////////////////
//' @name bmdsMCMC
//' @description run MCMC algorithm given in Oh and Raftery (2001) and return  posterior samples of
//' parameters as well as
//'  object configuration and other parameter estimates, for a given number of dimensions p
//' @title MCMC for Bayesian multidimensional scaling
//' @usage bmdsMCMC(DIST,p,nwarm = 1000,niter = 5000)
//' @param  DIST symmetric matrix of dissimilarity measures between objects
//' @param p number of dimensions of object configuration
//' @param nwarm number of iterations for  burn-in period in MCMC (default=1000)
//' @param niter number of MCMC iterations after burn-in period (default=5000)
//' @return  A list of  MCMC results
//' \describe{
//'   \item{x_bmds}{n by p matrix of object configuration that minimizes the sum of squares of residuals(SSR),
//'   where n is the number of objects, i.e., n=nrow(DIST)}
//'   \item{cmds}{n by p matrix of object configuration from the classical multidimensional scaling of Togerson(1952)}
//'   \item{minSSR}{minimum of sum of squares of residuals between the observed dissimilarities and
//'   the estimated Euclidean distances for pairs of objects}
//'   \item{minSSR_id}{index of the iteration corresponding to minimum SSR}
//'   \item{stress}{ STRESS computed from minSSR }
//'   \item{e_sigma}{ posterior mean of \eqn{\sigma^2}}
//'   \item{var_sigma}{ posterior variance of \eqn{\sigma^2}}
//'   \item{SSR.L}{niter dimensional vector of posterior samples of SSR}
//'   \item{lam.L}{niter by  p matrix of posterior samples of elements of \eqn{\Lambda}}
//'   \item{sigma.L}{niter dimensional vector of posterior samples of \eqn{\sigma^2} }
//'   \item{del.L}{niter by \eqn{n(n-1)/2} matrix of posterior samples of \eqn{\delta}, p-dimensional Euclidean distances
//'   between pairs of objects}
//'   }
//' @references Oh, M-S., Raftery A.E. (2001). Bayesian Multidimensional Scaling and Choice of Dimension,
//' Journal of the American Statistical Association, 96, 1031-1044.
//' @examples
//' \donttest{
//' data(cityDIST)
//' result=bmdsMCMC(cityDIST,p=3)
//' }
// [[Rcpp::export]]

extern "C" SEXP bmdsMCMC(SEXP DIST, SEXP p, int nwarm= 1000, int niter=5000){
  Rcpp::NumericMatrix DISTc(DIST);
  int pc = Rcpp::as<int> (p);

  int n = DISTc.nrow();
  int m = n*(n-1)/2;
  double scale = 2.38*2.38;
  int totiter = nwarm+niter;
  int xp = n*pc;

  Rcpp::Function RDist("dist");
  Rcpp::Function as_matrix("as.matrix");

  int i,j,iter,iterID;
  double s_dsq,xmean,sigma,stress,s_res;
  double pralpha_sig,prbeta_sig,alpha_lam,ssr;
  double s_sigma,sq_sigma,rmin_ssr,quad,quad_st,e_sigma;
  double t1, t1_st,s_temp,rl_gw,rl_fw,rl_f,rl_g,rl_w;

  double sig_old,cdvar_sig,sig_new,pst_alpha_lam;

  arma::vec pst_beta_lam(pc),s1(pc),rg(pc),lowertriD(m);
  arma::vec normR(1), unifR(1);
  Rcpp::NumericVector T1(n),T2(n),T3(n),T4(n);

  s_dsq=0;
  for(i=0; i<n; i++)
    s_dsq += sum(DISTc.row(i)*DISTc.row(i));
  arma::mat rI_p(pc,pc);
  rI_p.diag().fill(1);


  Rcpp::Environment stats("package:stats");
  Rcpp::Function cmdscale = stats["cmdscale"];
  Rcpp::NumericMatrix tmp_x = cmdscale(DISTc,Rcpp::Named("k",pc));
  Rcpp::NumericMatrix XX = Rcpp::clone(tmp_x);
  for(i=0; i < pc ;i++){
    xmean = mean(XX.column(i));
    XX.column(i) = XX.column(i)-xmean;
  }
  Rcpp::NumericMatrix XXX = Rcpp::clone(XX);
  arma::mat x(XXX.begin(),XXX.nrow(),XXX.ncol(),false);

  arma::mat Sig_x = x.t()*x/n;

  arma::mat eig_V=eigVec(Sig_x);

  x = x*eig_V;
  Sig_x = eig_V.t()*Sig_x*eig_V;

  XXX = Rcpp::clone(wrap(x));
  arma::mat x_int = x;

  // Rcpp::NumericMatrix delta = distRcpp(wrap(x_int));
  NumericMatrix delta = as_matrix(RDist(x_int));

    s_res =0;
  for(i=0; i<n; i++)
    s_res += sum((delta.row(i)- DISTc.row(i))*(delta.row(i)- DISTc.row(i)));
  s_res = s_res/2.0;
  sigma = s_res/m;
  stress = sqrt(s_res/s_dsq);


  //initial lambda
  arma::vec lambda = 1/Sig_x.diag();

  // set parameters of the prior

  pralpha_sig = 5.0;
  prbeta_sig = (pralpha_sig-1.0)*sigma;
  alpha_lam = 0.5;
  arma::vec beta_lam = 0.5*Sig_x.diag();

  // START MCMC
  s_sigma = 0;
  sq_sigma = 0;
  rmin_ssr = (stress*stress)*s_dsq;

  arma::mat x_keep = x;
  arma::vec SSR_samples(niter);
  arma::mat lam_samples(niter,pc);
  arma::vec sigma_samples(niter);
  arma::mat x_samples(niter,xp);

  arma::mat x_sv=x;
  arma::mat cd_var(pc,pc),cd_var_inv(pc,pc),cd_sig(pc,pc);
  arma::vec x_old(pc),cmd(pc),x_new(pc),rn(pc);
  for(iter=0; iter<totiter;iter++){
    x_sv = x;
    for(i=0;i<n; i++){
      cd_var = scale * sigma/(n-1)*rI_p;

      cd_var_inv.diag() =1/cd_var.diag();
      cd_sig = sqrt(scale*sigma/(n-1))*(rI_p);
      x_old = conv_to<vec>::from(x.row(i));

      cmd = x_old;
      rn = rnorm(pc);
      arma::vec yrn (rn.begin(),rn.size(),false);

      quad = 0;
      quad_st = 0;
      for(j=0;j<pc;j++){
        x_new(j) = yrn(j)*cd_sig.diag()(j)+cmd(j);
        quad += x_old(j)*x_old(j)*lambda(j);
        quad_st += x_new(j)*x_new(j)*lambda(j);
      }

      t1 = 0;
      t1_st = 0;
      s_temp = 0;

      arma::mat x_prop = x;
      x_prop.row(i) = x_new.as_row();
      arma::mat D_old = square(x.each_row() - x.row(i));
      arma::mat D_new = square(x_prop.each_row() - x_prop.row(i));
      arma::vec delta_i = sqrt(sum(D_old, 1));
      arma::vec delta_i_new = sqrt(sum(D_new, 1));

      Rcpp::NumericVector DISTc_row_i = DISTc(i, _);
      arma::rowvec DISTc_row_i_arma = arma::rowvec(DISTc_row_i.begin(), DISTc_row_i.size());

      t1 = (-0.5/sigma)*sum(square(delta_i.t() - DISTc_row_i_arma));
      t1_st = (-0.5/sigma)*sum(square(delta_i_new.t() - DISTc_row_i_arma));


      delta_i.shed_row(i);
      delta_i_new.shed_row(i);
      arma::vec a = log(arma::normcdf(delta_i_new/ sqrt(sigma)) / arma::normcdf(delta_i/ sqrt(sigma)));
      s_temp = sum(a);

      rl_gw = 0;
      rl_fw = t1_st - t1 - s_temp - 0.5*(quad_st-quad);
      unifR = Rcpp::runif(1);
      if(log(unifR(0)) < (rl_fw - rl_gw))
        x_old = x_new;
      x.row(i) = conv_to<rowvec>::from(x_old);

    }
    for(i=0; i<pc;i++)
      x.col(i) = x.col(i)-mean(x.col(i));

    Sig_x = x.t()*x/n;
    eig_V=eigVec(Sig_x);


    x = x*eig_V;
    Sig_x = eig_V.t()*Sig_x*eig_V;

    // generate sigma
    // delta = distRcpp(wrap(x));
    delta = as_matrix(RDist(wrap(x)));
    ssr =0;
    for(i=0; i<n; i++)
      ssr += sum((delta.row(i)- DISTc.row(i))*
        (delta.row(i)- DISTc.row(i)));
    ssr = ssr/2.0;
    sig_old = sigma;


    cdvar_sig = scale *2*s_res*s_res /((m-2)*(m-2)*(m-4));
    sig_new = -1;
    while(sig_new<0){
      normR = Rcpp::rnorm(1);
      sig_new = normR(0)*sqrt(cdvar_sig)+sig_old;
    }

    T1 = wrap(delta / sqrt(sig_new));
    T2 = wrap(delta / sqrt(sig_old));

    rl_f = -sum(log(pnorm(T1)/pnorm(T2)))/2-
      (0.5*s_res+prbeta_sig)*(1/sig_new-1/sig_old)-
      (0.5*m+pralpha_sig+1)*log(sig_new/sig_old);
    T3 = wrap(sig_new/sqrt(cdvar_sig));
    T4 = wrap(sig_old/sqrt(cdvar_sig));

    rl_g = sum(log(pnorm(T3)/pnorm(T4)));
    rl_w = rl_f-rl_g;

    unifR = Rcpp::runif(1);
    if(log(unifR(0))<=rl_w){
      sig_old = sig_new;
    }
    sigma = sig_old;
    //  generate Lambda
    for(i=0; i<pc; i++){
      s1(i)=0;
      for(j=0;j<n; j++)
        s1(i) += x(j,i)*x(j,i);
    }
    pst_beta_lam = 0.5*s1+beta_lam;
    pst_alpha_lam = 0.5*n+alpha_lam;

    for(i=0;i<pc;i++){
      rg(i) = rgamma(1,pst_alpha_lam)(0);
    }

    arma::vec yrg (rg.begin(),rg.size(),false);

    lambda=yrg/pst_beta_lam;
    if(iter >= nwarm){
      iterID = iter-nwarm;
      s_sigma = s_sigma+sigma;
      sq_sigma = sq_sigma+sigma*sigma;

      ssr =0;
      for(i=0; i<n; i++)
        ssr += sum((delta.row(i)- DISTc.row(i))*
          (delta.row(i)- DISTc.row(i)));
      ssr = ssr/2.0;

      if(ssr < rmin_ssr){
        rmin_ssr = ssr;
        x_keep = x;
      }

      SSR_samples(iterID) = ssr;
      arma::rowvec tempV(xp);
      for(i=0;i<pc;i++)
        for(j=0;j<n;j++)
          tempV(i*n+j) = x(j,i);

      x_samples.row(iterID)= conv_to<rowvec>::from(tempV);

      lam_samples.row(iterID)=conv_to<rowvec>::from(lambda);

      sigma_samples(iterID)=sigma;
      int tempcount=0;
      for(i=0; i<n;i++){
        for(j=(i+1);j<n;j++){
          lowertriD(tempcount) = delta(i,j);
          tempcount++;
        }
      }
    }
    // end MCMC

  }
  e_sigma = s_sigma/ niter;
  return Rcpp::List::create(
    Rcpp::Named("x_bmds") = x_keep,
    Rcpp::Named("e_sigma")= e_sigma);
}

