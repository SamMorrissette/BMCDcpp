#ifndef __UTILITIES__
#define __UTILITIES__

typedef struct Params Params;
typedef struct Theta_calcs Theta_calcs;

double CalcSSR(Rcpp::NumericMatrix mat1, Rcpp::NumericMatrix mat2);
arma::rowvec Arma_colSums(const arma::mat& x);

Rcpp::NumericMatrix distRcpp(Rcpp::NumericMatrix X);

void InitChains(arma::cube &X,
                arma::mat &p, arma::mat &z,
                arma::cube &means,
                arma::field<arma::cube> &covs);

arma::mat UpdateX(arma::mat obs_distances, arma::mat current_X,
                  double sigmasq, arma::rowvec z,
                  arma::mat means, arma::cube covs);

arma::mat TransformX(arma::mat X, arma::mat X_star);

double UpdateSigmasq(arma::mat obs_distances, arma::mat current_X, double s_res, double sigmasq);

arma::mat CalculateClassProbs(arma::rowvec p,
                              arma::mat mean,
                              arma::cube cov);

arma::rowvec UpdateClasses(arma::mat prob_mat);

arma::mat rdirichlet_cpp(int num_samples,
                         arma::vec alpha_m);

arma::rowvec UpdateWeights(arma::rowvec z);


// Label Switching Functions
Theta_calcs InitializeCenters(arma::mat p, arma::cube means, arma::field<arma::cube> covs);
arma::mat CostMatrix(Theta_calcs Theta, arma::rowvec p, arma::mat mean, arma::cube cov);
arma::uvec Hungarian(arma::mat cost_mat);
Theta_calcs UpdateCenters(Theta_calcs Theta, int t, arma::rowvec p, arma::mat mean, arma::cube cov);
void ChangeLabels(int t,
                  arma::uvec perms,
                  arma::mat &p,
                  arma::cube &means,
                  arma::field<arma::cube> &covs);

// Updating GMM parameters
Params UpdateTheta(arma::rowvec z, int modelIndex);
Params DrawUnequalUnrestricted(arma::rowvec z);
Params DrawUnequalSpherical(arma::rowvec z);
Params DrawUnequalDiagonal(arma::rowvec z);
Params DrawEqualUnrestricted(arma::rowvec z);
Params DrawEqualSpherical(arma::rowvec z);
Params DrawEqualDiagonal(arma::rowvec z);

double CalcBIC(arma::cube X, arma::mat p, arma::cube means, arma::field<arma::cube> covs,
               arma::cube class_probs,
               const int G,
               int modelIndex,
               int burn, int iters);
#endif // __UTILITIES__
