#' Run Bayesian Multidimensional Scaling
#' @importFrom foreach %dopar% foreach
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#'
#' @description Represent distances between objects in $p$-dimensional Euclidean space.
#' @param distances A (symmetric) matrix of distance values returned by \code{dist}, or a full matrix of values.
#' @param max_dim The maximum dimension to be estimated by the algorithm. The BMDS algorithm  will be
#'  run for all dimensions between 1 and \code{max_dim}, and the dimension resulting in the lowest MDSIC value will be chosen as the
#'  final object representation.
#' @param parallel Boolean indicating whether or not parallel computing should be utilized when running
#' multiple MCMC chains. If TRUE, the number of cores must be specified.
#' @param cores The number of CPU cores to be used when parallel is set to TRUE.
#'
#' @return A list containing the following items.
#' @return \item{X}{A list of matrices, with each element corresponding to a different dimensional representation of the objects.}
#' @return \item{sigma_sq}{A list of estimated measurement errors, with each element of the list corresponding to the representation of the objects given in the returned list \code{X}.}
#' @return \item{MDSIC}{A vector of Multidimensional Scaling Information Criterion values (MDSIC). Each element of the vector corresponds to the representation of objects given in the returned list \code{X}.}
#' @examples
#' BMDS(cityDIST, max_dim = 6)
#' @export

BMDS <- function(distances, max_dim, parallel = FALSE, cores) {
  if(parallel == TRUE & !getDoParRegistered()) {
    doParallel::registerDoParallel(cores=cores)
  }

  n <- nrow(distances)

  X <- vector("list", length = max_dim)
  sigma_sq <- vector("list", length = max_dim)

  if (parallel == FALSE) {
    bmds_burn = 1000
    bmds_iter = 5000
    for (i in 1:max_dim) {
      temp_bmds <- bmdsMCMC(DIST = distances, p = i, nwarm = bmds_burn, niter = bmds_iter)
      X[[i]] <- temp_bmds$x_bmds
      sigma_sq[[i]] <- temp_bmds$e_sigma
      print(i)
    }
  } else if (parallel == TRUE) {
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    out_list <- foreach::foreach(j=1:max_dim, .packages ="BMCDcpp", .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
      bmds_burn = 1000
      bmds_iter = 5000
      output <- bmdsMCMC(DIST = distances, p = j, nwarm = bmds_burn, niter = bmds_iter)
      list(output$x_bmds, output$e_sigma)
    }
    X <- out_list[[1]]
    sigma_sq <- out_list[[2]]
  }
  MDSIC <- MDSIC(distances, X)

  return(list(X = X,
              sigma_sq = sigma_sq,
              MDSIC = MDSIC))
}
