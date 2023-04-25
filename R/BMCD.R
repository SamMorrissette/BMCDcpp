#' Run Bayesian Model-based Clustering with Dissimilarities
#'
#' @description The main function of the package that runs both the Bayesian Multidimensional Scaling and Bayesian Model-Based Clustering
#' algorithm.
#' @param distances A (symmetric) matrix of distance values returned by \code{dist}, or a full matrix of values.
#' @param max_dim The maximum dimension to be estimated by Bayesian Multidimensional Scaling (BMDS). The BMDS algorithm  will be
#'  run for all dimensions between 1 and \code{max_dim}, and the dimension resulting in the lowest MDSIC value will be chosen as the
#'  final object representation.
#' @param G The number of clusters to be estimated. G can be a single value or a range of values.
#' @param burn The number of burn-in iterations for the main clustering MCMC algorithm. The burn-in iterations will be deleted
#' from the beginning of the MCMC chain.
#' @param iters The total number of iterations for the main clustering MCMC algorithm.
#' @param model_names The type of model(s) to be fit to the data. One or more of the following models may be specified: \itemize{
#' \item UU - Unequal Unrestricted
#' \item US - Unequal Spherical
#' \item UD - Unequal Diagonal
#' \item EU - Equal Unrestricted
#' \item ES - Equal Spherical
#' \item ED - Equal Diagonal}
#' @param parallel Boolean indicating whether or not parallel computing should be utilized when running
#' multiple MCMC chains. If TRUE, the number of cores must be specified.
#' @param cores The number of CPU cores to be used when parallel is set to TRUE.
#'
#' @return A list containing the following items. All MCMC iterations are given without burn-in iterations.
#' @return \item{BMDS_X}{The object configuration matrix estimated by Bayesian Multidimensional Scaling.}
#' @return \item{X}{Iterations of the object configuration matrix.}
#' @return \item{sigma_sq}{Iterations of the estimated measurement error.}
#' @return \item{means}{Iterations of the estimated component means. Each column of the matrix corresponds to a different component.}
#' @return \item{covs}{Iterations of the estimated component covariance matrices. Each matrix slice corresponds to a different component.}
#' @return \item{z}{Iterations of the object classes. Each column corresponds to a different data point, whereas the rows correspond to the iterations.}
#' @return \item{class_probs}{Iterations of the estimated class probabilities. Each matrix slice corresponds to a different iteration, and each row of the matrix corresponds to the data points.}
#' @return \item{p}{Iterations of the estimated mixing probabilities. The rows correspond to iterations and the columns correspond to the components.}
#' @return \item{BIC}{The Bayesian Information Criterion (BIC) of the chosen model.}
#' @return \item{all_BIC}{The Bayesian Information Criterion (BIC) of all candidate models. The row names correspond to the number of clusters.}
#' @return \item{model_name}{The name (type) of the chosen model.}
#' @return \item{G}{The number of components of the chosen model.}
#' @examples
#' BMCD(cityDIST, max_dim = 4, G = 2:6)
#'
#' @export
BMCD <- function(distances, max_dim, G,
                 burn = 1000, iters = 5000,
                 model_names = c("UU", "US", "UD", "EU", "ES", "ED"),
                 parallel = FALSE, cores) {
  match.arg(model_names, several.ok = TRUE)

  distances <- as.matrix(distances)
  if (!isSymmetric(distances)) {
    stop("The distance matrix should be symmetrical.")
  }

  if (iters <= burn) {
    stop("The total number of iterations must be higher than the burn-in.")
  }
  if ((isTRUE(parallel) & missing(cores))) {
    stop("Please provide the number of cores for parallelization (must be an integer greater than 0).")
  } else if (isFALSE(parallel) & !missing(cores)) {
    stop("Number of cores should not be specified when parallel is set to FALSE.")
  }

  if (!all(G > 0)) {
    stop("The number of clusters must be a value or range of values greater than 0.")
  }

  if (missing(cores)) {
    cores = 0
  }

  allModels <- c("UU", "US", "UD", "EU", "ES", "ED")
  modelIndices <- match(model_names, allModels)

  BMDS_out <- BMDS(distances, max_dim,
                   parallel = parallel, cores = cores)
  p <- which.min(BMDS_out$MDSIC)
  X_est <- BMDS_out$X[[p]]
  sigmasq_est <- BMDS_out$sigma_sq[[p]]

  BIC_table <- matrix(nrow = length(G), ncol = length(model_names))
  colnames(BIC_table) <- model_names
  rownames(BIC_table) <- G

  # Running BMCD
  winning_bic = Inf
  for (num_clust in G) {
    BMCD_out <- RunBMCD(distances, X_est, num_clust, sigmasq_est, burn, iters, modelIndices,
                        parallel = parallel, cores = cores)
    names(BMCD_out) <- model_names
    bics <- sapply(BMCD_out, FUN = function(x) {
      x$BIC
    })
    min_bic = min(bics)
    bic_ind = which.min(bics)
    if (min_bic < winning_bic) {
      winning_bic = min_bic
      winning_name <- model_names[bic_ind]
      winning_clust <- num_clust
      final_model <- BMCD_out[[bic_ind]]
    }
    BIC_table[as.character(num_clust), ] <- bics
  }
  final_model$all_BIC = BIC_table
  final_model$model_name = winning_name
  final_model$G = winning_clust
  final_model
}
