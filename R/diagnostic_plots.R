#' @importFrom coda mcmc
#' @importFrom mcmcplots mcmcplot

diagnostics <- function(BMCD_object) {
  output <- BMCD_object
  iters <- nrow(output$p)
  num_clusts <- output$G
  chosen_dim <- ncol(output$BMDS_X)

  p_mat <- matrix(NA, nrow = iters, ncol = num_clusts)
  pnames <- c(paste("weight[", 1:output$G, "]", sep=""))
  p_mat <- output$p
  colnames(p_mat) <- pnames

  #Means
  num_means <- num_clusts * chosen_dim
  mean_mat <- matrix(NA, nrow = iters, ncol = num_means)
  index <- 1
  mean_names <- c()
  for (i in 1:chosen_dim) {
    for (j in 1:num_clusts) {
      name <- paste("mean[", i, ",", j, "]", sep = "")
      elements <- output$means[i, j, , drop = TRUE]
      mean_mat[,index] <- elements
      mean_names <- c(mean_names, name)
      index <- index + 1
    }
  }
  colnames(mean_mat) <- mean_names

  #Covariances
  num_covs <- chosen_dim^2 * num_clusts
  cov_mat <- matrix(NA, nrow = iters, ncol = num_covs)
  cov_names <- c()
  index <- 1
  for(i in 1:num_clusts) {
    for (j in 1:chosen_dim) {
      for (k in 1:chosen_dim) {
        name <- paste("covariance[", i, ",", j, ",", k, "]", sep = "")
        elements <- sapply(output$covs, FUN = function(x) {x[,,i][j,k]})
        cov_mat[,index] <- elements
        cov_names <- c(cov_names, name)
        index <- index + 1
      }
    }
  }
  colnames(cov_mat) <- cov_names
  for(i in 1:num_covs) {
    if (all(cov_mat[,i] == 0)) {
      cov_mat[,i] <- NA
    }
  }
  colnames(cov_mat) <- cov_names

  new_cov_matrix <- cov_mat[, !colSums(is.na(cov_mat))]


  # Plot
  full_mat <- cbind(p_mat, mean_mat, new_cov_matrix)
  full_coda <- coda::mcmc(full_mat)
  mcmcplot(full_coda, heading = "Diagnostic Plots")
}
