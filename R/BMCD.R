BMCD <- function(distances, max_p, G,
                 burn = 1000, iters = 5000,
                 modelNames = c("UU", "US", "UD", "EU", "ES", "ED"),
                 parallel = FALSE, cores) {
  match.arg(modelNames, several.ok = TRUE)

  if ((isTRUE(parallel) & missing(cores))) {
    stop("Please provide the number of cores for parallelization (must be an integer greater than 0)")
  } else if (isFALSE(parallel) & !missing(cores)) {
    stop("Number of cores should not be specified when parallel is set to FALSE")
  }

  if (missing(cores)) {
    cores = 0
  }

  allModels <- c("UU", "US", "UD", "EU", "ES", "ED")
  modelIndices <- match(modelNames, allModels)

  BMDS_out <- BMDS(distances, max_p,
                   parallel = parallel, cores = cores)
  p <- which.min(BMDS_out$mdsics)
  X_est <- BMDS_out$X[[p]]
  sigmasq_est <- BMDS_out$sigma_sq[[p]]

  # Running BMCD
  BMCD_out <- RunBMCD(distances, X_est, G, sigmasq_est, burn, iters, modelIndices,
                      parallel = parallel, cores = cores)

  BMCD_out
}
