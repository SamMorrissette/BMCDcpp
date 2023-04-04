createOutput <- function(bmcd_obj, init_X, burn, iters) {
  ind <- (burn+1):iters
  out <- list(BMDS_X = init_X,
              X = dpobj$X[,,ind],
              sigmasq = bmcd_obj$sigmasq[ind,],
              means = bmcd_obj$means[,,ind],
              covs = bmcd_obj$covs[ind,],
              z = bmcd_obj$z[ind,],
              class_probs = bmcd_obj$class_probs[,,ind],
              p = bmcd_obj$p[ind,],
              perms = bmcd_obj$perms[ind,])
  return(out)
}


RunBMCD <- function(distances, init_X, init_sigmasq, burn, iters, modelIndices, parallel, cores) {
  if (parallel == FALSE) {
    # Run MCMC one model at a time
    allModels <- vector("list", length = length(modelIndices))
    for (i in 1:length(modelIndices)) {
      bmcd_obj <- BMCD_MCMC(distances, init_X, init_sigmasq, iters, modelIndices[i])
      allModels[[i]] <- createOutput(bmcd_obj, init_X, burn, iters)
    }
  } else if (parallel == TRUE) {
    if (!foreach::getDoParRegistered()) {
      doParallel::registerDoParallel(cores=cores)
    }
    allModels <- foreach::foreach(j=1:length(modelIndices), .packages ="BMCD") %dopar% {
      bmcd_obj <- BMCD_MCMC(distances, init_X, init_sigmasq, iters, modelIndices[j])
      output <- createOutput(bmcd_obj, init_X, burn, iters)
    }
  }

  return(allModels)
}
