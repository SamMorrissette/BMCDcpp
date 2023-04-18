#' @importFrom foreach %dopar% foreach
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#' @export

BMDS <- function(distances, max_p, parallel = FALSE, cores) {
  if(parallel == TRUE & !getDoParRegistered()) {
    doParallel::registerDoParallel(cores=cores)
  }
  n <- nrow(distances)

  X <- vector("list", length = max_p)
  sigma_sq <- vector("list", length = max_p)

  if (parallel == FALSE) {
    bmds_burn = 1000
    bmds_iter = 5000
    for (i in 1:max_p) {
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
    out_list <- foreach::foreach(j=1:max_p, .packages ="BMCD", .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
      bmds_burn = 1000
      bmds_iter = 5000
      output <- bmdsMCMC(DIST = distances, p = j, nwarm = bmds_burn, niter = bmds_iter)
      list(output$x_bmds, output$e_sigma)
    }
    X <- out_list[[1]]
    sigma_sq <- out_list[[2]]
  }
  mdsics <- MDSIC(distances, X)

  return(list(X = X,
              sigma_sq = sigma_sq,
              mdsics = mdsics))
}
