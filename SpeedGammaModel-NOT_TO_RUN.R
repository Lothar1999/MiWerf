#' This file contains that bit of code, that needs to be inserted 
#' in to the route inference script, when tracing the 
#' Estelle Metropolis function
#' 
#' Authors: Mike Werfeli, Peter Ranacher and Felix Liechti
#'
#' Just copy the code between the curly brackets {}
#' 

speedGammaModel<- function (beta, dt)
  {
  if (any(dt <= 0)) 
    stop("Data not ordered in time")
  if (!is.matrix(beta)) 
    beta <- t(beta)
  if (ncol(beta) == 2) {
    estelle.logpb <- function(x, z) {
      spd <- pmax.int(trackDist2(x, z), 1e-06)/dt
      bird_vec_angle <- directions(x, z)
      
      logpb_list <- wind_model(x,z, bird_vec_angle, dt)
      return(logpb_list)
    }
  }
  if (ncol(beta) == 3) {
    estelle.logpb <- function(x, z) {
      spd <- pmax.int(trackDist2(x, z), 1e-06)/dt
      angle <- trackBearingChange2(x, z)
      dgamma(spd, beta[, 1L], beta[, 2L], log = TRUE) + 
        dnorm(angle, 0, beta[, 3L], log = TRUE)
    }
  }
  stella.logpb <- function(x) {
    spd <- pmax.int(trackDist(x), 1e-06)/dt
    dgamma(spd, beta[, 1L], beta[, 2L], log = TRUE)
  }
  list(estelle.logpb = estelle.logpb, stella.logpb = stella.logpb, 
       beta = beta, dt = dt)
}


