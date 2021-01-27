function (model, proposal.x, proposal.z, x0 = NULL, z0 = NULL, 
        iters = 1000L, thin = 10L, chains = 1L, verbose = interactive()) 
{
  if (is.null(x0)) 
    x0 <- model$x0
  if (is.null(z0)) 
    z0 <- model$z0
  x0 <- rep(if (is.list(x0)) x0 else list(x0), length.out = chains)
  z0 <- rep(if (is.list(z0)) z0 else list(z0), length.out = chains)
  n <- nrow(x0[[1]])
  m <- ncol(x0[[1]])
  logpx <- model$logpx
  logpz <- model$logpz
  logpb <- model$estelle.logpb
  fixedx <- model$fixedx
  ch.xs <- vector(mode = "list", chains)
  ch.zs <- vector(mode = "list", chains)
  for (k1 in 1:chains) {
    ch.x <- array(0, c(n, m, iters))
    ch.z <- array(0, c(n - 1L, 2L, iters))
    x1 <- x0[[k1]]
    z1 <- z0[[k1]]
    dimnames(x1) <- NULL
    dimnames(z1) <- NULL
    logp.x1 <- logpx(x1)
    logp.z1 <- logpz(z1)
    logp.b1 <- logpb(x1, z1)
    logDF<<- logp.x1
    k2 <- 0
    if (verbose) {
      cat("iter ", sprintf("%6d", k2))
      flush.console()
    }
    for (k2 in 1:iters) {
      if (verbose && k2%%10 == 0) {
        cat("\b\b\b\b\b\b")
        cat(sprintf("%6d", k2))
        flush.console()
      }
      for (k3 in 1:thin) {
        x2 <- proposal.x(x1)
        x2[fixedx, ] <- x1[fixedx, ]
        logp.x2 <- logpx(x2)
        x <- x1
        x[c(1L, n), ] <- x2[c(1L, n), ]
        logp.b2 <- logpb(x, z1)
        if (!fixedx[1L]) {
          logp1 <- logp.x1[1L] + logp.b1[1L]
          logp2 <- logp.x2[1L] + logp.b2[1L]
          if (logp2 - logp1 > log(runif(1))) {
            x1[1L, ] <- x2[1L, ]
            logp.x1[1L] <- logp.x2[1L]
            logp.b1[1L] <- logp.b2[1L]
          }
        }
        if (!fixedx[n]) {
          logp1 <- logp.x1[n] + logp.b1[n - 1L]
          logp2 <- logp.x2[n] + logp.b2[n - 1L]
          if (logp2 - logp1 > log(runif(1))) {
            x1[n, ] <- x2[n, ]
            logp.x1[n] <- logp.x2[n]
            logp.b1[n - 1L] <- logp.b2[n - 1L]
          }
        }
        for (rb in 2:3) {
          is <- seq.int(rb, n - 1L, by = 2L)
          x <- x1
          x[is, ] <- x2[is, ]
          logp.b2 <- logpb(x, z1)
          logp1 <- logp.x1[is] + logp.b1[is - 1L] + 
            logp.b1[is]
          logp2 <- logp.x2[is] + logp.b2[is - 1L] + 
            logp.b2[is]
          accept <- is[logp2 - logp1 > log(runif(length(is)))]
          x1[accept, ] <- x[accept, ]
          logp.x1[accept] <- logp.x2[accept]
          logp.b1[accept] <- logp.b2[accept]
          logp.b1[accept - 1L] <- logp.b2[accept - 1L]
        }
        logDF <<- cbind(logDF, logp.x1)
        z2 <- proposal.z(z1)
        logp.z2 <- logpz(z2)
        logp.b2 <- logpb(x1, z2)
        logp1 <- logp.z1 + logp.b1
        logp2 <- logp.z2 + logp.b2
        accept <- logp2 - logp1 > log(runif(n - 1L))
        z1[accept, ] <- z2[accept, ]
        logp.z1[accept] <- logp.z2[accept]
        logp.b1[accept] <- logp.b2[accept]
      }
      ch.x[, , k2] <- x1
      ch.z[, , k2] <- z1
    }
    ch.xs[[k1]] <- ch.x
    ch.zs[[k1]] <- ch.z
    if (verbose) 
      cat("\n")
  }
  list(model = model, x = ch.xs, z = ch.zs)
}
