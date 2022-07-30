dyn.load("../vecchia")


vecchia_profbeta_loglik_grad_info_matern <- function(covparms, locs, y, X, NNarray){
  
  n <- length(y)
  m <- ncol(NNarray)
  dim <- ncol(locs)
  p <- ncol(X)
  nparms <-length(covparms)

  ll <- 0.0
  betahat <- rep(0,p)

  grad <- rep(0,nparms)
  info <- rep(0, nparms*nparms)
  
  a <- .C("vecchia_profbeta_likelihood_matern",
          ll = as.double(ll),
	  betahat = as.double(betahat),
          grad = as.double(grad),
	  info = as.double(info),
          as.double(covparms),
	  as.integer(nparms),
          as.double(y),
	  as.integer(n),
          as.double(locs),
          as.integer(dim),
	  as.double(X),
	  as.integer(p),
          as.integer(t(NNarray)),
          as.integer(m),
          NAOK = TRUE)

  info_mat <- matrix(a$info, nrow  = nparms, ncol = nparms)
  return(list(loglik = a$ll, grad = a$grad, info = info_mat, betahat = a$betahat))
}



vecchia_profbeta_loglik_grad_info_matern_multi <- function(covparms, locs, y, X, NNarray){
  
  n <- length(y)
  m <- ncol(NNarray)
  dim <- ncol(locs)
  p <- ncol(X)
  nparms <-length(covparms)

  ll <- 0.0
  betahat <- rep(0,p)

  grad <- rep(0,nparms)
  info <- rep(0, nparms*nparms)
  
  a <- .C("vecchia_profbeta_likelihood_matern_multi",
          ll = as.double(ll),
	  betahat = as.double(betahat),
          grad = as.double(grad),
	  info = as.double(info),
          as.double(covparms),
	  as.integer(nparms),
          as.double(y),
	  as.integer(n),
          as.double(locs),
          as.integer(dim),
	  as.double(X),
	  as.integer(p),
          as.integer(t(NNarray)),
          as.integer(m),
          NAOK = TRUE)

  info_mat <- matrix(a$info, nrow  = nparms, ncol = nparms)
  return(list(loglik = a$ll, grad = a$grad, info = info_mat, betahat = a$betahat))
}


