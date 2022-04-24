dyn.load("vecchia")


vecchia_meanzero_loglik_grad_info_matern <- function(covparms, locs, y, NNarray){
  
  n <- length(y)
  m <- ncol(NNarray)
  dim <- ncol(locs)
  nparms <-length(covparms)

  ll <- 0.0
  grad <- rep(0,nparms)
  info <- rep(0, nparms*nparms)
  
  a <- .C("vecchia_meanzero_likelihood_matern",
          ll = as.double(ll),
          grad = as.double(grad),
	  info = as.double(info),
          as.double(covparms),
	  as.integer(nparms),
          as.double(y),
	  as.integer(n),
          as.double(locs),
          as.integer(dim),
          as.integer(t(NNarray)),
          as.integer(m),
          NAOK = TRUE)

  info_mat <- matrix(a$info, nrow  = nparms, ncol = nparms)
  return(list(loglik = a$ll, grad = a$grad, info = info_mat))
}


vecchia_meanzero_loglik_grad_info_bivariate_matern <- function(covparms, locs, y, NNarray){
  
  n <- length(y)
  m <- ncol(NNarray)
  dim <- ncol(locs)
  nparms <-length(covparms)

  ll <- 0.0
  grad <- rep(0,nparms)
  info <- rep(0, nparms*nparms)
  
  a <- .C("vecchia_meanzero_likelihood_bivariate_matern",
          ll = as.double(ll),
          grad = as.double(grad),
	  info = as.double(info),
          as.double(covparms),
	  as.integer(nparms),
          as.double(y),
	  as.integer(n),
          as.double(locs),
          as.integer(dim),
          as.integer(t(NNarray)),
          as.integer(m),
          NAOK = TRUE)

  info_mat <- matrix(a$info, nrow  = nparms, ncol = nparms)
  return(list(loglik = a$ll, grad = a$grad, info = info_mat))
}




con <- function(d, covparms){
  
  nparms <-length(covparms)
  cons <- 0
  
  a <- .C("con",
	  cons = as.double(cons),
	  as.integer(d),
          as.double(covparms),
	  as.integer(nparms))

  return(a$cons)
}

