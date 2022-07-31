library(matrixcalc)

makeSymm <- function(m) {
   m[upper.tri(m)] <- t(m)[upper.tri(m)]
   return(m)
}

is.cnd <- function(M){
	p <- nrow(M)
	new <- matrix(0,p,p)
	for(i in 1:p){
	    for(j in 1:p){
		new[i,j] <- M[i,p] + M[p,j] - M[i,j] - M[p,p]
	    }
	}
	return(is.positive.semi.definite(new))
}


check_flex <- function(covparms,d, delta_A){

   ncomp <- ncomp_from_nparms(length(covparms))
   neach <- ncomp*(ncomp+1)/2
   var_inds <- 1:neach
   ran_inds <- (neach + 1):(2*neach)
   smo_inds <- (2*neach + 1):(3*neach)  	
   Z <- matrix(NA,ncomp,ncomp)	

   Z[upper.tri(Z, diag=T)] <- covparms[var_inds]
   var <- t(Z)
   var <- makeSymm(var)

   Z[upper.tri(Z, diag=T)] <- covparms[smo_inds]
   smo <- t(Z)
   smo <- makeSymm(smo)

   Z[upper.tri(Z, diag=T)] <- covparms[ran_inds]
   ran <- t(Z)
   ran <- makeSymm(ran)
   rann2 <- ran^(-2) 

   con <- matrix(0,ncomp, ncomp)
   for(i in 1:ncomp){
	for(j in 1:ncomp){
		con[i,j] <- ran[i,j]^(-delta_A - smo[i,i] - smo[j,j])*
			    gamma(smo[i,j] + d/2)/
			    gamma((smo[i,i] + smo[j,j])/2 + d/2)*gamma(smo[i,j])
	}
   }

  b1 <- is.positive.semi.definite(round(con*var,9))
  b2 <-  is.cnd(round(rann2,9))

  return(b1 && b2)
}


check_mflex <- function(covparms, beta){

   ncomp <- ncomp_from_nparms(length(covparms))
   neach <- ncomp*(ncomp+1)/2
   var_inds <- 1:neach
   ran_inds <- (neach + 1):(2*neach)
   smo_inds <- (2*neach + 1):(3*neach)  	
   Z <- matrix(NA,ncomp,ncomp)	

   Z[upper.tri(Z, diag=T)] <- covparms[var_inds]
   var <- t(Z)
   var <- makeSymm(var)

   Z[upper.tri(Z, diag=T)] <- covparms[smo_inds]
   smo <- t(Z)
   smo <- makeSymm(smo)

   Z[upper.tri(Z, diag=T)] <- covparms[ran_inds]
   ran <- t(Z)
   ran <- makeSymm(ran)
   rann2 <- ran^(-2) 

   con <- matrix(0,ncomp, ncomp)
   for(i in 1:ncomp){
	for(j in 1:ncomp){
		con[i,j] <- ((ran[i,j]^(-2))/beta)^(smo[i,j]) * exp(-smo[i,j])/gamma(smo[i,j])
	}
   }

  b0 <- is.cnd(smo) 
  b1 <- is.positive.semi.definite(con*var)
  b2 <- is.cnd(round(rann2 -beta*smo,9))

  return(b0 && b1 && b2 )
}
