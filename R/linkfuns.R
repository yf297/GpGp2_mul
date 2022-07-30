library(rootSolve)
source("helper.R")

get_linkfun <- function(covfun_name, locs){
    
    if( covfun_name == "matern_multi" ){
        link1 <- function(logparms){

	d <- ncol(locs)-1
	nlogparms <- length(logparms)
	ncomp <- ncomp_from_nlogparms(nlogparms)
	neach <- ncomp*(ncomp+1)/2
	
	#extract logparms
	log.sig <- logparms[1:neach]
	log.ran <- logparms[(neach+1): (2*neach + 1)] 
	log.smo <- logparms[(2*neach + 2) :(3*neach + 2) ]
	log.nug <- logparms[(3*neach + 3): length(logparms)]

	Z <- matrix(NA,ncomp,ncomp)	
	Z[upper.tri(Z, diag=T)] <- log.sig
	log.sig <- t(Z)

	Delta_B <- exp(log.ran[(neach+1)])
	log.ran <- log.ran[1:neach]
	Z[upper.tri(Z, diag=T)] <- log.ran
	log.ran <- t(Z)

	Delta_A <- exp(log.smo[(neach+1)])
	log.smo <- log.smo[1:neach]
	Z[upper.tri(Z, diag=T)] <- log.smo
	log.smo <- t(Z)

	Z[upper.tri(Z, diag=T)] <- log.nug
	log.nug <- t(Z)


	V <- f(t(log.sig)[upper.tri(log.sig)])
	B <- f(exp(t(log.ran)[upper.tri(log.ran)]))
        A <- f(exp(t(log.smo)[upper.tri(log.smo)]))
	S <- f(t(log.nug)[upper.tri(log.nug)])        
	
	marginal.sig <- exp(diag(log.sig))	
	marginal.ran <- exp(diag(log.ran))   	
	marginal.smo <- exp(diag(log.smo))   
	marginal.nug <- exp(diag(log.nug))
	
	ran <- matrix(NA,ncomp,ncomp)
	smo <- matrix(NA,ncomp,ncomp)	
	diag(ran) <- marginal.ran 
	diag(smo) <- marginal.smo        

	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		ran[i,j] <- (((marginal.ran[i]^(-2) + marginal.ran[j]^(-2))/2) + Delta_B*(1-B[i,j]))^(-1/2)
		smo[i,j] <-  ((marginal.smo[i] + marginal.smo[j])/2) + Delta_A*(1-A[i,j])
	    }
	}

	u <- matrix(NA,ncomp,ncomp)       
	for(i in 1:ncomp){
	    for(j in 1:i){
		u[i,j] <- (ran[i,j]^(2*Delta_A + smo[i,i] + smo[j,j]))*gamma(smo[i,j])*gamma((smo[i,i] + smo[j,j])/2 + d/2)/gamma(smo[i,j] + d/2) 
	    }
	}

	sig <- matrix(NA,ncomp,ncomp)
	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		sig[i,j] <- (marginal.sig[i]*marginal.sig[j])^(1/2) * (u[i,i]*u[j,j])^(-1/2) * V[i,j] * u[i,j]
	    }
	}
	diag(sig) <- marginal.sig 
	
	nug <- matrix(NA,ncomp,ncomp)
	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		nug[i,j] <- (marginal.nug[i]*marginal.nug[j])^(1/2) * S[i,j]
	    }
	}
	diag(nug) <- marginal.nug 

	sig <- t(sig)[upper.tri(sig, diag = T)]
	ran <- t(ran)[upper.tri(ran, diag = T)]    
	smo <- t(smo)[upper.tri(smo, diag = T)]     
	nug <- t(nug)[upper.tri(nug, diag = T)] 

	return(c(sig, ran, smo, nug))
	}
    	
	dlink1 <- function(logparms){ rootSolve::gradient(link1, logparms, pert = 1e-6)} 
	
	link2 <- function(logparms){

	d <- ncol(locs)-1
	ncomp <- ncomp_from_nlogparms(length(logparms) - 1)
	neach <- ncomp*(ncomp+1)/2
	
	#extract logparms
	log.sig <- logparms[1:neach]
	log.ran <- logparms[(neach+1): (2*neach + 1)] 
	log.smo <- logparms[(2*neach + 2) :(3*neach + 2) ]
	log.nug <- logparms[(3*neach + 3): (length(logparms)-1)]

	Z <- matrix(NA,ncomp,ncomp)	
	Z[upper.tri(Z, diag=T)] <- log.sig
	log.sig <- t(Z)

	Delta_B <- exp(log.ran[(neach+1)])
	log.ran <- log.ran[1:neach]
	Z[upper.tri(Z, diag=T)] <- log.ran
	log.ran <- t(Z)

	Delta_A <- exp(log.smo[(neach+1)])
	log.smo <- log.smo[1:neach]
	Z[upper.tri(Z, diag=T)] <- log.smo
	log.smo <- t(Z)

	Z[upper.tri(Z, diag=T)] <- log.nug
	log.nug <- t(Z)


	V <- f(t(log.sig)[upper.tri(log.sig)])
	B <- f(exp(t(log.ran)[upper.tri(log.ran)]))
        A <- f(exp(t(log.smo)[upper.tri(log.smo)]))
	S <- f(t(log.nug)[upper.tri(log.nug)])        
	beta <- exp(logparms[length(logparms)])

	marginal.sig <- exp(diag(log.sig))	
	marginal.ran <- exp(diag(log.ran))   	
	marginal.smo <- exp(diag(log.smo))   
	marginal.nug <- exp(diag(log.nug))
	
	ran <- matrix(NA,ncomp,ncomp)
	smo <- matrix(NA,ncomp,ncomp)	
	diag(ran) <- marginal.ran 
	diag(smo) <- marginal.smo        

	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		smo[i,j] <-  ((marginal.smo[i] + marginal.smo[j])/2) + Delta_A*(1-A[i,j])
		ran[i,j] <- (((marginal.ran[i]^(-2) + marginal.ran[j]^(-2))/2) + Delta_B*(1-B[i,j]) +
			     beta*(smo[i,j] - (marginal.smo[i] + marginal.smo[j])/2))^(-1/2)

	    }
	}

	u <- matrix(NA,ncomp,ncomp)       
	for(i in 1:ncomp){
	    for(j in 1:i){
		u[i,j] <- exp(smo[i,j]) * ran[i,j]^(2*smo[i,j]) * beta^(smo[i,j]) * gamma(smo[i,j]) 
	    }
	}

	sig <- matrix(NA,ncomp,ncomp)
	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		sig[i,j] <- (marginal.sig[i]*marginal.sig[j])^(1/2) * (u[i,i]*u[j,j])^(-1/2) * V[i,j] * u[i,j]
	    }
	}
	diag(sig) <- marginal.sig 
	
	nug <- matrix(NA,ncomp,ncomp)
	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		nug[i,j] <- (marginal.nug[i]*marginal.nug[j])^(1/2) * S[i,j]
	    }
	}
	diag(nug) <- marginal.nug 

	sig <- t(sig)[upper.tri(sig, diag = T)]
	ran <- t(ran)[upper.tri(ran, diag = T)]    
	smo <- t(smo)[upper.tri(smo, diag = T)]     
	nug <- t(nug)[upper.tri(nug, diag = T)] 

	return(c(sig, ran, smo, nug))
	}
    	
	dlink2 <- function(logparms){ rootSolve::gradient(link2, logparms, pert = 1e-6)} 

	link3 <- function(logparms){

	d <- ncol(locs)-1
	ncomp <- ncomp_from_nlogparms_pars(length(logparms) )
	neach <- ncomp*(ncomp+1)/2
	
	#extract logparms
	log.sig <- logparms[1:neach]
	log.ran <- logparms[neach + 1] 
	log.smo <- logparms[(neach + 2) :(neach + 2 + (ncomp-1)) ]
	log.nug <- logparms[(neach + 2 + (ncomp-1) + 1): length(logparms)]

	Z <- matrix(NA,ncomp,ncomp)	
	Z[upper.tri(Z, diag=T)] <- log.sig
	log.sig <- t(Z)

	Z[upper.tri(Z, diag=T)] <- log.nug
	log.nug <- t(Z)

	V <- f(t(log.sig)[upper.tri(log.sig)])
	S <- f(t(log.nug)[upper.tri(log.nug)])        

	marginal.sig <- exp(diag(log.sig))	
	marginal.ran <- exp(log.ran)   	
	marginal.smo <- exp(log.smo)   
	marginal.nug <- exp(diag(log.nug))
	
	ran <- matrix(NA,ncomp,ncomp)
	smo <- matrix(NA,ncomp,ncomp)	
	diag(ran) <- marginal.ran
	diag(smo) <- marginal.smo   

	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		smo[i,j] <-  (marginal.smo[i] + marginal.smo[j])/2
		ran[i,j] <- marginal.ran

	    }
	}


	sig <- matrix(NA,ncomp,ncomp)
	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		sig[i,j] <- (marginal.sig[i]*marginal.sig[j])^(1/2)*
			    V[i,j]*(gamma(smo[i,i] + d/2)^(1/2))/(gamma(smo[i,i])^(1/2))*
			           (gamma(smo[j,j] + d/2)^(1/2))/(gamma(smo[j,j])^(1/2))* 
				    gamma(smo[i,j])/gamma(smo[i,j] + d/2)
	    }
	}
	diag(sig) <- marginal.sig 
	
	nug <- matrix(NA,ncomp,ncomp)
	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		nug[i,j] <- (marginal.nug[i]*marginal.nug[j])^(1/2) * S[i,j]
	    }
	}
	diag(nug) <- marginal.nug 

	sig <- t(sig)[upper.tri(sig, diag = T)]
	ran <- t(ran)[upper.tri(ran, diag = T)]    
	smo <- t(smo)[upper.tri(smo, diag = T)]     
	nug <- t(nug)[upper.tri(nug, diag = T)] 

	return(c(sig, ran, smo, nug))
	}
    	
	dlink3 <- function(logparms){ rootSolve::gradient(link3, logparms, pert = 1e-6)}

	link4 <- function(x){
    	
	ncomp <- ncomp_from_nparms( length(x) )
    	y <- rep(NA,length(x))
    	for(j1 in 1:ncomp){
            for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            y[ii$range] <- exp(x[ii$range])
            y[ii$smoothness] <- exp(x[ii$smoothness])
            if(j1 == j2){
                y[ii$variance] <- exp(x[ii$variance])
                y[ii$nugget] <- exp(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                y[ii$variance] <- exp(0.5*(x[ii1$variance]+x[ii2$variance]))*exp2(x[ii$variance])
                y[ii$nugget] <- exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*exp2(x[ii$nugget])
                }
         	}
 		}
    	return(y)
	}
 
	dlink4 <- function(x){
    	# figure out indices of cross covariances
    	ncomp <- ncomp_from_nparms( length(x) )
    	d <- diag(length(x))
    	for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            d[ii$range,ii$range] <- exp(x[ii$range])
            d[ii$smoothness,ii$smoothness] <- exp(x[ii$smoothness])
            if(j1 == j2){
                d[ii$variance,ii$variance] <- exp(x[ii$variance])
                d[ii$nugget,ii$nugget] <- exp(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                d[ii$variance,ii$variance] <-
                    exp(0.5*(x[ii1$variance]+x[ii2$variance]))*d_exp2(x[ii$variance])
                d[ii$variance,ii1$variance] <- 
                    0.5*exp(0.5*(x[ii1$variance]+x[ii2$variance]))*exp2(x[ii$variance])
                d[ii$variance,ii2$variance] <- 
                    0.5*exp(0.*(x[ii1$variance]+x[ii2$variance]))*exp2(x[ii$variance])
                d[ii$nugget,ii$nugget] <-
                    exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*d_exp2(x[ii$nugget])
                d[ii$nugget,ii1$nugget] <- 
                    0.5*exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*exp2(x[ii$nugget])
                d[ii$nugget,ii2$nugget] <- 
                    0.5*exp(0.5*(x[ii1$nugget]+x[ii2$nugget]))*exp2(x[ii$nugget])
            }
        }
    }
    return(d)
}
  invlink4 <- function(x){
    # figure out indices of cross covariances
    ncomp <- ncomp_from_nparms( length(x) )
    y <- rep(NA,length(x))
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            ii <- multi_matern_parm_index(ncomp, j1, j2 )
            y[ii$range] <- log(x[ii$range])
            y[ii$smoothness] <- log(x[ii$smoothness])
            if(j1 == j2){
                y[ii$variance] <- log(x[ii$variance])
                y[ii$nugget] <- log(x[ii$nugget])
            } else {
                ii1 <- multi_matern_parm_index(ncomp, j1, j1)
                ii2 <- multi_matern_parm_index(ncomp, j2, j2)
                y[ii$variance] <- inv_exp2(x[ii$variance]/sqrt(x[ii1$variance]*x[ii2$variance]))
                y[ii$nugget] <- inv_exp2(x[ii$nugget]/sqrt(x[ii1$nugget]*x[ii2$nugget]))
            }
        }
    }
    return(y)
}
 

    }
	return(list(link1 = link1, dlink1 = dlink1, link2 = link2, dlink2 = dlink2, link3 = link3, dlink3 = dlink3,link4 = link4, dlink4 = dlink4, invlink4 = invlink4 ))
}



exp2 <- function(x){ exp(x) - exp(-x) }
d_exp2 <- function(x){ exp(x) + exp(-x) }
inv_exp2 <- function(x){
    y <- rep(NA, length(x))
    y[ x >= 0] <-  log( ( x[ x>=0] + sqrt(x[ x>=0]^2+4))/2 )
    y[ x < 0 ] <- -log( (-x[ x<0 ] + sqrt(x[ x<0 ]^2+4))/2 )
    return(y)
}
