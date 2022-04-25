source("R_wrapper.R")
source("fisher.R")
source("linkfuns.R")
source("penalties.R")
source("start_parms.R")
source("nn.R")
library(magic)

fit_bivariate_matern <- function(y1, y2, X1, X2, locs0 ){
	
	ncomp <- 2
	pars <- matrix(rep(0,8), nrow = ncomp)

	# NNarray is same for both components
	NNarray <- find_ordered_nn(locs0, 30)


	# get link funs, same for both components
	linkfuns <- get_linkfun("matern_isotropic")
	link <- linkfuns$link
	dlink <- linkfuns$dlink
	invlink <- linkfuns$invlink

	# penalty for component 1
	penalty <- get_penalty(y1, X1, locs0, "matern_isotropic") 
	pen <- penalty$pen
	dpen <- penalty$dpen
	ddpen <-  penalty$ddpen        


	cat("Fitting component 1\n")
	Sys.sleep(3)
	# component 1 loglikelihood to be maximized
	likfun <- function(vartheta){
                likobj <- vecchia_profbeta_loglik_grad_info_matern(link(vartheta), locs0, y1, X1, NNarray) 
                likobj$loglik <- likobj$loglik + pen(link(vartheta))
                likobj$grad <- t(likobj$grad %*%dlink(vartheta) + dpen(link(vartheta))%*%dlink(vartheta))
                likobj$info <- t(dlink(vartheta))%*%likobj$info%*%dlink(vartheta) - t(dlink(vartheta))%*%ddpen(link(vartheta))%*%dlink(vartheta) 
                return(likobj)
          }


	# component 1 starting parameter
	theta <- get_start_parms(y1, X1, locs0, "matern_isotropic")
	vartheta <- invlink(theta)

	#fit component 1
	f1 <- fisher_opt(likfun, vartheta, link)
	#store parameters
	pars[1,] <- as.vector(f1$covparms)

	# penalty for component 2
	penalty <- get_penalty(y2, X2, locs0, "matern_isotropic") 
	pen <- penalty$pen
	dpen <- penalty$dpen
	ddpen <-  penalty$ddpen        

	cat("Fitting component 2\n")
	Sys.sleep(3)

	# component 2 loglikelihood to be maximized         
	likfun <- function(vartheta){
                likobj <- vecchia_profbeta_loglik_grad_info_matern(link(vartheta), locs0, y2, X2, NNarray) 
			#vecchia_profbeta_loglik_grad_info(link(vartheta),"matern_isotropic", y,X,locs0,NNarray)
                likobj$loglik <- likobj$loglik + pen(link(vartheta))
                likobj$grad <- t(likobj$grad %*%dlink(vartheta) + dpen(link(vartheta))%*%dlink(vartheta))
                likobj$info <- t(dlink(vartheta))%*%likobj$info%*%dlink(vartheta) - t(dlink(vartheta))%*%ddpen(link(vartheta))%*%dlink(vartheta) 
                return(likobj)
            }


	# component 2 starting parameter       
	theta <- get_start_parms(y2, X2, locs0, "matern_isotropic")
	vartheta <- invlink(theta)

	#fit component 2
	f2 <- fisher_opt(likfun, vartheta, link)
	#store parameters
	pars[2,] <- as.vector(f2$covparms)

	var <- c(pars[1,1], 0, pars[2,1])
	ran <- c(pars[1,2], 0, pars[2,2])    
	smo <- c(pars[1,3], 0, pars[2,3])    
	nug <- c(pars[1,4],    pars[2,4])    


	# Now to fit cross components

	# stack components
	y <- as.vector(y_mat)
	dim <- ncol(locs0)
	locs <- cbind(locs0, 1)
	colnames(locs)[dim+1] <- "comp"		
	for(i in 2:ncomp){
		locs_temp <- cbind(locs0,i)
		colnames(locs_temp)[dim+1] <- "comp"   
		locs <- rbind(locs,locs_temp)
	}
	
	
	X <- adiag(X1, X2)	

	# new NNarray
	NNarray <- find_ordered_nn(locs, 30)

	linkfuns <- get_linkfun("matern_multi_reparam", pars, dim)
	link <- linkfuns$link
	dlink <- linkfuns$dlink
	invlink <- linkfuns$invlink

	penalty <- get_penalty(y, X, locs, "matern_isotropic_reparam") 
	pen <- penalty$pen
	dpen <- penalty$dpen
	ddpen <-  penalty$ddpen        

	cat("Fitting cross components\n")
	Sys.sleep(3)

	likfun <- function(vartheta){
		theta <- link(vartheta)
		covparms <- c(var[1], theta[1], var[3], ran[1], theta[2], ran[3], smo[1], theta[3], smo[3], nug[1], nug[2]) 
		likobj <- vecchia_profbeta_loglik_grad_info_bivariate_matern(covparms, locs, y, X, NNarray)
		loglik <- likobj$loglik
		grad <- likobj$grad
		info <- likobj$info
		grad <- grad[c(2,5,8)]
		info <- info[c(2,5,8), c(2,5,8)] 

		likobj$info <- t(dlink(vartheta)) %*% info %*% dlink(vartheta)  - ddpen(vartheta)
		likobj$grad <- t(grad%*%dlink(vartheta) + dpen(vartheta))
		likobj$loglik <- loglik + pen(vartheta)
	
		return(likobj)

	}
	

        #starting param
	theta <- get_start_parms(y, X, locs, "matern_multi_reparam", pars) 
	vartheta <- invlink(theta)
	#fit cross
	f3 <- fisher_opt(likfun, vartheta, link)
	theta <- f3$covparms



	cat("Fitting all parameters simultaneously\nusing previous fits as starting parameter\n")
	Sys.sleep(5)
	likfun <- function(covparms){
		likobj <- vecchia_profbeta_loglik_grad_info_bivariate_matern(covparms, locs, y, X, NNarray)
		return(likobj)
	}

	covparms <- c(var[1], theta[1], var[3], ran[1], theta[2], ran[3], smo[1], theta[3], smo[3], nug[1], nug[2]) 

	link <- function(x){x}
	f4 <- fisher_opt(likfun, covparms, link)
	covparms <- f4$covparms
	loglik <- f4$loglik
	betahat <- f4$betahat	
	return(list(betahat = betahat, covparms = covparms, loglik = loglik))
}	


	

