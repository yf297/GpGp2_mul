get_start_parms <- function(y, locs, covfun_name, pars = NULL, n1 = NULL){

    if(covfun_name == "matern_isotropic"){
    	start_var <- var(y)
   	start_smooth <- 0.8
    	start_nug <- 0.1
    	n <- length(y)

    	randinds <- sample(1:n, min(n,200))
    	dmat <- fields::rdist(locs[randinds,])

    	start_range <- mean( dmat )/4
    	start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }

    if(covfun_name == "matern_multi_reparam"){
	var <- c(pars[1,1], 0, pars[2,1])
	ran <- c(pars[1,2], 0, pars[2,2])    
	smo <- c(pars[1,3], 0, pars[2,3])    
	nug <- c(pars[1,4],    pars[2,4])    
	n <- length(y)
	n1 <- n/2
	start_parms <- c(var(y[1:n1], y[(n1+1): (2*n1)]),  1/sqrt( (((1/ran[1]^2) + (1/ran[3]^2))/2) + 1e-8), ((smo[1] + smo[3])/2) + 1e-8 )
	
    }	

    return(  start_parms  )
    
}

