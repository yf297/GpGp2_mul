get_start_parms <- function(y,X,locs,covfun_name){

    fitlm <- stats::lm(y ~ X - 1 )
    start_var <- summary(fitlm)$sigma^2
    start_smooth <- 0.8
    start_nug <- 0.1
    n <- length(y)

    randinds <- sample(1:n, min(n,200))
    dmat <- fields::rdist(locs[randinds,])

    if(covfun_name == "matern"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }
    if(covfun_name == "matern_multi"){
	comps <- to.list(y,X,locs)
	ly <- comps$ly
	lX <- comps$lX
	llocs <- comps$llocs
	ncomp <- length(ly)
	d <- ncol(locs)-1

	marginal.pars <- matrix(0, ncomp, nrow = ncomp, ncol = 4)
	for(i in 1:ncomp){
    	    start_parms <- get_start_parms(ly[[i]],lX[[i]],llocs[[i]],"matern")$start_parms
            # modify the nugget
            start_parms[4] <- start_parms[1]*start_parms[4]
	    marginal.pars[i,] <- start_parms
	}
	marginal.sig <- marginal.pars[,1]
	marginal.ran <- marginal.pars[,2]
	marginal.smo <- marginal.pars[,3]
	marginal.nug <- marginal.pars[,4]
	
	sig <- matrix(NA,ncomp,ncomp) 
	ran <- matrix(NA,ncomp,ncomp) 
	smo <- matrix(NA,ncomp,ncomp)
	nug <- matrix(NA,ncomp,ncomp)    

	for(i in 2:ncomp){
	    for(j in 1:(i-1)){
		sig[i,j] <- 0
		ran[i,j] <- ((marginal.ran[i]^(-2) + marginal.ran[j]^(-2))/2)^(-1/2)  #irrelevant since sig[i,j] = 0 
		smo[i,j] <- ((marginal.smo[i] + marginal.smo[j])/2)                   #irrelevant since sig[i,j] = 0 
		nug[i,j] <- 0
	    }
	}
	diag(sig) <- marginal.sig
	diag(ran) <- marginal.ran
 	diag(smo) <- marginal.smo          
	diag(nug) <- marginal.nug
	
	sig <- t(sig)[upper.tri(sig, diag = T)]
	ran <- t(ran)[upper.tri(ran, diag = T)]    
	smo <- t(smo)[upper.tri(smo, diag = T)]     
	nug <- t(nug)[upper.tri(nug, diag = T)] 
	
	start_parms  <- c(sig, ran, smo, nug)
    }
    return( list( start_parms = start_parms ) )
}
