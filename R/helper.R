
ncomp_from_nparms <- function( nparms ){
    ncomp <- NA
    for(j in 1:20){
        if( nparms == 4*j*(j+1)/2 ){
            ncomp <- j
        }
    }
    if( is.na(ncomp) ){ stop("could not determine number of components") }
    return(ncomp)
}


ncomp_from_nlogparms <- function( nlogparms ){
    ncomp <- NA
    for(j in 1:20){
        if( nlogparms == (2 + 4*(j*(j+1)/2)) ){
            ncomp <- j
        }
    }
    if( is.na(ncomp) ){ stop("could not determine number of components") }
    return(ncomp)
}



multi_parms_mat <- function( covparms ){
    ncomp <- ncomp_from_nparms( length(covparms) )

    vmat <- matrix(NA, ncomp, ncomp)
    vst <- 0
    rmat <- matrix(NA, ncomp, ncomp)
    rst <- ncomp*(ncomp+1)/2
    smat <- matrix(NA, ncomp, ncomp)
    sst <- 2*rst
    nmat <- matrix(NA, ncomp, ncomp)
    nst <- 3*rst
    cnt <- 0
    for(j1 in 1:ncomp){
        for(j2 in 1:j1){
            k <- (j1-1)*j1/2 + j2         
            vmat[j1,j2] <- covparms[vst + k]
            vmat[j2,j1] <- covparms[vst + k]
            rmat[j1,j2] <- covparms[rst + k]
            rmat[j2,j1] <- covparms[rst + k]
            smat[j1,j2] <- covparms[sst + k]
            smat[j2,j1] <- covparms[sst + k]
            nmat[j1,j2] <- covparms[nst + k]
            nmat[j2,j1] <- covparms[nst + k]
        }
    }
    return( list( variance=vmat, range=rmat, smoothness=smat, nugget = nmat ) )
}

multi_parms_vec <- function( parm_mats ){

    ncomp <- nrow( parm_mats$variance )
    nmultiparms <- 3*ncomp*(ncomp+1)/2
    inds_mat <- matrix( 1:(ncomp^2), ncomp, ncomp )
    mat_inds <- 1
    for(j in 2:ncomp){ mat_inds <- c(mat_inds, inds_mat[j,1:j]) }
    
    parm_vec <- rep(NA, nmultiparms)
    inds <- 1:(ncomp*(ncomp+1)/2)

    cur_ind <- 0
    parm_vec[ cur_ind + inds ] <- parm_mats$variance[ mat_inds ]

    cur_ind <- ncomp*(ncomp+1)/2
    parm_vec[ cur_ind + inds ] <- parm_mats$range[ mat_inds ]

    cur_ind <- 2*ncomp*(ncomp+1)/2
    parm_vec[ cur_ind + inds ] <- parm_mats$smoothness[ mat_inds ]

    return( parm_vec )
}


multi_matern_parm_index <- function( ncomp, j1, j2 ){

    # assumes j1 >= j2, fix if not
    if( j2 > j1 ){
        tmp <- j2
        j2 <- j1
        j1 <- tmp
    }

    v1 <- 0
    r1 <- ncomp*(ncomp+1)/2
    s1 <- 2*r1
    n1 <- 3*r1

    k <- (j1-1)*j1/2 + j2
    index <- list()
    index$variance   <- v1 + k
    index$range      <- r1 + k
    index$smoothness <- s1 + k
    index$nugget     <- n1 + k
    
    return(index)
}

to.list <- function(y, X, locs){
	
    ncomp <- length(unique(locs[,ncol(locs)]))  
    ly <- list()
    lX <- list()
    llocs <- list()
    for(i in 1:ncomp){
	  inds <- which(locs[, ncol(locs)] == i) 
	  ly[[i]] <- y[inds]
	  lX[[i]] <- X[inds,]
	  llocs[[i]] <- locs[inds, 1:(ncol(locs)-1)]
    }
	
    return(list(ly = ly, lX  = lX, llocs = llocs))

}


f <- function(theta){

    ncomp <- NA
    for(j in 1:20){
        if( length(theta) == j*(j-1)/2 ){
            ncomp <- j
        }
    }

    L <- matrix(0,ncomp,ncomp)
    L[upper.tri(L, diag=FALSE )] <- theta
    diag(L) <- 1
    L <- t(L)
    norms <- apply(L, 1, function(x) norm(x, type = "2"))
    L <- L / norms
    M <- L%*%t(L)
}


finv <- function(M){
	
    ncomp <- nrow(M)
    tol <- 1e-8
    eigenvalues <- eigen(M, only.values = TRUE)$values

    for ( i in 1: ncomp ) {
        if ( abs( eigenvalues[i] ) < tol ) {
            eigenvalues[i] <- 0
        }
    }    
    if ( any( eigenvalues < 0 ) ) {
        theta <- rep(NA, ncomp*(ncomp-1)/2)
    }else{
    	L <- t(chol(M))
    	L <- L * kronecker(t(rep(1, ncomp)), rep(1, ncomp)/diag(L))
    	theta <- t(L)[upper.tri(L)]
    }
}
