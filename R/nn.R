find_ordered_nn_brute <- function( locs, m ){
     # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
     # by convention, this includes locs[j,], which is distance 0
     n <- dim(locs)[1]
     m <- min(m,n-1)
     NNarray <- matrix(NA,n,m+1)
     for(j in 1:n ){
         distvec <- c(fields::rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
         NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
     }
     return(NNarray)
}



find_ordered_nn <- function(locs,m, lonlat = FALSE, st_scale = NULL){
    
    # if locs is a vector, convert to matrix
    if( is.null(ncol(locs)) ){
        locs <- as.matrix(locs)
    }

    # number of locations
    n <- nrow(locs)
    m <- min(m,n-1)
    mult <- 2
    
    # FNN::get.knnx has strange behavior for exact matches
    # so add a small amount of noise to each location
    ee <- min(apply( locs, 2, stats::sd ))
    locs <- locs + matrix( ee*1e-4*stats::rnorm(n*ncol(locs)), n, ncol(locs) )    
    
    if(lonlat){ # convert lonlattime to xyztime or lonlat to xyz
        lon <- locs[,1]
        lat <- locs[,2]
        lonrad <- lon*2*pi/360
        latrad <- (lat+90)*2*pi/360
        x <- sin(latrad)*cos(lonrad)
        y <- sin(latrad)*sin(lonrad)
        z <- cos(latrad)
        if(ncol(locs)==3){
            time <- locs[,3]
            locs <- cbind(x,y,z,time)
        } else {
            locs <- cbind(x,y,z)
        }
    }

    
    if( !is.null(st_scale) ){ 
        d <- ncol(los)-1
        locs[ , 1:d] <- locs[ , 1:d]/st_scale[1]
        locs[ , d+1] <- locs[ , d+1]/st_scale[2]
    }

    # to store the nearest neighbor indices
    NNarray <- matrix(NA,n,m+1)

    # to the first mult*m+1 by brutce force
    maxval <- min( mult*m + 1, n )
    NNarray[1:maxval,] <- find_ordered_nn_brute(locs[1:maxval,,drop=FALSE],m)

    query_inds <- min( maxval+1, n):n
    data_inds <- 1:n

    msearch <- m

    while( length(query_inds) > 0 ){
        msearch <- min( max(query_inds), 2*msearch )
        data_inds <- 1:min( max(query_inds), n )
        NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
        less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
        sum_less_than_k <- apply(less_than_k,1,sum)
        ind_less_than_k <- which(sum_less_than_k >= m+1)

        NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))

        NNarray[ query_inds[ind_less_than_k], ] <- NN_m

        query_inds <- query_inds[-ind_less_than_k]

    }

    return(NNarray)
}


# get nearest neighbors regardless of component
nearest_multi_any <- function( locs, m ){
    locs_for_NN <- locs[, -ncol(locs), drop=FALSE ]
    NNarray <- find_ordered_nn(locs_for_NN, m=m, lonlat = FALSE, st_scale = NULL)
    return(NNarray)
}

nearest_multi_pref_this <- function(locs, m ){
    return( nearest_multi_pref(locs, m, sum_to_m_pref ) )
}

nearest_multi_balanced <- function(locs, m ){
    return( nearest_multi_pref(locs, m, sum_to_m_balanced) )
}

# pass in a function to determine how to allocate the m neighbors to the components
# see functions below for examples
nearest_multi_pref <- function( locs, m, alloc_fun ){

    comp <- as.numeric(as.factor( locs[, ncol(locs) ] ))
    ncomp <- max(locs[,ncol(locs)])
    locs_for_NN <- locs[, -ncol(locs), drop=FALSE ]
    n <- nrow(locs)

    # get fac*m nearest neighbors
    fac <- 5
    bigNN <- find_ordered_nn( locs_for_NN, m=5*m, lonlat=FALSE, st_scale=NULL )

    # initialize with values in bigNN. This sets it up properly
    NNarray <- bigNN[ , 1:(m+1) ]
    
    # loop over the rows, attempt to get a balanced number from each component
    for(j in (m+2):n){

        # get the component numbers of the neighbors
        comp_nn <- comp[ bigNN[j,2:ncol(bigNN)] ] 
        # count the number of each component
        table_comp <- rep(0, ncomp)
        for(k in 1:ncomp){
            table_comp[k] <- sum( comp_nn == k, na.rm = TRUE )
        }

        # figure out how many to take from each component
        # this is a little complicated because we might not have enough from each
        # the component of this observation
        num_from_comp <- alloc_fun( table_comp, m, comp[ NNarray[j,1] ] )

        # figure out which num_from_comp[k] neighbors to grab from component k
        # assumes distances to observations NNarray[j,] are sorted, nearest first
        inds_row <- c()
        for(k in 1:ncomp){
            if( num_from_comp[k] > 0 ){
                inds_row <- c(inds_row, which( comp_nn == k )[1:num_from_comp[k]])
            }
        }
        inds_row <- sort(inds_row)
        # add 1 because comp_nn constructed from columns 2 through ncol(bigNN)
        NNarray[j,2:(m+1)] <- bigNN[j, inds_row+1 ]
    }
    return(NNarray)
}


sum_to_m_balanced <- function( counts, m, pref_comp ){

    if( sum(counts) < m ){ stop("insufficient counts") }
    if( max( abs( counts - round(counts) ) ) > 1e-8 ){ stop("need integer counts") }

    ncomp <- length(counts)

    # initialize with zero
    num_from_comp <- rep(0, ncomp)

    # get the ordering for how we will loop over components
    not_pref <- (1:ncomp)[-pref_comp]
    comp_ord <- c( pref_comp, not_pref[ sample.int(length(not_pref)) ] )

    # try setting count for component l to k
    for(k in 1:max(counts)){
        for(l in comp_ord){
            # try adding one from this component
            if( sum( num_from_comp ) < m ){
                if( num_from_comp[l] < counts[l] ){
                    num_from_comp[l] <- num_from_comp[l] + 1
                }
            }
        }
        if( sum(num_from_comp) == m ){ break }
        if( sum(num_from_comp) > m ){ stop("exceeded max # of neighbors") }
    }
    return(num_from_comp)
}
    

sum_to_m_pref <- function( counts, m, pref_comp ){

    if( sum(counts) < m ){ stop("insufficient counts") }
    if( max( abs( counts - round(counts) ) ) > 1e-8 ){ stop("need integer counts") }

    ncomp <- length(counts)

    # initialize with zero
    num_from_comp <- rep(0, ncomp)

    # get the ordering for how we will loop over components
    not_pref <- (1:ncomp)[-pref_comp]
    comp_ord <- c( rep(pref_comp,ncomp), not_pref[ sample.int(length(not_pref)) ] )

    # try setting count for component l to k
    for(k in 1:max(counts)){
        for(l in comp_ord){
            # try adding one from this component
            if( sum( num_from_comp ) < m ){
                if( num_from_comp[l] < counts[l] ){
                    num_from_comp[l] <- num_from_comp[l] + 1
                }
            }
        }
        if( sum(num_from_comp) == m ){ break }
        if( sum(num_from_comp) > m ){ stop("exceeded max # of neighbors") }
    }
    return(num_from_comp)
}
    
