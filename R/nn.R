

# functions to find nearest neighbors and do the grouping operations


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



c
