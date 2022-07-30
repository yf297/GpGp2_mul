order_completely_random <- function(locs){
    set.seed(1)
    return( sample( 1:nrow(locs) ) )
}

order_bycomponent_random <- function(locs){
    set.seed(1)
    # order by component, but then randomly within component

    # get indices of each component
    unique_comps <- unique( locs[, ncol(locs)] )
    ncomp <- length( unique_comps )
    inds_of_comp <- list()
    for(j in 1:ncomp){
        inds_of_comp[[ unique_comps[j] ]] <- which( locs[, ncol(locs)] == unique_comps[j] )
    }

    # reorder within component
    for(j in 1:ncomp){
        inds_of_comp[[ unique_comps[j] ]] <- sample( inds_of_comp[[ unique_comps[j] ]] )
    }

    # unlist and return
    return( unlist( inds_of_comp ) )

}

order_cycle_component_random <- function(locs){
    set.seed(1)
    # cycle through components, each time picking a random next observation
    
    # get indices of each component
    unique_comps <- unique( locs[, ncol(locs)] )
    ncomp <- length( unique_comps )
    inds_of_comp <- list()
    for(j in 1:nomp){
        inds_of_comp[[ unique_comps[j] ]] <- which( locs[, ncol(locs)] == unique_comps[j] )
    }

    # reorder within component
    for(j in 1:ncomp){
        inds_of_comp[[ unique_comps[j] ]] <- sample( inds_of_comp[[ unique_comps[j] ]] )
    }

    # put in a matrix
    comp_lengths <- sapply( inds_of_comp, length )
    ord_mat <- matrix(NA, max( comp_lengths ), ncomp )
    for(j in 1:ncomp){
        ord_mat[ 1:comp_lengths[j], j ] <- inds_of_comp[[ j ]]
    }

    # take the transpose and coerce to vector (essentially, do row-major ordering)
    ord_with_NA <- c( t(ord_mat) )

    ord <- ord_with_NA[ !is.na(ord_with_NA) ]
    return( ord )

}
c
