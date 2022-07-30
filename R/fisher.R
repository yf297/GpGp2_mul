
test_likelihood_object <- function(likobj){
    
    pass <- TRUE
    allvals <- c( likobj$loglik, likobj$grad, c(likobj$info) )
    if( sum(is.na(allvals)) > 0  ||  sum( abs(allvals) == Inf ) > 0 ){
        pass <- FALSE
    }
    return(pass)
}


condition_number <- function(info){
    # assumes that information matrix has finite numbers in it
    ee <- eigen(info)
    return( max(ee$values)/min(ee$values) )
}    




check_step <- function(step,silent){

    # if step size large, then make it smaller
    if (sum(step^2) > 1) {
        if(!silent) cat("step too large\n")
        step <- step/sqrt(sum(step^2))
    }
    return(step)

}

check_likelihood <- function(likobj, likobj0){

    pass <- TRUE
    allvals <- c( likobj$loglik, likobj$grad, c(likobj$info) )
    if( sum(is.na(allvals)) > 0  ||  sum( abs(allvals) == Inf ) > 0 ){
        pass <- FALSE
    }
    tol <- 2
    if( likobj$loglik > likobj0$loglik + tol ){
        pass <- FALSE 
    }
    return(pass)

}


fisher_scoring_multi <- function( likfun, link, 
    silent = FALSE, convtol = 1e-4, max_iter = 40, start_logparms, active ){
    
    # evaluate function at initial values
    logparms <- start_logparms[active]
    likobj <- likfun(logparms)
    
    # test likelihood object    
    if( !test_likelihood_object(likobj) ){
        stop("Invalid Starting values")
        logparms <- 0.1*logparms
        likobj <- likfun(logparms)
    }
    
    # assign loglik, grad, and info
    loglik <- likobj$loglik        
    grad <- likobj$grad
    info <- as.matrix(likobj$info)
    # add a small amount of regularization
    info <- info + 0.1*min(diag(info))*diag(nrow(info))

    # print some stuff out
    lp <- rep(NA,length(start_logparms))
    lp[active] <- logparms
    lp[!active] <- start_logparms[!active]
    
    if(!silent){
        cat(paste0("Iter ",0,": \n"))
        cat("pars = ",  paste0(round(link(lp),4)), "  \n" )
        cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
       cat("grad = ")
        cat(as.character(round(-grad,3)))
        cat("\n\n")
    }
    
    iter <- 0
    for(j in seq_len(max_iter)){

        iter <- j
        
        likobj0 <- likobj
        
        # if condition number of info matrix large, regularize
        tol <- 1e-10
        if ( 1/condition_number(info) < tol) {
            if (!silent) cat("Cond # of info matrix > 1/tol \n")
            # regularize
            ee <- eigen(info)
            ee_ratios <- ee$values/max(ee$values)
            ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
            ee$values <- max(ee$values)*ee_ratios
            info <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
        }

        # calculate fisher step 
        step <- -solve(info, grad)
        step <- check_step(step,silent)

        # take step and calculate loglik, grad, and info
        newlogparms <- logparms + step
        likobj <- likfun(newlogparms)
        
        # take smaller step along fisher direction if bad step
        if( !check_likelihood(likobj, likobj0) ){
            if (!silent) cat("increase or inf or na or nan in likobj\n")
            step <- 0.25 * step
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        # take step along gradient if bad step
        if( !check_likelihood(likobj, likobj0) ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- -0.25*solve(info0,grad)
            step <- check_step(step,silent)
            if(!silent){print("moving along gradient")}              
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        # take step along gradient if bad step
        if( !check_likelihood(likobj, likobj0) ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- -0.05*solve(info0,grad)
            step <- check_step(step,silent)
            if(!silent){print("moving along gradient")}              
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        # take step along gradient if bad step
        if( !check_likelihood(likobj, likobj0) ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- -0.01*solve(info0,grad)
            step <- check_step(step,silent)
            if(!silent){print("moving along gradient")}              
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
        }

        no_decrease <- FALSE
        # set no_decrease to true if we can't move
        if( !check_likelihood(likobj, likobj0) ){
            if(!silent){print("could not increase along gradient")}              
            no_decrease <- TRUE
        }

        stepgrad <- c(crossprod(step,grad))
        
        # redefine logparms, loglik, grad, info
        logparms <- logparms + step
        loglik <- likobj$loglik        
        grad <- likobj$grad
        info <- likobj$info

        lp[active] <- logparms
        # print some stuff out
        if(!silent){
            cat(paste0("Iter ",j,": \n"))
            cat("pars = ",  paste0(round(link(lp),4)), "  \n" )
            cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
            cat("grad = ")
            cat(as.character(round(-grad,4)),"\n")
            cat("step dot grad = ",stepgrad,"\n")
            cat("\n")
        }
        
        # if gradient is small, then break and return the results        
        if( abs(stepgrad) < convtol || no_decrease ){
	    break
        }
    }

    # collect information to return
    betahat <- as.vector(likobj$betahat)

    ret <- list(
        covparms = link(lp), 
        logparms = lp,
        betahat = betahat, 
        loglik = -loglik,
        grad = likobj$grad,
        info = likobj$info,
	iter = iter
    )
    return(ret)
}
 
