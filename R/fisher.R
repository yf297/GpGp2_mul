fisher_opt <- function( likfun, vartheta, link, convtol = 1e-2, max_iter = 40 ){
     

    likobj <- likfun(vartheta)
 
    # assign loglik, grad, and info
    loglik <- likobj$loglik        
    grad <- likobj$grad
    info <- likobj$info
 	
   # add a small amount of regularization
    diag(info) <- diag(info) + 0.1*min(diag(info))
 
    cat(paste0("Iter ",0,": \n"))
    cat("pars = ",  paste0(round(link(vartheta),4)), "  \n" )
    cat(paste0("loglik = ", round(loglik,6),         "  \n"))
    cat("grad = ") 
    cat(as.character(round(grad,4)),"\n")
    cat("norm_grad = ", paste0(round(norm(grad, type = "2"), 4)), " \n" )
    cat("\n\n")
    
    for(j in 1:max_iter){
        
	likobj0 <- likobj
	
     	tol <- 1e-5
        if (kappa(info) > 1 / tol) {
            #info <- 1.0*max(likobj0$info)*diag(nrow(likobj0$info))
            # regularize
            ee <- eigen(info)
            ee_ratios <- ee$values/max(ee$values)
            ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
            ee$values <- max(ee$values)*ee_ratios
            info <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
            
            #diag(info) <- diag(info) + tol*max(diag(info))
        }

	# calculate fisher step 
	step <- solve(info, grad)
	
        # if step size large, then make it smaller
        if (mean(step^2) > 1) {
            step <- step/sqrt(mean(step^2))
        }

	newvartheta <- vartheta + step
	
	likobj <- likfun(newvartheta)

        
        # check if log likelihood did not increase
        # take smaller stepsize
        if( likobj$loglik < likobj0$loglik ){
            step <- 0.25*step
            newvartheta <- vartheta + step
            likobj <- likfun(newvartheta)
        }

        # check again, move along gradient
        if( likobj$loglik < likobj0$loglik ){
            info0 <- diag( rep(mean(diag(info)),nrow(info)) )
            step <- solve(info0,grad)
            newvartheta <- vartheta + step
            likobj <- likfun(newvartheta)
        }
            
        # check once move, take smaller step along gradient
        if( likobj$loglik < likobj0$loglik ){
            info0 <- diag( rep(max(diag(info)),nrow(info)) )
            step <- solve(info0,grad)
            newvartheta <- vartheta + step
            likobj <- likfun(newvartheta)
        }



	vartheta <- vartheta + step
        loglik <- likobj$loglik        
        grad <- likobj$grad
        info <- likobj$info



	stepgrad <- c(crossprod(step,grad))

	# print some stuff out
        cat(paste0("Iter ",j,": \n"))
        cat("pars = ",  paste0(round( link(vartheta),4)), "  \n" )
        cat(paste0("loglik = ", round(loglik,6),         "  \n"))
        cat("grad = ")
        cat(as.character(round(grad,4)),"\n")
	cat("norm_grad = ", paste0(round(norm(grad, type = "2"), 4)), " \n" )
	cat("step_grad = ", paste0(round(abs(stepgrad), 4)), " \n" )
        cat("\n")
       


	
	if( abs(stepgrad) < convtol){
            break
        }

    }


    ret <- list(
        covparms = link(vartheta), 
        loglik = loglik,
        grad = likobj$grad,
        info = likobj$info,
	norm_grad = norm(grad, type = "2")
			 
    )
    
    return(ret)
}
