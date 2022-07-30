library(rootSolve)
source("linkfuns.R")

expit <- function(x){ exp(x)/(1+exp(x)) }
intexpit <- function(x){ log(1+exp(x)) }
pen_hi <- function(x,tt,aa){ -tt*intexpit(x-aa) }
dpen_hi <- function(x,tt,aa){ -tt*expit(x-aa) }
ddpen_hi <- function(x,tt,aa){ -tt*expit(x-aa)/(1+exp(x-aa)) }
pen_lo <- function(x,tt,aa){ -tt*intexpit(-x+aa) }
dpen_lo <- function(x,tt,aa){ +tt*expit(-x+aa) }
ddpen_lo <- function(x,tt,aa){ -tt*expit(-x+aa)/(1+exp(-x+aa)) }

pen_loglo <- function(x,tt,aa){ 
    if(x==0){ return(0.0) 
    } else { 
        return( pen_lo(log(x),tt,aa) )
    }
}

dpen_loglo <- function(x,tt,aa){ 
    if( x==0 ){
        return(0.0) 
    } else {
        return( dpen_lo(log(x),tt,aa)/x )
    }
}


ddpen_loglo <- function(x,tt,aa){ 
    if( x==0 ){
        return( 0.0 )
    } else {
        return( ddpen_lo(log(x),tt,aa)/x^2 - dpen_lo(log(x),tt,aa)/x^2 )
    }
}



get_penalty <- function(y, X, locs,covfun_name){

    fitlm <- stats::lm(y ~ X - 1 )
    vv <- summary(fitlm)$sigma^2
    # by default, no penalty
    pen <- function(x) 0.0
    dpen <- function(x) rep(0,length(x))
    ddpen <- function(x) matrix(0,length(x),length(x))
    # nugget penalty
    pen_nug <- function(x,j){ pen_loglo(x[j],.1,log(0.001)) }
    dpen_nug <- function(x,j){
        dpen <- rep(0,length(x))
        dpen[j] <- dpen_loglo(x[j],.1,log(0.001))
        return(dpen)
    }
    ddpen_nug <- function(x,j){
        ddpen <- matrix(0,length(x),length(x))
        ddpen[j,j] <- ddpen_loglo(x[j],.1,log(0.001))
        return(ddpen)
    }
    # smoothness penalty
    pen_sm <- function(x,j){ pen_loglo(x[j],.1,log(0.2)) }
    dpen_sm <- function(x,j){
        dpen <- rep(0,length(x))
        dpen[j] <- dpen_loglo(x[j],.1,log(0.2))
        return(dpen)
    }
    ddpen_sm <- function(x,j){
        ddpen <- matrix(0,length(x),length(x))
        ddpen[j,j] <- ddpen_loglo(x[j],.1,log(0.2))
        return(ddpen)
    }
    # variance penalty
    # dangerous because vv could get redefined
    pen_var <- function(x,j){ pen_hi(x[j]/vv,1,6) }
    dpen_var <- function(x,j){
        dpen <- rep(0,length(x))
        dpen[j] <- 1/vv*dpen_hi(x[j]/vv,1,6)
        return(dpen)
    }
    ddpen_var <- function(x,j){
        ddpen <- matrix(0,length(x),length(x))
        ddpen[j,j] <- 1/vv^2*ddpen_hi(x[j]/vv,1,6)
        return(ddpen)
    }

    # penalty on large smoothness parameters
    pen_sm_hi <- function(x,j){sm <-8.0; bb<-0.5; tt<-3.0; pen_hi(x[j]/bb,tt,sm) }
    dpen_sm_hi <- function(x,j){
        sm <- 8.0
        bb <- 0.5
        tt <- 3.0
        dpen <- rep(0,length(x))
        dpen[j] <- 1/bb*dpen_hi(x[j]/bb,tt,sm)
        return(dpen)
    }
    ddpen_sm_hi <- function(x,j){
        sm <- 8.0
        bb <- 0.5
        tt <- 3.0
        ddpen <- matrix(0,length(x),length(x))
        ddpen[j,j] <- 1/bb^2*ddpen_hi(x[j]/bb,tt,sm)
        return(ddpen)
    }


    if(covfun_name == "matern"){
        pen <- function(x){
              pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1) + pen_sm_hi(x,3)
          }
        dpen <- function(x){
             dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1) + dpen_sm_hi(x,3)
         }
        ddpen <- function(x){
            ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) + ddpen_sm_hi(x,3) }
    }
    

    if(covfun_name == "matern_multi"){
	pen <- function(x){
#	    dB <- (logparms[2*neach + 1])^2
#   	    dA <- (logparms[3*neach + 2])^2 
#	    link <- get_linkfun("matern_multi", locs)$link
#	    x <- link(logparms)
	    p <- 0      
	    pmats <- multi_parms_mat( x )
    	    ncomp <- nrow(pmats$variance)
	    for(j in 1:ncomp){
                these_parms <- c(pmats[[1]][j,j], pmats[[2]][j,j], pmats[[3]][j,j], pmats[[4]][j,j])
                ii <- locs[, ncol(locs)] == j
                pen_obj = get_penalty(
                    y[ii],
                    X[ii, , drop=FALSE],
                    locs[ii, -ncol(locs), drop=FALSE ],
                    "matern"
                )
                p <- p + pen_obj$pen( these_parms )
            }
            return(p)
	    }

	dpen <- function(x){

            dp <- rep(0, length(x))
            pmats <- multi_parms_mat( x )
            ncomp <- nrow(pmats$variance)

	        for(j in 1:ncomp){
                these_parms <- c(pmats[[1]][j,j], pmats[[2]][j,j], pmats[[3]][j,j], pmats[[4]][j,j])
                ii <- locs[, ncol(locs)] == j
                pen_obj = get_penalty(
                    y[ii],
                    X[ii, , drop=FALSE],
                    locs[ii, -ncol(locs), drop=FALSE ],
                    "matern"
                )
                inds <- unlist( multi_matern_parm_index( ncomp, j, j ) )
                dp[inds] <- dp[inds] + pen_obj$dpen( these_parms )
            }
            return(dp)
        }

        ddpen <- function(x){

            ddp <- matrix(0, length(x), length(x) )
            pmats <- multi_parms_mat( x )
            ncomp <- nrow(pmats$variance)

	        for(j in 1:ncomp){
                these_parms <- c(pmats[[1]][j,j], pmats[[2]][j,j], pmats[[3]][j,j], pmats[[4]][j,j])
                ii <- locs[, ncol(locs)] == j
                pen_obj = get_penalty(
                    y[ii],
                    X[ii, , drop=FALSE],
                    locs[ii, -ncol(locs), drop=FALSE ],
                    "matern"
                )
                inds <- unlist( multi_matern_parm_index( ncomp, j, j ) )
                ddp[inds,inds] <- ddp[inds,inds] + pen_obj$ddpen( these_parms )
            }
            #ddpm <- ddpen_logdet_cross_spec(x, effrange)-ddpen_logvar_marg_spec(x,effrange)
            #ddp <- ddp+fac*ddpm
            return(ddp)
        }


#        dpen <- function(logparms){rootSolve::gradient(pen,logparms)}
#	 ddpen <- function(logparms){rootSolve::hessian(pen,logparms)}       

   }
    return( list( pen = pen, dpen = dpen, ddpen = ddpen ) )
} 
