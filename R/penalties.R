library(rootSolve)
# penalty functions


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


    if(covfun_name == "matern_isotropic"){
          pen <- function(x){
              pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1) + pen_sm_hi(x,3)
          }
         dpen <- function(x){
             dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1) + dpen_sm_hi(x,3)
         }
        ddpen <- function(x){
            ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) + ddpen_sm_hi(x,3) }
    }
    

    if(covfun_name == "matern_multi_reparam"){
    	pen <- function(vartheta){
		eta <- h(vartheta)
		delta_B <- eta[2]
		delta_A <- eta[3]
		return( -(delta_A^2  + delta_B^2 ))
	}

        dpen <- function(vartheta){gradient(pen, vartheta)}
	ddpen <- function(vartheta){hessian(pen, vartheta)}

    }

    return( list( pen = pen, dpen = dpen, ddpen = ddpen ) )

}

