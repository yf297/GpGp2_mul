source("R_wrapper.R")
library(rootSolve)

get_linkfun <- function(covfun_name, pars = NULL, dim = NULL){

	if(covfun_name == "matern_isotropic"){
        	link <- function(x){ c(exp(x[1]), exp(x[2]), exp(x[3]), exp(x[4])) }
        	dlink <- function(x){ gradient(link, x) } 
        	ddlink <- function(x){ hessian(link, x) }
		invlink <- function(x){ c(log(x[1]), log(x[2]), log(x[3]), log(x[4]) ) } 
        }
	
	if(covfun_name == "matern_multi_reparam"){
 		var <- c(pars[1,1], 0, pars[2,1])
		ran <- c(pars[1,2], 0, pars[2,2])    
		smo <- c(pars[1,3], 0, pars[2,3])    
		nug <- c(pars[1,4],    pars[2,4])    


		g <- function(eta){
			rho_V   <- eta[1]
			delta_B <- eta[2]
			delta_A <- eta[3]

			smo12 <- ((smo[1] + smo[3])/2) + delta_A
		        ran12 <-  1/sqrt( (((1/ran[1]^2) + (1/ran[3]^2))/2) + delta_B)  	
			covparms <- c(ran[1], delta_B, ran[3], smo[1], delta_A, smo[3])
			t <- con(dim,covparms)
			var12 <- rho_V*(var[1]*var[3]* t)^0.5
			return(c(var12, ran12, smo12))
		
		}

		ginv <- function(theta){
			var12  <- theta[1]
			ran12  <- theta[2]
			smo12  <- theta[3]

			delta_A <- smo12 -  ((smo[1] + smo[3])/2)
		        delta_B <- (1/ran12^2) - ( ((1/ran[1]^2) + (1/ran[3]^2))/2) 	
			covparms <- c(ran[1], delta_B, ran[3], smo[1], delta_A, smo[3])
			t <- con(dim,covparms)
			rho_V <- var12 /(var[1]*var[3]*t)^0.5
			return(c(rho_V, delta_B, delta_A))
		
		}


		h <- function(vartheta){
			h1  <- tanh(vartheta[1])
			h2 <- exp(vartheta[2])
			h3 <- exp(vartheta[3])
			return(c(h1,h2,h3))
		}
	
		hinv <- function(vartheta){
			h1  <- atanh(vartheta[1])
			h2 <-  log(vartheta[2])
			h3 <-  log(vartheta[3])
			return(c(h1,h2,h3))
		}


		link <- function(x){ g(h(x)) }
        	dlink <- function(x){ gradient(link, x) } 
        	ddlink <- function(x){ hessian(link, x) }
		invlink <- function(x){ hinv(ginv(x))}

	}

	return(list(
        link = link, dlink = dlink, invlink = invlink))

}




