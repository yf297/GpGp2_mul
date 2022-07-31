library(GpGpm, include.only = "vecchia_profbeta_loglik_grad_info")
source("R_wrapper.R")
source("nn.R")
source("fisher.R")
source("penalties.R")
source("linkfuns.R")
source("start_parms.R")
source("helper.R")
source("check.R")
source("ord.R")
library(rcompanion)
# weather example
data <- get(load("../data/b2.RData"))
nugb <- T
# load data
y <- as.vector(data[,1])

locs <- as.matrix(data[,2:ncol(data)])
locs[,ncol(locs)] <- as.numeric( as.factor( locs[,ncol(locs)] ) )
X <- model.matrix(lm( y ~ -1 + as.factor(locs[,ncol(locs)])))    
ord <- order_completely_random(locs)

y <- y[ord]
locs <- locs[ord,,drop=FALSE]
X <- as.matrix( X[ord,,drop=FALSE] )
NNarray <- nearest_multi_any(locs, 30) 

# some info
ncomp <- length(unique(locs[,ncol(locs)]))
neach <- ncomp*(ncomp+1)/2
d <- ncol(locs) - 1
M <- matrix(0.5, ncomp, ncomp)
diag(M) <- 1

# start marginal parms and logparms
start_parms <- get_start_parms(y, X , locs, "matern_multi")$start_parms
start_logparms <- log(start_parms)
start_logparms <- append(start_logparms, 0, 2*neach)
start_logparms <- append(start_logparms, 0, 3*neach+1)


# penalty
penalty <- get_penalty(y,X,locs,"matern_multi") 
pen <- penalty$pen
dpen <- penalty$dpen
ddpen <- penalty$ddpen

# logparms indices
inds <- matrix(FALSE, ncomp, ncomp)    
inds[upper.tri(inds, diag = FALSE)] <- TRUE
inds <- t(inds)

log_cross_var_inds <- which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)   
log_cross_ran_inds <-    neach +     log_cross_var_inds
log_delta_B_ind    <-  2*neach + 1
log_cross_smo_inds <-  2*neach + 1 + log_cross_var_inds
log_delta_A_ind    <-  3*neach + 2
log_cross_nug_inds <-  3*neach + 2 + log_cross_var_inds     
log_beta_ind       <-  4*neach + 3
cross_nug_inds     <-  3*neach     + log_cross_var_inds
smo_inds <- (2*neach + 1):(3*neach)

############################################################ Independent Model 
start_logparms[log_cross_var_inds] <- 0
start_logparms[log_cross_ran_inds]  <- 1 #irrelavent
start_logparms[log_delta_B_ind] <- 1 #irrelavent
start_logparms[log_cross_smo_inds]  <- 1 #irrelavent
start_logparms[log_delta_A_ind] <- 1 #irrelavent
start_logparms[log_cross_nug_inds] <- 0 

active <- rep(TRUE, length(start_logparms))
active[c(log_cross_var_inds,
	 log_cross_ran_inds,
	 log_delta_B_ind,
	 log_cross_smo_inds,
	 log_delta_A_ind,
	 log_cross_nug_inds)] <- FALSE


linkfuns <- get_linkfun("matern_multi", locs)
link <- linkfuns$link1
dlink <- linkfuns$dlink1


likfun <- function(logparms){

    lp <- rep(NA,length(start_logparms))
    lp[active] <- logparms
    lp[!active] <- start_logparms[!active]

    likobj <- vecchia_profbeta_loglik_grad_info(link(lp),"matern_multi",y,X,locs,NNarray) 
	    #vecchia_profbeta_loglik_grad_info_matern_multi(link(lp),locs, y, X, NNarray)
    likobj$loglik <- -likobj$loglik - pen(link(lp))
    likobj$grad <-  as.vector(-likobj$grad%*%dlink(lp) - dpen(link(lp))%*%dlink(lp)) 
    likobj$info <- t(dlink(lp))%*%likobj$info%*%dlink(lp) - t(dlink(lp))%*%ddpen(link(lp))%*%dlink(lp)
    likobj$grad <- likobj$grad[active]
    likobj$info <- likobj$info[active,active]
    
    return(likobj)    
}

fit_ind <- fisher_scoring_multi( 
	    likfun = likfun,
            link = link,
	    start_logparms = start_logparms, 
	    active = active)


####################################################### Flex, deltas 0
start_logparms <- fit_ind$logparms
start_logparms[log_cross_ran_inds]  <- log(finv(M)) #irrelavent
start_logparms[log_delta_B_ind] <- log(0)
start_logparms[log_cross_smo_inds]  <-  log(finv(M)) #irrelavent      
start_logparms[log_delta_A_ind] <- log(0) 
start_logparms[log_cross_nug_inds] <- 0

active <- rep(TRUE, length(start_logparms))
active[c(log_cross_ran_inds,
	 log_delta_B_ind,
	 log_cross_smo_inds,
	 log_delta_A_ind,
	 log_cross_nug_inds)] <- FALSE
active[log_cross_nug_inds] <- nugb

fit_flex0 <- fisher_scoring_multi( 
	    likfun = likfun,
            link = link,
	    start_logparms = start_logparms, 
	    active = active)

###################################################### Flex deltas > 0
start_logparms <- fit_ind$logparms
start_logparms[log_cross_ran_inds]  <- log(finv(M)) 
start_logparms[log_delta_B_ind] <- log(0.01)
start_logparms[log_cross_smo_inds]  <-  log(finv(M))  
start_logparms[log_delta_A_ind] <- log(0.01) 
start_logparms[log_cross_nug_inds] <- 0

active <- rep(TRUE, length(start_logparms))
active[c(log_cross_nug_inds)] <- nugb
if(ncomp==2){
	active[log_cross_ran_inds] <- FALSE
	active[log_cross_smo_inds] <- FALSE  
}

fit_flex1 <- fisher_scoring_multi( 
	    likfun = likfun,
            link = link,
	    start_logparms = start_logparms, 
	    active = active)

##################################################### Flex

fit_flex <- fit_flex0
if(fit_flex1$loglik > fit_flex0$loglik){
fit_flex <- fit_flex1
}

# change link and likfun
linkfuns <- get_linkfun("matern_multi", locs)
link <- linkfuns$link4
dlink <- linkfuns$dlink4
invlink  <- linkfuns$invlink4

likfun <- function(logparms){

    lp <- rep(NA,length(start_logparms))
    lp[active] <- logparms
    lp[!active] <- start_logparms[!active]
    
    likobj <- vecchia_profbeta_loglik_grad_info(link(lp),"matern_multi",y,X,locs,NNarray)      

    likobj$loglik <- -likobj$loglik - pen(link(lp))
    likobj$grad <-  as.vector(-likobj$grad%*%dlink(lp) - dpen(link(lp))%*%dlink(lp)) 
    likobj$info <- t(dlink(lp))%*%likobj$info%*%dlink(lp) - t(dlink(lp))%*%ddpen(link(lp))%*%dlink(lp)
    likobj$grad <- likobj$grad[active]
    likobj$info <- likobj$info[active,active]
    
    return(likobj)    
}


start_logparms <- invlink(fit_flex$covparms)
active <- rep(TRUE, length(start_logparms))
active[c(cross_nug_inds)] <- nugb
active[smo_inds] <- FALSE

fit_flex_final <- fisher_scoring_multi( 
	    likfun = likfun,
            link = link,
	    start_logparms = start_logparms, 
	    active = active)

delta_A <- exp(fit_flex$logparms[log_delta_A_ind])
c1 <- check_flex(fit_flex_final$covparms,d,delta_A)
if(c1 && fit_flex_final$loglik > fit_flex$loglik){
fit_flex <- fit_flex_final
}






# change link and likfun
linkfuns <- get_linkfun("matern_multi", locs)
link <- linkfuns$link2
dlink <- linkfuns$dlink2


likfun <- function(logparms){

    lp <- rep(NA,length(start_logparms))
    lp[active] <- logparms
    lp[!active] <- start_logparms[!active]
    
    likobj <- vecchia_profbeta_loglik_grad_info(link(lp),"matern_multi",y,X,locs,NNarray)      

    likobj$loglik <- -likobj$loglik - pen(link(lp))
    likobj$grad <-  as.vector(-likobj$grad%*%dlink(lp) - dpen(link(lp))%*%dlink(lp)) 
    likobj$info <- t(dlink(lp))%*%likobj$info%*%dlink(lp) - t(dlink(lp))%*%ddpen(link(lp))%*%dlink(lp)
    likobj$grad <- likobj$grad[active]
    likobj$info <- likobj$info[active,active]
    
    return(likobj)    
}


##################################################### mFlex deltas > 0
start_logparms <- fit_ind$logparms
start_logparms <- append(start_logparms, 0, length(start_logparms))

start_logparms[log_cross_ran_inds]  <- log(finv(M)) 
start_logparms[log_delta_B_ind] <- log(5)
start_logparms[log_cross_smo_inds]  <-  log(finv(M)) 
start_logparms[log_delta_A_ind] <- log(0.01) 
start_logparms[log_cross_nug_inds] <- 0
start_logparms[log_beta_ind] <- log(5)

active <- rep(TRUE, length(start_logparms))
active[c(log_cross_nug_inds)] <- nugb
if(ncomp==2){
	active[log_cross_ran_inds] <- FALSE
	active[log_cross_smo_inds] <- FALSE  
}


fit_mflex1 <- fisher_scoring_multi( 
	    likfun = likfun,
            link = link,
	    start_logparms = start_logparms, 
	    active = active)


##################################################### mFlex

fit_mflex <- fit_mflex1

# change link and likfun
linkfuns <- get_linkfun("matern_multi", locs)
link <- linkfuns$link4
dlink <- linkfuns$dlink4
invlink  <- linkfuns$invlink4

likfun <- function(logparms){

    lp <- rep(NA,length(start_logparms))
    lp[active] <- logparms
    lp[!active] <- start_logparms[!active]
    
    likobj <- vecchia_profbeta_loglik_grad_info(link(lp),"matern_multi",y,X,locs,NNarray)      

    likobj$loglik <- -likobj$loglik - pen(link(lp))
    likobj$grad <-  as.vector(-likobj$grad%*%dlink(lp) - dpen(link(lp))%*%dlink(lp)) 
    likobj$info <- t(dlink(lp))%*%likobj$info%*%dlink(lp) - t(dlink(lp))%*%ddpen(link(lp))%*%dlink(lp)
    likobj$grad <- likobj$grad[active]
    likobj$info <- likobj$info[active,active]
    
    return(likobj)    
}


start_logparms <- invlink(fit_mflex$covparms)
active <- rep(TRUE, length(start_logparms))
active[c(cross_nug_inds)] <- nugb

fit_mflex_final <- fisher_scoring_multi( 
	    likfun = likfun,
            link = link,
	    start_logparms = start_logparms, 
	    active = active)

beta <- exp(fit_mflex$logparms[log_beta_ind])
c2 <- check_mflex(fit_mflex_final$covparms,beta)
if(c2 && fit_mflex_final$loglik > fit_mflex$loglik){
fit_mflex <- fit_mflex_final
}


# change link and likfun
linkfuns <- get_linkfun("matern_multi", locs)
link <- linkfuns$link3
dlink <- linkfuns$dlink3


likfun <- function(logparms){

    lp <- rep(NA,length(start_logparms))
    lp[active] <- logparms
    lp[!active] <- start_logparms[!active]
    
    likobj <- vecchia_profbeta_loglik_grad_info(link(lp),"matern_multi",y,X,locs,NNarray)      
	    #vecchia_profbeta_loglik_grad_info_matern_multi(link(lp),locs, y, X, NNarray)

    likobj$loglik <- -likobj$loglik - pen(link(lp))
    likobj$grad <-  as.vector(-likobj$grad%*%dlink(lp) - dpen(link(lp))%*%dlink(lp)) 
    likobj$info <- t(dlink(lp))%*%likobj$info%*%dlink(lp) - t(dlink(lp))%*%ddpen(link(lp))%*%dlink(lp)
    likobj$grad <- likobj$grad[active]
    likobj$info <- likobj$info[active,active]
    
    return(likobj)    
}



####################################################### pars
start_logparms <- fit_ind$logparms
inds <- matrix(FALSE, ncomp, ncomp)    
diag(inds) <- TRUE
log_marginal_smo_inds <- (2*neach + 1) +  which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)
marginal_ran_inds <- neach + which(t(inds)[upper.tri(inds, diag = TRUE)] == TRUE)  
log_cross_nug_inds <- (neach + ncomp + 1) + log_cross_var_inds 
a <- log(mean( fit_ind$covparms[marginal_ran_inds]))

start_logparms <- c(start_logparms[1:(neach)],
		    a,
		   start_logparms[log_marginal_smo_inds],
		   start_logparms[(3*neach + 3): length(start_logparms)])

active <- rep(TRUE, length(start_logparms))
active[log_cross_nug_inds] <- nugb


fit_pars <- fisher_scoring_multi( 
	    likfun = likfun,
            link = link,
	    start_logparms = start_logparms, 
	    active = active)



parms <- rbind(fit_ind$covparms,fit_pars$covparms, fit_flex$covparms, fit_mflex$covparms)
lls <- c(fit_ind$loglik,fit_pars$loglik, fit_flex$loglik, fit_mflex$loglik)
dat <- round(cbind(parms, lls),2)
write.csv(dat,file =  "we.csv")
