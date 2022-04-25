#include <stdlib.h>
#include "covfun.h"
#include "linalg.h"
#include "copy_in.h"
#include <string.h>
#include <omp.h>

#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })



extern "C" {
void vecchia_profbeta_likelihood_matern(
			double*  ll,
			double*  betahat,
                        double*  grad,
			double*  info,
                        double*  covparms,
                        int*     nparms_t,
			double*  y,
			int*     n_t,
                        double*  locs,
                        int*     dim_t,
			double*  X,
			int*     p_t,
                        int*     NNarray,
                        int*     m_t){


    int nparms = *nparms_t;
    int n   = *n_t;
    int dim = *dim_t;
    int p   = *p_t;
    int m  = *m_t;

    double* space1    = (double*) malloc((p*p + p + nparms + (p*p*nparms) + (p*nparms) + nparms)*sizeof(double));
    double* XSX       = space1;
    double* ySX       = XSX  + p*p; 
    double* dySy      = ySX  + p;
    double* dXSX      = dySy + nparms;
    double* dySX      = dXSX + p*p*nparms;
    double* dlogdet   = dySX + p*nparms;

    double ySy;
    double logdet;

    memset(space1, 0, (p*p + p + nparms + (p*p*nparms) + (p*nparms) + nparms)*sizeof(double));
    ySy  = 0.0;
    logdet = 0.0;


#pragma omp parallel
{   
    double* space3     = (double*) malloc( ((m*dim) + m + (m*p) + m + (m*m) + (nparms*m*m) + m + (m*nparms) + p + p ) * sizeof(double));
    double* locsub     =  space3; 
    double* ysub       =  locsub  + m*dim;
    double* Xsub       =  ysub    + m;
    double* choli2     =  Xsub    + m*p;
    double* covmat     =  choli2  + m;
    double* dcovmat    =  covmat  + m*m;
    double* LidSLi3    =  dcovmat + nparms*m*m; 
    double* LidSLi2    =  LidSLi3 + m;
    double* v1         =  LidSLi2 + m*nparms ; 
    double* LiX0i2     =  v1 + p;

#pragma omp for reduction(+:ySy) reduction(+:XSX[:p*p]) reduction(+:ySX[:p]) reduction(+:logdet) \
                reduction(+:dySy[:nparms]) reduction(+:dXSX[:nparms*p*p]) reduction(+:dySX[:nparms*p])\
                reduction(+:dlogdet[:nparms]) 
    
    for(int i = 0; i < n; ++i){
      
        int bsize = MIN(i+1, m);
        
        memset(choli2, 0.0, bsize*sizeof(double));      
        choli2[bsize-1] = 1.0;

        // Compute covariance and cholesky
        copy_in(ysub, locsub, Xsub, bsize, y, locs, X, n, dim, p, NNarray, i*m);
        matern_isotropic(covmat, dcovmat, bsize, covparms, locsub, dim);
        cholesky(covmat, bsize);          
        // Extract last row of inverse cholesky
        transpose_triangular_solve(covmat, choli2, bsize, 1);
     
        // ysub becomes Li*y
        // Xsub becomes Li*X
       triangular_solve(covmat, ysub, bsize, 1);
       triangular_solve(covmat, Xsub, bsize, p);
        
        // Get last row of LiX0
        memset(LiX0i2, 0.0, p*sizeof(double));
        for(int i = 0; i <p; i++){
            LiX0i2[i] = Xsub[i*bsize + bsize-1];
        }

        double Liy0i2 = ysub[(bsize-1)] ;

        outer_product(1, LiX0i2, LiX0i2, XSX, p, p);
        
	vector_add(LiX0i2, Liy0i2, ySX, p);
        
        logdet += 2*log(covmat[(bsize-1)*(bsize) + bsize-1]);
        ySy += pow( Liy0i2, 2);

        if(bsize > 1){
            for(int j = 0; j < nparms; j++){
                // Compute LidLi3 L*dS*Li
                symmetric_matrix_vector(dcovmat + j*bsize*bsize, choli2, bsize, LidSLi3);
                triangular_solve(covmat, LidSLi3, bsize, 1);
                double LidSLi3i2 = LidSLi3[bsize-1];
                
                // Compute s1 = <LidSLi3, ysub>
                double s1 = dot(LidSLi3, ysub, bsize);
                // Compute vec1 = (LiX0)'LidDLi3
                transpose_matrix_vector(Xsub, LidSLi3, bsize, p, v1);
                

                dySy[j] += 2*s1*Liy0i2 - LidSLi3i2*Liy0i2*Liy0i2;
                outer_product(1, v1, LiX0i2, dXSX + j*p*p, p, p);
                outer_product(1, LiX0i2, v1, dXSX + j*p*p, p, p);
                outer_product(-1*LidSLi3i2, LiX0i2, LiX0i2, dXSX + j*p*p, p, p);
                scale(v1, Liy0i2, p);
                vector_add(LiX0i2, s1, v1, p);
                vector_add(LiX0i2, -1*LidSLi3i2*Liy0i2, v1, p);
                vector_add(v1, 1, dySX + j*p, p);
                dlogdet[j] += LidSLi3i2;
		copy(LidSLi3, LidSLi2 + j*bsize, bsize);

	    }   
		
	    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    info[j*nparms + i] += (dot(LidSLi2 + j*bsize, LidSLi2 + i*bsize, bsize) - 0.5*LidSLi2[j*bsize + bsize-1]*LidSLi2[i*bsize + bsize-1]);
		    info[i*nparms + j] = info[j*nparms + i];}}

        }else{
            for(int j = 0; j < nparms; j++){
                double LidSLi = *(dcovmat + j*bsize*bsize) * (1/pow(covmat[0],2));
                cblas_dger(CblasColMajor, p, p , LidSLi, Xsub, 1, Xsub, 1, dXSX + j*p*p, p); 
                dySy[j] += Liy0i2 * LidSLi * Liy0i2;
                vector_add(Xsub, Liy0i2 * LidSLi, dySX + j*p, p);
                dlogdet[j] += LidSLi;
               // LidSLi2[j] = LidSLi;
            }

	    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                info[j*nparms + i] += 0.5*LidSLi2[i]*LidSLi2[j];
		info[i*nparms + j] = info[j*nparms + i];}}
            
            }

    }

    free(space3);
}    
    //compute betahat
    copy(ySX, betahat, p);
    symmetric_solve(XSX, betahat, p,1);
    
    // compute log-likelihood
    double temp = dot(ySX, betahat, p);
    double sig2 = ySy - temp;
    *ll = -0.5* (n*log(2.0 * M_PI) + logdet + sig2);
    for(int j=0; j<nparms; j++){
        grad[j] = 0.0;
        grad[j] -= 0.5*dlogdet[j];
        grad[j] += 0.5*dySy[j];
        grad[j] -= 1.0*dot(betahat, dySX + j*p, p );
        symmetric_matrix_vector(dXSX + j*p*p, betahat, p, ySX);
        grad[j] += 0.5*dot(betahat, ySX, p );
    }

    free(space1);
}
}



extern "C" {
void vecchia_profbeta_likelihood_bivariate_matern(
			double*  ll,
			double*  betahat,
                        double*  grad,
			double*  info,
                        double*  covparms,
                        int*     nparms_t,
			double*  y,
			int*     n_t,
                        double*  locs,
                        int*     dim_t,
			double*  X,
			int*     p_t,
                        int*     NNarray,
                        int*     m_t){


    int nparms = *nparms_t;
    int n   = *n_t;
    int dim = *dim_t;
    int p   = *p_t;
    int m  = *m_t;

    double* space1    = (double*) malloc((p*p + p + nparms + (p*p*nparms) + (p*nparms) + nparms)*sizeof(double));
    double* XSX       = space1;
    double* ySX       = XSX  + p*p; 
    double* dySy      = ySX  + p;
    double* dXSX      = dySy + nparms;
    double* dySX      = dXSX + p*p*nparms;
    double* dlogdet   = dySX + p*nparms;

    double ySy;
    double logdet;

    memset(space1, 0, (p*p + p + nparms + (p*p*nparms) + (p*nparms) + nparms)*sizeof(double));
    ySy  = 0.0;
    logdet = 0.0;


#pragma omp parallel
{   
    double* space3     = (double*) malloc( ((m*dim) + m + (m*p) + m + (m*m) + (nparms*m*m) + m + (m*nparms) + p + p ) * sizeof(double));
    double* locsub     =  space3; 
    double* ysub       =  locsub  + m*dim;
    double* Xsub       =  ysub    + m;
    double* choli2     =  Xsub    + m*p;
    double* covmat     =  choli2  + m;
    double* dcovmat    =  covmat  + m*m;
    double* LidSLi3    =  dcovmat + nparms*m*m; 
    double* LidSLi2    =  LidSLi3 + m;
    double* v1         =  LidSLi2 + m*nparms ; 
    double* LiX0i2     =  v1 + p;

#pragma omp for reduction(+:ySy) reduction(+:XSX[:p*p]) reduction(+:ySX[:p]) reduction(+:logdet) \
                reduction(+:dySy[:nparms]) reduction(+:dXSX[:nparms*p*p]) reduction(+:dySX[:nparms*p])\
                reduction(+:dlogdet[:nparms]) 
    
    for(int i = 0; i < n; ++i){
      
        int bsize = MIN(i+1, m);
        
        memset(choli2, 0.0, bsize*sizeof(double));      
        choli2[bsize-1] = 1.0;

        // Compute covariance and cholesky
        copy_in(ysub, locsub, Xsub, bsize, y, locs, X, n, dim, p, NNarray, i*m);
        matern_isotropic_bivariate(covmat, dcovmat, bsize, covparms, locsub, dim);
        cholesky(covmat, bsize);          
        // Extract last row of inverse cholesky
        transpose_triangular_solve(covmat, choli2, bsize, 1);
     
        // ysub becomes Li*y
        // Xsub becomes Li*X
       triangular_solve(covmat, ysub, bsize, 1);
       triangular_solve(covmat, Xsub, bsize, p);
        
        // Get last row of LiX0
        memset(LiX0i2, 0.0, p*sizeof(double));
        for(int i = 0; i <p; i++){
            LiX0i2[i] = Xsub[i*bsize + bsize-1];
        }

        double Liy0i2 = ysub[(bsize-1)] ;

        outer_product(1, LiX0i2, LiX0i2, XSX, p, p);
        
	vector_add(LiX0i2, Liy0i2, ySX, p);
        
        logdet += 2*log(covmat[(bsize-1)*(bsize) + bsize-1]);
        ySy += pow( Liy0i2, 2);

        if(bsize > 1){
            for(int j = 0; j < nparms; j++){
                // Compute LidLi3 L*dS*Li
                symmetric_matrix_vector(dcovmat + j*bsize*bsize, choli2, bsize, LidSLi3);
                triangular_solve(covmat, LidSLi3, bsize, 1);
                double LidSLi3i2 = LidSLi3[bsize-1];
                
                // Compute s1 = <LidSLi3, ysub>
                double s1 = dot(LidSLi3, ysub, bsize);
                // Compute vec1 = (LiX0)'LidDLi3
                transpose_matrix_vector(Xsub, LidSLi3, bsize, p, v1);
                

                dySy[j] += 2*s1*Liy0i2 - LidSLi3i2*Liy0i2*Liy0i2;
                outer_product(1, v1, LiX0i2, dXSX + j*p*p, p, p);
                outer_product(1, LiX0i2, v1, dXSX + j*p*p, p, p);
                outer_product(-1*LidSLi3i2, LiX0i2, LiX0i2, dXSX + j*p*p, p, p);
                scale(v1, Liy0i2, p);
                vector_add(LiX0i2, s1, v1, p);
                vector_add(LiX0i2, -1*LidSLi3i2*Liy0i2, v1, p);
                vector_add(v1, 1, dySX + j*p, p);
                dlogdet[j] += LidSLi3i2;
		copy(LidSLi3, LidSLi2 + j*bsize, bsize);

	    }   
		
	    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    info[j*nparms + i] += (dot(LidSLi2 + j*bsize, LidSLi2 + i*bsize, bsize) - 0.5*LidSLi2[j*bsize + bsize-1]*LidSLi2[i*bsize + bsize-1]);
		    info[i*nparms + j] = info[j*nparms + i];}}

        }else{
            for(int j = 0; j < nparms; j++){
                double LidSLi = *(dcovmat + j*bsize*bsize) * (1/pow(covmat[0],2));
                cblas_dger(CblasColMajor, p, p , LidSLi, Xsub, 1, Xsub, 1, dXSX + j*p*p, p); 
                dySy[j] += Liy0i2 * LidSLi * Liy0i2;
                vector_add(Xsub, Liy0i2 * LidSLi, dySX + j*p, p);
                dlogdet[j] += LidSLi;
               // LidSLi2[j] = LidSLi;
            }

	    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                info[j*nparms + i] += 0.5*LidSLi2[i]*LidSLi2[j];
		info[i*nparms + j] = info[j*nparms + i];}}
            
            }

    }

    free(space3);
}    
    //compute betahat
    copy(ySX, betahat, p);
    symmetric_solve(XSX, betahat, p,1);
    
    // compute log-likelihood
    double temp = dot(ySX, betahat, p);
    double sig2 = ySy - temp;
    *ll = -0.5* (n*log(2.0 * M_PI) + logdet + sig2);
    for(int j=0; j<nparms; j++){
        grad[j] = 0.0;
        grad[j] -= 0.5*dlogdet[j];
        grad[j] += 0.5*dySy[j];
        grad[j] -= 1.0*dot(betahat, dySX + j*p, p );
        symmetric_matrix_vector(dXSX + j*p*p, betahat, p, ySX);
        grad[j] += 0.5*dot(betahat, ySX, p );
    }

    free(space1);
}
}


