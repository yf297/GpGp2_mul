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
void vecchia_meanzero_likelihood_matern(double*  ll,
                        double*  grad,
			double*  info,
                        double*  covparms,
                        int*     nparms_t,
			double*  y,
			int*     n_t,
                        double*  locs,
                        int*     dim_t,
                        int*     NNarray,
                        int*     m_t){


    int nparms = *nparms_t;
    
    int n = *n_t;
    int dim = *dim_t;
    int m  = *m_t;

    double* dySy = (double*) malloc( 2*nparms * sizeof(double) );
    memset(dySy, 0, 2*nparms*sizeof(double));

    double* dlogdet = dySy + nparms;

    double ySy = 0;
    double logdet = 0;


#pragma omp parallel
{   
    double* space     = (double*) malloc( ((m*dim) + m + m + (m*m) + (nparms*m*m) + m + m*nparms) * sizeof(double));
    double* locsub     =  space; 
    double* ysub       =  locsub  + m*dim;
    double* choli2     =  ysub    + m;
    double* covmat     =  choli2  + m;
    double* dcovmat    =  covmat  + m*m;
    double* LidSLi3    =  dcovmat + nparms*m*m; 
    double* LidSLi2    =  LidSLi3 + m;

#pragma omp for reduction(+:ySy) reduction(+:logdet) reduction(+:dySy[:nparms]) reduction(+:dlogdet[:nparms]) reduction(+:info[:nparms*nparms])
    
    for(int i = 0; i < n; ++i){
      
        int bsize = MIN(i+1, m);
        
        memset(choli2, 0.0, bsize*sizeof(double));      
        choli2[bsize-1] = 1.0;

        // Compute covariance and cholesky
        copy_in(ysub, locsub, bsize, y, locs, n, dim, NNarray, i*m);
        matern_isotropic(covmat, dcovmat, bsize, covparms, locsub, dim);
        cholesky(covmat, bsize);          
        
	// Extract last row of inverse cholesky
        transpose_triangular_solve(covmat, choli2, bsize, 1);
     
        // ysub becomes Li*y
        triangular_solve(covmat, ysub, bsize, 1);
        
        double Liy0i2 = ysub[(bsize-1)] ;
        
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
                
		dySy[j] += 2*s1*Liy0i2 - LidSLi3i2*Liy0i2*Liy0i2;
                dlogdet[j] += LidSLi3i2;
		copy(LidSLi3, LidSLi2 + j*bsize, bsize);


	    }
		for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    info[j*nparms + i] += (dot(LidSLi2 + j*bsize, LidSLi2 + i*bsize, bsize) - 0.5*LidSLi2[j*bsize + bsize-1]*LidSLi2[i*bsize + bsize-1]);
		    info[i*nparms + j] = info[j*nparms + i];

            }}

            
        }else{
	    for(int j = 0; j < nparms; j++){
                double LidSLi = *(dcovmat + j*bsize*bsize) * (1/pow(covmat[0],2));
                dySy[j] += Liy0i2 * LidSLi * Liy0i2;
                dlogdet[j] += LidSLi; 
		LidSLi2[j] = LidSLi;
	    }


		for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    info[j*nparms + i] += 0.5*LidSLi2[i]*LidSLi2[j];
		    info[i*nparms + j] = info[j*nparms + i];

            }}
            
            }

    }

    free(space);
}    
    
    // compute log-likelihood, gradient
    *ll = -0.5* (n*log(2.0 * M_PI) + logdet + ySy);
    for(int i=0; i<nparms; i++){
        grad[i] = -0.5*dlogdet[i] + 0.5*dySy[i];
    }

}
}

extern "C" {
void vecchia_meanzero_likelihood_bivariate_matern(double*  ll,
                        double*  grad,
			double*  info,
                        double*  covparms,
                        int*     nparms_t,
			double*  y,
			int*     n_t,
                        double*  locs,
                        int*     dim_t,
                        int*     NNarray,
                        int*     m_t){


    int nparms = *nparms_t;
    
    int n = *n_t;
    int dim = *dim_t;
    int m  = *m_t;

    double* dySy = (double*) malloc( 2*nparms * sizeof(double) );
    memset(dySy, 0, 2*nparms*sizeof(double));

    double* dlogdet = dySy + nparms;

    double ySy = 0;
    double logdet = 0;


#pragma omp parallel
{   
    double* space     = (double*) malloc( ((m*dim) + m + m + (m*m) + (nparms*m*m) + m + m*nparms) * sizeof(double));
    double* locsub     =  space; 
    double* ysub       =  locsub  + m*dim;
    double* choli2     =  ysub    + m;
    double* covmat     =  choli2  + m;
    double* dcovmat    =  covmat  + m*m;
    double* LidSLi3    =  dcovmat + nparms*m*m; 
    double* LidSLi2    =  LidSLi3 + m;

#pragma omp for reduction(+:ySy) reduction(+:logdet) reduction(+:dySy[:nparms]) reduction(+:dlogdet[:nparms]) reduction(+:info[:nparms*nparms])
    
    for(int i = 0; i < n; ++i){
      
        int bsize = MIN(i+1, m);
        
        memset(choli2, 0.0, bsize*sizeof(double));      
        choli2[bsize-1] = 1.0;

        // Compute covariance and cholesky
        copy_in(ysub, locsub, bsize, y, locs, n, dim, NNarray, i*m);
        matern_isotropic_bivariate(covmat, dcovmat, bsize, covparms, locsub, dim);
        cholesky(covmat, bsize);          
        
	// Extract last row of inverse cholesky
        transpose_triangular_solve(covmat, choli2, bsize, 1);
     
        // ysub becomes Li*y
        triangular_solve(covmat, ysub, bsize, 1);
        
        double Liy0i2 = ysub[(bsize-1)] ;
        
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
                
		dySy[j] += 2*s1*Liy0i2 - LidSLi3i2*Liy0i2*Liy0i2;
                dlogdet[j] += LidSLi3i2;
		copy(LidSLi3, LidSLi2 + j*bsize, bsize);


	    }
		for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    info[j*nparms + i] += (dot(LidSLi2 + j*bsize, LidSLi2 + i*bsize, bsize) - 0.5*LidSLi2[j*bsize + bsize-1]*LidSLi2[i*bsize + bsize-1]);
		    info[i*nparms + j] = info[j*nparms + i];

            }}

            
        }else{
	    for(int j = 0; j < nparms; j++){
                double LidSLi = *(dcovmat + j*bsize*bsize) * (1/pow(covmat[0],2));
                dySy[j] += Liy0i2 * LidSLi * Liy0i2;
                dlogdet[j] += LidSLi; 
		LidSLi2[j] = LidSLi;
	    }


		for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                    info[j*nparms + i] += 0.5*LidSLi2[i]*LidSLi2[j];
		    info[i*nparms + j] = info[j*nparms + i];

            }}
            
            }

    }

    free(space);
}    
    
    // compute log-likelihood, gradient
    *ll = -0.5* (n*log(2.0 * M_PI) + logdet + ySy);
    for(int i=0; i<nparms; i++){
        grad[i] = -0.5*dlogdet[i] + 0.5*dySy[i];
    }

}
}

