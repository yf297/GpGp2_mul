#ifndef COVFUN_H
#define COVFUN_H

#include <armadillo>
#include <vector>
#include <algorithm> // for max and min
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace arma;
using namespace boost::math;

 
void matern_isotropic(double*   covmat,
                      double*   dcovmat,
                      int       bsize,
                      double*   covparms, 
                      double*   locsub,
                      int       dim){
 
  double c0 = covparms[0];
  double c1 = covparms[1];
  double c2 = covparms[2];
  double c3 = covparms[3];
  double nugget = c0*c3;

  double normcon = c0/(pow(2.0, c2-1)*tgamma(c2));
  double eps  = 1e-8;
  double normconeps = c0/(pow(2.0, c2 + eps-1)*tgamma(c2+eps));

  double cov;

  for(int j = 0; j < bsize; ++j){
    covmat[j*bsize + j] = c0 + nugget;
    dcovmat[0*bsize*bsize + j*bsize + j] = 1 + c3;
    dcovmat[1*bsize*bsize + j*bsize + j] = 0.0;
    dcovmat[2*bsize*bsize + j*bsize + j] = 0.0;
    dcovmat[3*bsize*bsize + j*bsize + j] = 0.0 + c0;

    for(int i = (j+1); i < bsize; ++i){
      double d = 0.0;
      for(int k = 0; k < dim; ++k){
        d += (locsub[k*bsize + j] - locsub[k*bsize + i]) * (locsub[k*bsize + j] - locsub[k*bsize + i]) ;
     }
      d = sqrt(d);
      cov =  normcon * pow(d/c1, c2) * boost::math::cyl_bessel_k(c2,d/c1);
      covmat[j*bsize + i] = cov;
      dcovmat[0*bsize*bsize + j*bsize + i] = cov/c0;
      dcovmat[1*bsize*bsize + j*bsize + i] = normcon*pow(d/c1,c2)*cyl_bessel_k(c2-1,d/c1)*d/(c1*c1);
      dcovmat[2*bsize*bsize + j*bsize + i] = (normconeps*pow(d/c1,c2+eps)*cyl_bessel_k(c2+eps,d/c1) - cov)/eps;
      dcovmat[3*bsize*bsize + j*bsize + i] = 0.0;

    }
  }
}




arma::mat matern_multi(arma::vec covparms, arma::mat locs ){

    // figure out the number of components from the number of parameters
    int nparm = covparms.n_elem;
    int p;
    int maxp = 10;
    for(int i=1; i<maxp+1; i++){
	if( pow( nparm - 1.5*i*i - 2.5*i, 2.0 ) < 0.01 ){
	    p = i;
	}
    }

	
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;



    // define matrices for storing the covariance parameters
    int neach = (p*(p+1))/2;
    arma::mat cov0(p,p);
    arma::mat range(p,p);
    arma::mat smooth(p,p);
    arma::mat normcon(p,p);
    arma::vec nugget(p);
	

    int st1 = 0;
    int st2 = neach;
    int st3 = 2*neach;
    int st4 = 3*neach;
    for(int i1=0; i1<p; i1++){
        for(int i2=0; i2<=i1; i2++){
	    int i = (i1*(i1+1))/2 + i2;
	    //Rcout << i1 << "," << i2 << ":" << i << endl;
	    cov0(i1,i2) = covparms(st1 + i);
	    cov0(i2,i1) = cov0(i1,i2);
	    range(i1,i2) = covparms(st2 + i);
	    range(i2,i1) = range(i1,i2);
	    smooth(i1,i2) = covparms(st3 + i);
	    smooth(i2,i1) = smooth(i1,i2);
	    normcon(i1,i2) = cov0(i1,i2)/(pow(2.0,smooth(i1,i2)-1.0) * tgamma( smooth(i1,i2) ));
	    normcon(i2,i1) = normcon(i1,i2);
	}
	nugget(i1) = cov0(i1,i1)*covparms(st4 + i1);
    }
	
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
	    int c1 = locs(i1,dim) - 1;
	    int c2 = locs(i2,dim) - 1;
	    double rng = range(c1,c2);
	    double sm  = smooth(c1,c2);
	    
            // calculate distance
            double d = 0.0;
            for(int j=0; j<dim; j++){
                d += pow( ( locs(i1,j)/rng) - (locs(i2,j)/rng ), 2.0 );
            }
            d = pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = cov0(c1,c2);
            } else {
                // calculate covariance            
                covmat(i2,i1) = normcon(c1,c2)*pow( d, sm )*cyl_bessel_k(sm, d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget(c1); } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }
    }
    return covmat;
}



arma::cube d_matern_multi(arma::vec covparms, arma::mat locs ){

    // figure out the number of components from the number of parameters
    int nparm = covparms.n_elem;
    int p;
    int maxp = 10;
    for(int i=1; i<maxp+1; i++){
	if( pow( nparm - 1.5*i*i - 2.5*i, 2.0 ) < 0.01 ){
	    p = i;
	}
    }

    //Rcout << "p = " << p << endl;

    // fail-safe to prevent large smoothness values
    //covparms(2) = std::min( covparms(2), 8.0 );
	
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;
    double eps = 1e-8;

    //Rcout << "n = " << n << endl;

    // define matrices for storing the covariance parameters
    int neach = p*(p+1)/2;
    arma::mat cov0(p,p);
    arma::mat range(p,p);
    arma::mat smooth(p,p);
    arma::mat normcon(p,p);
    arma::mat normconeps(p,p);
    arma::vec nugget(p);

    int st1 = 0;
    int st2 = neach;
    int st3 = 2*neach;
    int st4 = 3*neach;

    for(int i1=0; i1<p; i1++){
        for(int i2=0; i2<=i1; i2++){
	    int i = (i1*(i1+1))/2 + i2;
	    //Rcout << i1 << "," << i2 << ":" << i << endl;
	    cov0(i1,i2) = covparms(st1 + i);
	    cov0(i2,i1) = cov0(i1,i2);
	    range(i1,i2) = covparms(st2 + i);
	    range(i2,i1) = range(i1,i2);
	    smooth(i1,i2) = covparms(st3 + i);
	    smooth(i2,i1) = smooth(i1,i2);
	    normcon(i1,i2) = cov0(i1,i2)/(pow(2.0,smooth(i1,i2)-1.0)*tgamma(smooth(i1,i2) ));
	    normcon(i2,i1) = normcon(i1,i2);
	    normconeps(i1,i2) = cov0(i1,i2)/(pow(2.0,smooth(i1,i2)+eps-1.0)*tgamma(smooth(i1,i2)+eps ));
	    normconeps(i2,i1) = normconeps(i1,i2);
	}
	nugget(i1) = cov0(i1,i1)*covparms(st4 + i1);
    }
	

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){

	//Rcout << i1 << "," << i2 << " : "; 
        int c1 = locs(i1,dim) - 1;
        int c2 = locs(i2,dim) - 1;
        //Rcout << c1 << "," << c2 << endl;

	// need to look up which parameters in covparms that parm(c1,c2) refer to
	int m1 = std::max(c1,c2); 
	int m2 = std::min(c1,c2); 
	int i = (m1*(m1+1))/2 + m2;

        double c0 = cov0(c1,c2);
        double rng = range(c1,c2);
        double sm  = smooth(c1,c2);
        double nc  = normcon(c1,c2);
        double nce  = normconeps(c1,c2);
        	    
        // calculate distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( ( locs(i1,j) - locs(i2,j) )/rng, 2.0 );
        }
        d = pow( d, 0.5 );
            
        double cov;        
        if( d == 0.0 ){
            cov = c0;
            dcovmat(i1,i2,st1 + i) += 1.0;
        } else {
            cov = nc*pow( d, sm )*cyl_bessel_k(sm, d);
            // variance parameter
            dcovmat(i1,i2,st1+i) += cov/covparms(st1+i);
            // range parameter
            dcovmat(i1,i2,st2+i) += nc*pow(d,sm)*cyl_bessel_k(sm-1.0, d)*d/rng;
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,st3+i) += (nce*pow(d,sm+eps)*cyl_bessel_k(sm+eps,d) - cov)/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,st1+i) += covparms(st4+c1);
            dcovmat(i1,i2,st4+c1) += covparms(st1+i); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}


double B(double a, double b){
	double c =  tgamma(a)*tgamma(b)/tgamma(a + b);     
	return(c);
}

void matern_isotropic_bivariate(double*   covmat,
                                double*   dcovmat,
                                int       bsize,
                                double*   covparms, 
                                double*   locsub,
                                int       dim){

	
	int nparms = 11;
	arma::vec a_covparms(covparms, nparms);
	arma::mat a_locs(locsub, bsize, dim);

	arma::mat a_covmat;
	a_covmat = matern_multi(a_covparms, a_locs);
	
	arma::cube a_dcovmat(bsize, bsize, nparms);
	a_dcovmat = d_matern_multi(a_covparms, a_locs);

	for(int i = 0; i < bsize; ++i){
            for(int j = 0; j < bsize; ++j){ 
    	        covmat[j*bsize + i] =  a_covmat(i,j);
		for(int p = 0; p < nparms; ++p){
		    dcovmat[p*bsize*bsize + j*bsize + i] = a_dcovmat(i,j,p);
	    	}
  	    }
	}

}




extern "C" {
void con(double* cons, int* d_t, double* covparms, int* nparms_t){
	
	int nparms = *nparms_t;
	int d = *d_t;

	arma::vec a_covparms(covparms, nparms);
				
	double delta_A = a_covparms(4);
	double delta_B = a_covparms(1);
	
	double nu11 = a_covparms(3);
	double nu22 = a_covparms(5);
	double nu12 = delta_A + (nu11 + nu22)/2;
	
	double a11 = 1/a_covparms(0);
	double a22 = 1/a_covparms(2);
	double a12 = pow(delta_B + (pow(a11,2) + pow(a22,2))/2, 0.5) ;	
	
	double tau1 = pow(B(nu12, d/2), 2) / pow(B( (nu11 + nu22)/2 , d/2), 2);
	double tau2 = pow(a11*a22/pow(a12,2), 2*delta_A);
	double tau3 = (pow(tgamma( (nu11 + nu22)/2), 2)/ pow(a12, 2*(nu11 + nu22) ) ) / ( (tgamma(nu11)/pow(a11,2*nu11)) *  (tgamma(nu22)/pow(a22,2*nu22)) ); 

	*cons = tau1*tau2*tau3;	


	}
}




#endif
