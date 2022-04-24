#ifndef LINALG_H
#define LINALG_H

#include <cblas.h>
#include <lapacke.h>
  
void cholesky(double*  a, 
              int      n){
  
    
   LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L', n, a, n);
    
}


void triangular_solve(double* a, 
                      double*  b, 
           	      int      n,
                      int      p){
  
        
   LAPACKE_dtrtrs(LAPACK_COL_MAJOR,'L', 'N', 'N', n, p, a,n, b, n);

}

void symmetric_solve(double* a, 
                     double*  b, 
                     int      n,
                     int      p){
  
  int ipiv[n];
        
   LAPACKE_dsysv(LAPACK_COL_MAJOR,'L', n, p, a,n,ipiv, b, n);

}


void transpose_triangular_solve(double* a, 
             			double*  b, 
             			int      n,
	     			int      p){
  
        
   LAPACKE_dtrtrs(LAPACK_COL_MAJOR,'L', 'T', 'N', n, p, a,n, b, n);

}


double dot(double* a,
         double* b,
         int     n){
    
   return cblas_ddot(n, a, 1, b, 1);
    
}


void symmetric_matrix_vector(double* A,
           double* x,
           int M,
           double* y){
    
   cblas_dsymv(CblasColMajor, CblasLower, M, 1.0, A, M, x,1,0.0,y,1); 
    
}


void transpose_matrix_vector(double* A,
   		       	     double* x,
           		     int M,
                             int N,
                             double* y){
    
   cblas_dgemv(CblasColMajor,CblasTrans, M, N, 1.0, A, M, x,1,0.0,y,1); 
    
}


void outer_product(double s,
           double* x,
           double* y,
           double* A,
           int p1,
           int p2){
    
  cblas_dger(CblasColMajor, p1, p2 , s, x, 1, y, 1, A, p1); 
    
}


void scale(double* x,
           double y,
           int p){
    cblas_dscal(p,y,x,1);
    
}


void copy(double* x,
          double* y,
          int p){
    
    cblas_dcopy(p,x,1,y,1);
    
}


void vector_add(double* x,
         double a,
         double* y,
         int p){
    
    cblas_daxpy(p,a,x,1,y,1);
    
    }

#endif
