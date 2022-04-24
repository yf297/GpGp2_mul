#ifndef COPY_IN_H
#define COPY_IN_H
void copy_in(double*  ysub,
             double*  locsub,
             int      bsize,
             double*  y,
             double*  locs,
             int      n,
             int      dim,
             int*     inds,
             int      first_ind){
  
  

    int ii;
  
    for(int j = bsize-1; j >= 0; --j){
      ii = inds[first_ind + j];
      ysub[bsize-1-j] = y[ii - 1];
  }
  
  for(int k = 0 ; k < dim; ++k){
    for(int j = bsize-1; j >= 0; j--){
      ii = inds[first_ind + j];
      locsub[bsize-1-j + k*bsize] = locs[(ii - 1) + k*n];
    }  
  }
  

}



#endif
