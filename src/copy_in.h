void copy_in(double*  ysub,
             double*  locsub,
             double*  Xsub,
             int      bsize,
             double*  y,
             double*  locs,
             double*  X,
             int      n,
             int      dim,
             int      p,
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
  
  for(int k = 0 ; k < p; ++k){
    for(int j = bsize-1; j >= 0; j--){
      ii = inds[first_ind + j];
      Xsub[bsize-1-j + k*bsize] = X[(ii - 1) + k*n];
    }  
  }

}
