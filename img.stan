functions {
  // IMG coalescent frequency spectrum
  vector fs_coalescent(int K, real rho, real theta){
    // declare
    vector[K] fs;                       // freq spectrum
    vector[K] fsp;                      // freq product

    // init
    fs = rep_vector(0, K);              // init all 0
    for(k in 1:K)                       // cash prod series
      fsp[k] = (K+1 - k)/(K+rho - k);

    for(k in 1:K)                       // summarize series with length 1,2,...,k
      fs[k] = (theta/k) * prod(segment(fsp, 1, k));

    return fs;
  }
}
data {
  int K;
  real rho;
  real theta;
}
model {
  vector[K] fs;
  fs = fs_coalescent(K,rho,theta);
  print("Frequeny Spectrum: ", fs);
}
