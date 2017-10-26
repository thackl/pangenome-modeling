// number of genomes == number if components == K
// number of genes (COGs) == N
// number genomes gene was seen in: observation
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
  // observations
  int N;                                // n/o genes
  real y[N];                             // n/o genomes gene_i is in (observation)
  // mixture components
  int K;                                // n/o genomes (== components)
  // int sizes[K];                         // size of each component
  // vector[K] probs;                      // probability of each component
}
parameters {
  real<lower=0> rho;                    // rate of loosing one gene
  real<lower=0> theta;                  // number of genes gained in lineage
}
model {
  vector[K] fs;

  // priors
  rho ~ cauchy(0.3, 10);
  theta ~ cauchy(1000,1000);

  // generator
  fs = fs_coalescent(K,rho,theta);
  y ~ cauchy(fs, 2);
}
