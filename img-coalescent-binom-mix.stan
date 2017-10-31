// number of genomes == number of components == K
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
  int y[N];                             // n/o genomes gene_i is in (observation)
  // mixture components
  int K;                                // n/o genomes (== components)
  int sizes[K];                         // size of each component
  vector[K] probs;                      // probability of each component
}
parameters {
  real<lower=0> rho;                    // rate of loosing one gene
  //real<lower=0> theta;                  // number of genes gained in lineage
}
model {
  vector[K] fs;
  vector[K] lambda;
  vector[K] log_lambda;
  real theta = 1000;
  vector[K] lps;
  
  // priors
  rho ~ cauchy(0.3, 10);
  // theta ~ cauchy(1000,1000);

  // likelihood
  fs = fs_coalescent(K,rho,theta);
  lambda = fs/sum(fs);                  // normalize to one (simplex)
  log_lambda = log(lambda);

  for(i in 1:N) {                       // for each observation
    lps = log_lambda;                   // init ll with log_lambda
    for(k in 1:K) {                     // for each component
      if (y[i] > sizes[k])              // can't have more successes than trials
        lps[k] = negative_infinity();
      else                              // add ll
        lps[k] = lps[k]
          + binomial_lpmf(y[i] | sizes[k], probs[k]);
    }
    target += log_sum_exp(lps);         // update ll density
  }
}

