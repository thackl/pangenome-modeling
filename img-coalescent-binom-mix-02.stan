functions {
  // IMG coalescent relative frequency spectrum (sum=1)
  // no need for tau, effect is removed through normalization
  vector rfs_coalescent(int K, real rho){
    vector[K] fs = rep_vector(0, K);    // freq spectrum
    vector[K] fsp;                      // freq product
    for(k in 1:K)                       // cache prod series
      fsp[k] = (K+1 - k)/(K+rho - k);
    for(k in 1:K)                       // summarize series with length 1,2,...,k
      fs[k] = prod(segment(fsp, 1, k)) / k;
    return fs/sum(fs);
  }
}
data {
  // observations
  int N;                      // n/o genes
  int K;                      // n/o genomes (== components)
  int x[N];                   // gene counts (1..N) (genomes per gene)
  real c[K];                  // completeness
}                             
transformed data{             
  int y[K];                   // gene count frequency levels (1..K)
  real C;                     // total observed counts
  for (k in 1:K) y[k] = k;    // 'cause stan can't assign 1..K directly yet
  C = sum(x) * mean(c);
}                             
parameters {                  
  real<lower=0> rho;          // rate of loosing a single gene, time unit 2N_e
}
model {
  vector[K] lambda;
  vector[K] log_lambda;
  //real tau = 1000;
  vector[K] lps;
  
  // priors
  rho ~ cauchy(0.3, 10);

  // likelihood
  lambda = rfs_coalescent(K,rho);
  log_lambda = log(lambda);

  for(i in 1:N) {                       // for each gene count
    lps = log_lambda;                   // init ll with log_lambda
    for(k in 1:K) {                     // for each possible frequency (1..K)
      if (x[i] > y[k])                  // can't have more successes than trials
        lps[k] = negative_infinity();
      else                              // add ll
        lps[k] = lps[k]
          + binomial_lpmf(x[i] | y[k], c[k]);
    }
    target += log_sum_exp(lps);         // update ll density
  }
}
generated quantities{
  int sigma;                  // total expected counts for complete data
  real<lower=0> tau;          // average number of genes gained in 2N_e
  sigma = poisson_rng(C);
  tau = rho * sigma / K;      // reparameterized through sigma
}
