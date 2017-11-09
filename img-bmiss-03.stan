functions {
  // IMG coalescent frequency spectrum
  vector img_fs(int K, real rho, real theta){
    vector[K] fs = rep_vector(0, K);    // freq spectrum
    vector[K] fsp;                      // freq product

    for(k in 1:K)                       // cash prod series
      fsp[k] = (K+1 - k)/(K+rho - k);

    for(k in 1:K)                       // summarize series with length 1,2,...,k
      fs[k] = (theta/k) * prod(segment(fsp, 1, k));

    return fs;
  }

  // convert img_fs into a prob mass spectrum (0-trunc !)
  vector img_lpms(int K, real rho, real p){
    // high theta for good res for low vals at high rho
    vector[K] fs = img_fs(K, rho, 100000);
    vector[K] lw = log_softmax(log(fs));     // log weights
    vector[K] lp = rep_vector(negative_infinity(), K); // X=1..K
    real lp0 = negative_infinity(); // X=0

    if(p == 1) return lw;
    
    for(k in 1:K){
      lp0 = log_sum_exp(lp0, binomial_lpmf(0 | k, p) + lw[k]); // 0 outcome
      for(j in 1:k) // for each possible outcome (1..K)
        lp[j] = log_sum_exp(lp[j], binomial_lpmf(j | k, p) + lw[k]);
    }
    lp = lp - log1m_exp(lp0); // zero-truncate !!!
    lp = log_softmax(lp);     // not really necessary
    return lp;
  }

  // img log probability mass function
  real img_lpmf(int[] x, int K, real rho, real p){
    vector[K] lpms = img_lpms(K, rho, p); // precompute prob mass spectrum
    int N = num_elements(x);
    real lp = 0;
    for (n in 1:N)
      lp = lp + lpms[x[n]];
    return lp;
  }
}
data {
  int K;                      // n/o genomes
  int N;                      // n/o genes
  int x[N];                   // gene count (1..N)
  real<lower=0, upper=1> p;   // completeness (=detection probability)
}
parameters {
  real<lower=0> rho;                    // rate of loosing one gene
  // real<lower=0> theta;               // number of genes gained in lineage
}
transformed parameters{
  real theta = (rho * sum(x)) / (p * K);
}
model {
  // priors
  rho ~ cauchy(0.5, 10);
  // likelihood
  x ~ img(K, rho, p);
}

