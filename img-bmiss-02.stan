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
  // "miss" data through bernoulli/binomial
  vector bmiss(vector x, real p){
    int K = num_elements(x);
    vector[K] log_r = log(x/sum(x));
    vector[K] lp = rep_vector(negative_infinity(), K);
    for(k in 1:K){
      for(j in 1:k)
        lp[j] = log_sum_exp(lp[j], binomial_lpmf(j | k, p) + log_r[k]);
    }
    return exp(lp) * sum(x);
  }
}
data {
  int K;                      // n/o genomes
  int M;                      // total number of gene family members (M/p)
  real y[K];                  // gene count frequency levels (1..K)
  real<lower=0, upper=1> p;   // completeness (=~ detection probability)
}
parameters {
  real<lower=0> rho;                    // rate of loosing one gene
  real<lower=0> theta;                  // number of genes gained in lineage
}
model {
  vector[K] fs;
  vector[K] fs_miss;
  real theta_mu;
  
  // priors
  rho ~ cauchy(0.3, 10);
  // theta_mu = rho * (M/p) / K;
  // theta ~ normal(theta_mu, sqrt(theta_mu));
  // print("r: ", rho, ", t: ", theta);

  // generator
  fs = img_fs(K,rho,theta);
  // print("Full FS: ", round(fs));
  if(p < 1){
    fs_miss = bmiss(fs, p);
  }else{
    fs_miss = fs;
  }
  // print("Miss FS: ", round(fs_miss));
  y ~ normal(fs_miss, 2);
}
