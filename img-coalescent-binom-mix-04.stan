functions {
  // IMG coalescent frequency spectrum
  vector fs_coalescent(int K, real rho, real tau){
    vector[K] fs = rep_vector(0, K);    // freq spectrum
    vector[K] fsp;                      // freq product
    for(k in 1:K)                       // cache prod series
      fsp[k] = (K+1 - k)/(K+rho - k);
    for(k in 1:K)                       // summarize series with length 1,2,...,k
      fs[k] = tau/k * prod(segment(fsp, 1, k));
    return fs;
  }
  // lpmf for a binomial mixture spectrum (1..N trial classes)
  real binom_mix_lpmf(int[] n, vector lambda, real theta){
    int nmax = size(n);
    int Nmax = num_elements(lambda);
    int Omax = Nmax + 1;                // include 0 in outcomes
    int o;
    real LP[Omax, Nmax];
    vector[Omax] lps;
    real lp[nmax];
    vector[Nmax] std_lambda;
    vector[Nmax] log_lambda;
    std_lambda = lambda / sum(lambda);  // normalize to simplex
    log_lambda = log(std_lambda);       // cache log
    // get a weighted log pmf matrix for all n/o trials vs outcomes
    for(oi in 1:Omax){                  // for all possible outcomes
      o = oi - 1;                       // actual outcome (incl. 0)
      for(N in 1:Nmax){                 // for all possible number of trials
        if(o > N){ // cannot have more successes than trials
          LP[oi,N] = negative_infinity();
        }else{
          LP[oi,N] = binomial_lpmf(o | N, theta) + log_lambda[N];
        }
      }
    }
    // summarize over outcomes (log space)
    for(oi in 1:Omax)                   // o is actually oi-1
      lps[oi] = log_sum_exp(LP[oi]);
    // summarize lps of observed data
    for (i in 1:nmax)
      lp[i] = lps[n[i]+1];
    return sum(lp);
  }
}
data {
  int K;                      // n/o genomes
  int N;                      // number of genes (families)
  int x[N];                   // gene count observations (1..N)
  // int y[K];                   // gene count frequency levels (1..K)
  real<lower=0, upper=1> c;   // completeness
}
parameters {                  
  real<lower=0.1> rho;        // rate of loosing a single gene, time unit 2N_e
  vector<lower=0>[K] y_exp;
  // real<lower=0> tau;          // average number of genes gained in 2N_e
}
transformed parameters {
  vector[K] y_img;
  real tau = 200;
  y_img = fs_coalescent(K, rho, tau);
  print("y_img: ", y_img);
}
model {
  // priors
  rho ~ cauchy(0.5, 10);
  y_exp ~ normal(y_img, 1);   // expected complete freqs normal around img
  print("y_exp: ", y_exp);
  // tau ~ cauchy(1000, 500);
  print("rho: ", rho);
  x ~ binom_mix_lpmf(y_exp, c);
}
