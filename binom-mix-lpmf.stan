functions {
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
model {
  int x[3] = {1,2,3};
  vector[3] weights = [0,0,1]';
  real c = 0.8;
  real ll;
  real lls[3];
  real lx;
  real lxs[3];
  //
  ll = binomial_lpmf(x | 3, c);
  print("");
  print("-- Binomial LPMF --");
  print("ll: ", ll);

  print("");
  print("-- Binom Mix LPMF --");
  lx = binom_mix_lpmf(x | weights, c);
  print("lx: ", lx);
}
