data {
  // observations
  int N;                                // n/o observations
  int y[N];                             // observations
  // mixture components
  int K;                                // n/o components
  int sizes[K];                         // size of each component
  vector[K] probs;                      // probability of each component
  // dirichlet prior params
  simplex[K] phi;                       // expected mean
  real<lower=0> kappa;                  // strength as count
}
parameters {
  simplex[K] theta;                     // mixture weights to estimate
}
transformed parameters {
  vector[K] alpha = kappa * phi;        // reparameterized dirichlet alpha
}
model {
  // declare
  real log_theta[K];
  real lps[K];

  // priors
  // theta ~ dirichlet(rep_vector(2.0, K));
  // Reparameterization of Dirichlet Priors (Stan Manual p282)
  theta ~ dirichlet(alpha);

  // cache log calculation, slight speed-up from
  // (Estimating Parameters of a Mixture, Stan Manual p192)
  for (k in 1:K)
    log_theta[k] = log(theta[k]);

  // likelihood
  for(i in 1:N) {                       // for each observation
    lps = log_theta;                    // init ll with log_theta
    for(k in 1:K) {                     // for each component
      if (y[i] > sizes[k])              // can't have more successes than trials
        lps[k] = negative_infinity();
      else                              // add ll
        lps[k] = lps[k] + binomial_lpmf(y[i] | sizes[k], probs[k]);
    }
    target += log_sum_exp(lps);         // update ll density
  }
}
