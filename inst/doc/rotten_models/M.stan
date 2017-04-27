data {
  int<lower=0> N ; // number of obs (trees)
  vector[N] dbh ; // trees dbh
  int<lower=0, upper=1> rotten[N] ; // outcome (rotten)
}
parameters {
  real beta_0 ;
  real beta_1 ;
}
model {
  rotten ~ bernoulli_logit(beta_0 + beta_1 * dbh) ;
}
