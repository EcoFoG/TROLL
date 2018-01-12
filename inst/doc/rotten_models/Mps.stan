data {
  int<lower=0> N ; // number of obs (trees)
  int<lower=0> P ; // number of plots
  int plot[N] ; // plots id
  int<lower=0> S ; // number of species
  int species[N] ; // species id
  real dbh[N] ; // trees dbh
  int<lower=0, upper=1> rotten[N] ; // outcome (rotten)
}
parameters {
  real beta_0 ;
  real beta_1 ;
  real beta_2p[P] ;
  real sigma_p ;
  real beta_3s[S] ;
  real sigma_s ;
}
model {
  beta_2p ~ normal(0, sigma_p) ;
  beta_3s ~ normal(0, sigma_s) ;
  for(n in 1:N)
    rotten[n] ~ bernoulli_logit(beta_0 + beta_1 * dbh[n] + beta_2p[plot[n]] + beta_3s[species[n]]) ;
}

