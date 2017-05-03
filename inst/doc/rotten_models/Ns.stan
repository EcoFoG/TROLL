data {
  int<lower=0> N ;
  real<lower=0> dbh[N] ;
  int<lower=0,upper=1> Pr[N] ;
  real<lower=0> Vf[N] ;
  int<lower=0> S ;
  int<lower=0> species[N] ;
}
parameters {
  real<lower=0, upper=1> rho ;
  real<lower=0> beta_s[S] ;
  real<lower=0> sigma_s ;
  real<lower=-min(beta_s)> beta ;
  real<lower=0> sigma ;
}
model {
  beta_s ~ normal(0, sigma_s) ;
  for(n in 1:N)
    Vf[n] ~ lognormal(log(((beta + beta_s[species[n]])*dbh[n]*dbh[n])*(1-Pr[n]*rho)), sigma) ;
}
