data {
  int<lower=0> N ;
  real<lower=0> dbh[N] ;
  int<lower=0,upper=1> Pr[N] ;
  real<lower=0> Vf[N] ;
  int<lower=0> P ;
  int<lower=0, upper=P> plot[N] ;
}
parameters {
  real<lower=0> beta ;
  real<lower=0> sigma ;
  real<lower=0, upper=1> rho ;
  real<lower=-beta> beta_p[P] ;
  real<lower=0> sigma_p ;
}
model {
  beta_p ~ normal(0, sigma_p) ;
  for(n in 1:N)
    Vf[n] ~ lognormal(log(((beta + beta_p[plot[n]])*dbh[n]*dbh[n])*(1-Pr[n]*rho)), sigma) ;
}