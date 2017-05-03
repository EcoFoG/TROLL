data {
  int<lower=0> N ;
  real<lower=0> dbh[N] ;
  int<lower=0,upper=1> Pr[N] ;
  real<lower=0> Vf[N] ;
}
parameters {
  real<lower=0> beta ;
  real<lower=0, upper=1> rho ;
  real<lower=0> sigma ;
}
model {
  for(n in 1:N)
    Vf[n] ~ lognormal(log((beta*dbh[n]*dbh[n])*(1-Pr[n]*rho)), sigma) ;
}
