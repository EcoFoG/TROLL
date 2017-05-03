data {
  int<lower=0> N ;
  real<lower=0> dbh[N] ;
  int<lower=0,upper=1> Pr[N] ;
  real<lower=0> Vf[N] ;
  int<lower=0> P ;
  int<lower=0> plot[N] ;
  int<lower=0> S ;
  int<lower=0> species[N] ;
}
parameters {
  real<lower=0, upper=1> theta_s[S] ;
  real<lower=0> sigma_rs ;
  real<lower=-min(theta_s), upper=1-max(theta_s)> theta ;
  real<lower=0> beta_p[P] ;
  real<lower=0> sigma_p ;
  real<lower=0> beta_s[S] ;
  real<lower=0> sigma_s ;
  real<lower=-min(beta_p)-min(beta_s)> beta ;
  real<lower=0> sigma ;
}
transformed parameters{
  real<lower=0, upper=1> rho[N] ;
  for(n in 1:N)
    rho[n] = (theta + theta_s[species[n]])*dbh[n]*dbh[n] ;
}
model {
  theta_s ~ normal(0, sigma_rs) ;
  beta_p ~ normal(0, sigma_p) ;
  beta_s ~ normal(0, sigma_s) ;
  for(n in 1:N)
    Vf[n] ~ lognormal(log(((beta + beta_p[plot[n]] + beta_s[species[n]])*dbh[n]*dbh[n])*(1-Pr[n]*rho[n])), sigma) ; 
}
