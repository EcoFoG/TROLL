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
  real<lower=0> beta ;
  real<lower=0> sigma ;
  real<lower=-beta> beta_p[P] ;
  real<lower=0> sigma_p ;
  real<lower=-beta-min(beta_p)> beta_s[S] ;
  real<lower=0> sigma_s ;
  real<lower=0,upper=1> theta ;
  real<lower=-theta, upper=1-theta> theta_p[P] ;
  real<lower=0> sigma_rp ;
}
transformed parameters{
  real<lower=0, upper=1> rho[N] ;
  for(n in 1:N)
    rho[n] = theta*dbh[n]*dbh[n] + theta_p[plot[n]] ;
}
model {
  theta_p ~ normal(0, sigma_rp) ;
  beta_p ~ normal(0, sigma_p) ;
  beta_s ~ normal(0, sigma_s) ;
  for(n in 1:N)
    Vf[n] ~ lognormal(log(((beta + beta_p[plot[n]] + beta_s[species[n]])*dbh[n]*dbh[n])*(1-Pr[n]*rho[n])), sigma) ; 
}
