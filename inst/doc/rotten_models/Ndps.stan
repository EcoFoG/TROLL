data {
  int<lower=0> N ;
  int<lower=0> P ;
  int<lower=0> plot[N] ;
  int<lower=0> S ;
  int<lower=0> species[N] ;
  real dbh[N] ;
  int<lower=0, upper=1> Pr[N] ;
  real Vf[N] ;
}
parameters {
  vector[2] beta ;
  vector[2] theta ;
  real sigma ;
  vector[P] theta_p ;
  real sigma_p ;
  vector[S] theta_s ;
  real sigma_s ;
}
model {
  theta_p ~ normal(0, sigma_p) ;
  theta_s ~ normal(0, sigma_s) ;
  for(n in 1:N)
    Vf[n] ~ normal((beta[1] + beta[2]*dbh[n]*dbh[n])*(1-Pr[n]*(theta[1] + theta[2]*dbh[n]*dbh[n] + theta_p[plot[n]] + theta_s[species[n]])), sigma*dbh[n]) ;
}
