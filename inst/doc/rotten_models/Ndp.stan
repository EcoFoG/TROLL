data {
  int<lower=0> N ;
  int<lower=0> P ;
  int<lower=0> plot[N] ;
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
}
model {
  theta_p ~ normal(0, sigma_p) ;
  for(n in 1:N)
    Vf[n] ~ normal((beta[1] + beta[2]*dbh[n]*dbh[n])*(1-Pr[n]*(theta[1] + theta[2]*dbh[n]*dbh[n] + theta_p[plot[n]])), sigma*dbh[n]) ;
}
