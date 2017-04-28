data {
  int<lower=0> N ;
  real dbh[N] ;
  int<lower=0, upper=1> Pr[N] ;
  real Vf[N] ;
}
parameters {
  vector[2] beta ;
  real sigma ;
  vector[1] theta ;
}
model {
  for(n in 1:N)
    Vf[n] ~ normal((beta[1] + beta[2]*dbh[n]*dbh[n])*(1-Pr[n]*theta[1]), sigma*dbh[n]) ;
}
