data {
  int<lower=0> N ;
  real dbh[N] ;
  real Vf[N] ;
}
parameters {
  vector[2] beta;
  real sigma;
}
model {
  for(n in 1:N)
    Vf[n] ~ normal(beta[1] + beta[2]*dbh[n]*dbh[n], sigma*dbh[n]) ;
}
