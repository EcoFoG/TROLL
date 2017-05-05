data {
  int<lower=0> N ; // # obs
  vector<lower=0>[N] d ; // pred
  vector<lower=0,upper=1>[N] S ; // pred sylv or not
  int<lower=0,upper=1> M[N] ; // obs
}
parameters {
  real theta ;
  real beta ;
}
model {
  for(n in 1:N)
    M[n] ~ bernoulli_logit(theta + S[n]*beta*exp(-(0.001 + d[n]))) ;
}
