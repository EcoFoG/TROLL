data {
  int<lower=0> N ; // # obs
  int<lower=0,upper=1> M[N] ; // obs
}
parameters {
  real theta ;
}
model {
  M ~ bernoulli_logit(theta) ;
}
