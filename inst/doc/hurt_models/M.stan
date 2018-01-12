data {
  int<lower=0> N ; // # obs
  vector<lower=0>[N] d ; // pred
  int<lower=0,upper=1> M[N] ; // obs
}
parameters {
  real theta ;
  real beta ;
  real alpha ;
}
model {
  M ~ bernoulli_logit(theta + beta*exp(-alpha*d)) ;
}
