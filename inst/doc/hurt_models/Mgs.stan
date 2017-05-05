data {
  int<lower=0> N ; // # obs
  vector<lower=0>[N] d_gaps ; // pred
  int<lower=0,upper=1> M[N] ; // obs
}
parameters {
  real alpha ;
}
transformed parameters {
  vector<lower=0,upper=1>[N] p ;
  p = exp(-alpha*d_gaps) ;
}
model {
  M ~ bernoulli(p) ;
}
