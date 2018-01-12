data {
  int<lower=0> N;
  vector[N] dbh;
  int<lower=0,upper=1> rotten[N];
}
parameters {
  vector[2] beta;
}
model {
  rotten ~ bernoulli_logit(beta[1] + beta[2] * dbh);
}
