data {
  int N ;
  int SP ;
  int sp[N] ;
  real dbh[N] ;
  int rotten[N] ;
}
parameters {
  real beta_0 ;
  real beta_0sp[SP] ;
  real beta_3 ;
  real sigma_sp ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] =
     beta_0sp[sp[n]] +
     beta_3 * dbh[n]
     ;
   }
   beta_0sp ~ normal(beta_0, sigma_sp) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
