data {
  int N ;
  int SP ;
  int sp[N] ;
  real dbh[N] ;
  int rotten[N] ;
}
parameters {
  real beta_0 ;
  real beta_1 ;
  real beta_2 ;
  real beta_0sp[SP] ;
  real beta_1sp[SP] ;
  real beta_2sp[SP] ;
  real sigma_0 ;
  real sigma_1 ;
  real sigma_2 ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] =
     beta_0sp[sp[n]] +
     beta_1sp[sp[n]] +
     beta_2sp[sp[n]] * dbh[n]
     ;
   }
   beta_0sp ~ normal(beta_0, sigma_0) ;
   beta_1sp ~ normal(beta_1, sigma_1) ;
   beta_2sp ~ normal(beta_2, sigma_2) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
