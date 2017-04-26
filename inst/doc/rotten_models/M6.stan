data {
  int N ;
  int SP ;
  int sp[N] ;
  real dbh[N] ;
  int rotten[N] ;
}
parameters {
  real beta_2 ;
  real beta_2sp[SP] ;
  real sigma_2 ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] =
     beta_2sp[sp[n]] * dbh[n]
     ;
   }
   beta_2sp ~ normal(beta_2, sigma_2) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
