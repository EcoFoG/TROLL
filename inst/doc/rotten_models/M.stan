data {
  int N ;
  real dbh[N] ;
  int rotten[N] ;
}
parameters {
  real beta_3 ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] =
     beta_3 * dbh[n]
     ;
   }
   rotten ~ bernoulli(inv_logit(theta)) ;
}
