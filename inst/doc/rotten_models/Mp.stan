data {
  int N ;
  int P ;
  int plot[N] ;
  real dbh[N] ;
  int rotten[N] ;
}
parameters {
  real beta_0 ;
  real beta_0p[P] ;
  real beta_3 ;
  real sigma_p ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] =
     beta_0p[plot[n]] +
     beta_3 * dbh[n]
     ;
   }
   beta_0p ~ normal(beta_0, sigma_p) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
