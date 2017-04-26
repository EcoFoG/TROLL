data {
  int N ;
  int P ;
  int S ;
  int plot[N] ;
  int species[N] ; 
  real dbh[N] ;
  int rotten[N] ;
}
parameters {
  real beta_0 ;
  real beta_1 ;
  real beta_2 ;
  real beta_3 ;
  real beta_1p[P] ;
  real beta_2s[S] ;
  real sigma_p ;
  real sigma_s ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] = beta_0 +
     beta_1p[plot[n]] +
     beta_2s[species[n]] +
     beta_3 * dbh[n]
     ;
   }
   beta_1p ~ normal(beta_1, sigma_p) ;
   beta_2s ~ normal(beta_2, sigma_s) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
