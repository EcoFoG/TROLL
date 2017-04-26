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
  real beta_0p[P] ;
  real beta_1s[S] ;
  real beta_2 ;
  real sigma_p ;
  real sigma_s ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] =
     beta_0p[plot[n]] +
     beta_1s[species[n]] +
     beta_2 * dbh[n]
     ;
   }
   beta_0p ~ normal(beta_0, sigma_p) ;
   beta_1s ~ normal(0, sigma_s) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
