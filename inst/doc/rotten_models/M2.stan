data {
  int N ;
  int P ;
  int plot[N] ;
  real wsg[N] ; 
  real dbh[N] ;
  int rotten[N] ;
}
parameters {
  real beta_0p[P] ;
  real beta_1p[P] ;
  real beta_2p[P] ;
  real beta_0 ;
  real beta_1 ;
  real beta_2 ;
  real sigma_0 ;
  real sigma_1 ;
  real sigma_2 ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] = beta_0p[plot[n]] +
     beta_1p[plot[n]] * wsg[n] +
     beta_2p[plot[n]] * dbh[n]
     ;
   }
   beta_0p ~ normal(beta_0, sigma_0) ;
   beta_1p ~ normal(beta_1, sigma_1) ;
   beta_2p ~ normal(beta_2, sigma_2) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
