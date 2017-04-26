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
  real beta_0 ;
  real beta_1 ;
  real beta_2 ;
  real sigma ;
}
model {
   real theta[N] ;
   for(n in 1:N){
     theta[n] = 
     beta_0p[plot[n]] +
     beta_1 * wsg[n] +
     beta_2 * dbh[n]
     ;
   }
   beta_0p ~ normal(beta_0, sigma) ;
   rotten ~ bernoulli(inv_logit(theta)) ;
}
