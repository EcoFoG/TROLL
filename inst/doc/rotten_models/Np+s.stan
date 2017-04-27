data {
  int N ;
  int P ;
  int S ;
  int plot[N] ;
  int species[N] ; 
  int P_r[N] ;
  real dbh[N] ;
  real V_f[N] ;
}
parameters {
  real beta_0 ;
  real beta_0p[P] ;
  real beta_1s[S] ;
  real beta_2 ;
  real sigma ;
  real sigma_p ;
  real sigma_s ;
  real V_r[N] ;
}
model {
  real mu[N] ;
  for(n in 1:N){
     mu[n] =
     beta_0p[plot[n]] +
     beta_1s[species[n]] +
     beta_2 * dbh[n] -
     V_r[n]*P_r[n]
     ;
   }
   beta_0p ~ normal(beta_0, sigma_p) ;
   beta_1s ~ normal(0, sigma_s) ;
   V_r ~ normal(mu, sigma) ;
}
