data {
  int S ;
  int N ;
  int J[N] ;
  real LL[N] ;
  real LMA[N] ;
  real Nmass[N] ;
  real wsg[N] ;
}
parameters {
  real<lower=0,upper=5> sigma ;
  real<lower=0> beta_1 ;
  real<lower=0> beta_1s[S] ;
  real<lower=0,upper=50> sigma_1 ;
  real<lower=0> beta_2 ;
  real<lower=0> beta_2s[S] ;
  real<lower=0,upper=50> sigma_2 ;
  real<lower=0> beta_3 ;
}
model {
   real f[N] ;
   for(n in 1:N){
     f[n] =
     LMA[n]*beta_1s[J[n]] - 
     Nmass[n]*beta_2s[J[n]] +
     wsg[n]*beta_3
     ;
   }
   beta_1s ~ normal(beta_1, sigma_1) ;
   beta_2s ~ normal(beta_2, sigma_2) ;
   LL ~ lognormal(f, sigma) ;
}