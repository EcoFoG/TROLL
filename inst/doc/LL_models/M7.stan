data {
  int S ;
  int N ;
  int J[N] ;
  real LL[N] ;
  real LMA[N] ;
}
parameters {
  real<lower=0,upper=5> sigma ;
  real<lower=0,upper=50> beta_0 ;
  real<lower=0> beta_1 ;
  real<lower=0> beta_1s[S] ;
  real<lower=0,upper=50> sigma_1 ;
}
model {
   real f[N] ;
   for(n in 1:N){
     f[n] = beta_0 + LMA[n]*beta_1s[J[n]] ;
   }
   beta_1s ~ normal(beta_1, sigma_1) ;
   LL ~ lognormal(f, sigma) ;
}
