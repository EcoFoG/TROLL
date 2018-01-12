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
  real<lower=0,upper=3> beta_3 ;
  real<lower=0,upper=3> beta_3s[S] ;
  real<lower=0,upper=50> sigma_3 ;
}
model {
   real f[N] ;
   for(n in 1:N){
     f[n] = beta_0 +
     pow(LMA[n], beta_3s[J[n]])
     ;
   }
   beta_3s ~ normal(beta_3, sigma_3) ;
   LL ~ lognormal(f, sigma) ;
}
