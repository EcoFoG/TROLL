data {
  int S ;
  int N ;
  int J[N] ;
  real LL[N] ;
  real LMA[N] ;
  real Nmass[N] ;
}
parameters {
  real<lower=0,upper=5> sigma ;
  real<lower=0,upper=50> beta_0 ;
  real<lower=0> beta_1 ;
  real<lower=0> beta_1s[S] ;
  real<lower=0,upper=50> sigma_1 ;
  real<lower=0> beta_2 ;
  real<lower=0> beta_2s[S] ;
  real<lower=0,upper=50> sigma_2 ;
  real<lower=0,upper=3> beta_3 ;
  real<lower=0,upper=3> beta_3s[S] ;
  real<lower=0,upper=50> sigma_3 ;
  real<lower=0,upper=3> beta_4 ;
  real<lower=0,upper=3> beta_4s[S] ;
  real<lower=0,upper=50> sigma_4 ;
}
model {
   real f[N] ;
   for(n in 1:N){
     f[n] = beta_0 +
     beta_1s[J[n]]*pow(LMA[n], beta_3s[J[n]]) - 
     beta_2s[J[n]]*pow(Nmass[n], beta_4s[J[n]])
     ;
   }
   beta_1s ~ normal(beta_1, sigma_1) ;
   beta_2s ~ normal(beta_2, sigma_2) ;
   beta_3s ~ normal(beta_3, sigma_3) ;
   beta_4s ~ normal(beta_4, sigma_4) ;
   LL ~ lognormal(f, sigma) ;
}
