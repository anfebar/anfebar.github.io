/*
*Simple normal regression example
*/

data {
  int N; //the number of observations
  int K; //the number of columns in the model matrix
  real y[N]; //the response
  matrix[N,K] X; //the model matrix
}
parameters {
  vector[K-1] xi; //the regression parameters
  //vector[K] xi; //the regression parameters
  real<lower=0> sigma; //the standard deviation
  real mu; //the standard deviation
  real<lower=0>  lambda; //the standard deviation
}
transformed parameters {
  vector[K] xi_centered;
  //xi_centered = append_row(xi,-sum(xi));
  xi_centered = append_row(0,xi);
}
model {  

  vector[N] linpred;
  linpred = X*xi_centered;

  xi ~ normal(mu,sigma*sqrt(lambda));//prior for the slopes following Gelman 2008
  //xi ~ double_exponential(mu,sigma*sqrt(lambda));//prior for the slopes following Gelman 2008
  lambda ~ cauchy(0,1);//prior for the slopes following Gelman 2008
  mu ~ normal(0,3);
  
  y ~ normal(linpred,sigma);
}

