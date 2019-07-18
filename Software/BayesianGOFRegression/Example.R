#####################
### Required libraries
#####################
require(invgamma)
require(MASS)
require(Matrix)
require(mvtnorm)
require(splines)

#####################
### Read source
#####################
source(url("https://raw.githubusercontent.com/anfebar/anfebar.github.io/master/Software/BayesianGOFRegression/MCMCAlgorithm.R"))

##########################################
##########################################
# EXAMPLE WITH A MODEL CORRECTLY SPECIFIED
##########################################
##########################################


#####################
### Simulate Data
#####################
## Predictors
n = 100 # sample size
p1 = 5  # number of continuous predictors
p2 = 5  # number of binary predictors
S.X = matrix(1,p1,p1)
for(j1 in 1:p1)
{
  for(j2 in 1:p1) S.X[j1,j2] = 0.7^abs(j1-j2)
}
S.X[1,-1] = 0
S.X[-1,1] = 0
X = cbind(1,rmvnorm(n,rep(0,p1),S.X))
X = cbind(X,sapply(rep(0.5,len=p2),function(p)rbinom(n,1,p)))
## Response variable
p = ncol(X)
m = -2+X%*%rep(0.15,p)
mu = 1/(1+exp(m))
phi= 1.5
y=rbeta(n,mu*phi,(1-mu)*phi)
for(i in 1:n)
{
  while(y[i]==1 | y[i]==0 | is.na(y[i]))
  y[i]=rbeta(1,mu[i]*phi,(1-mu[i])*phi)
}

#####################
### Fit right model
#####################
theta=c(-1.75,rep(0.15,p-1),phi)
log_like=function(theta,y,x)
{
   m = X%*%theta[1:p]
   mu = 1/(1+exp(m))
   phi=theta[p+1]
    -sum(dbeta(y,mu*phi,(1-mu)*phi,log=T))
}
theta=optim(theta,log_like,y=y,x=x)$par;theta


#####################
### Compute universal residuals
#####################
m = X%*%theta[1:p]
mu = 1/(1+exp(m))
phi=theta[p+1]
  
u1 = pbeta(y,mu*phi,(1-mu)*phi)
 
u1[u1==1] = 0.999
u1[u1==0] = 0.001
u = qnorm(u1)
  

#####################
### Transform continous predictors
### to increase atoms' flexilibity
#####################
## Transform continous predictors
aux = cbind(rep(1,n))
for(j in 2:(p1+1)) aux = cbind(aux,bs(X[,j]))
Xbs = cbind(aux,X[,(p1+2):p])


#####################
## Create dataframe -- named column containing the universal residuals as "u"
#####################
Data = data.frame(u=u,X=Xbs)


#####################
### Set mcmc parameters --  H: DPM truncation level
#####################
MCMC = list(nburn = 5000, nskip = 10, nsave = 1000, H = 15)


#####################
### Compute Bayes factor
#####################
BF = f.MCMC(Data,MCMC)
BF # BF should be larger as n increases



##########################################
##########################################
# EXAMPLE WITH A MODEL CORRECTLY SPECIFIED
##########################################
##########################################

#####################
### Simulate Data
#####################
## Predictors
n = 100 # sample size
p1 = 5  # number of continuous predictors
p2 = 5  # number of binary predictors
S.X = matrix(1,p1,p1)
for(j1 in 1:p1)
{
  for(j2 in 1:p1) S.X[j1,j2] = 0.7^abs(j1-j2)
}
S.X[1,-1] = 0
S.X[-1,1] = 0
X = cbind(1,rmvnorm(n,rep(0,p1),S.X))
X = cbind(X,sapply(rep(0.5,len=p2),function(p)rbinom(n,1,p)))
## Response variable
p = ncol(X)
m = -2+X%*%rep(0.15,p)
mu = 1/(1+exp(m))
m = X[,c(1,2,ncol(X))]%*%rep(1,3)
phi= 1.5+(exp(m))
y=rbeta(n,mu*phi,(1-mu)*phi)
for(i in 1:n)
{
  while(y[i]==1 | y[i]==0 | is.na(y[i]))
  y[i]=rbeta(1,mu[i]*phi,(1-mu[i])*phi)
}


#####################
### Fit wrong model
#####################
theta=c(-1.75,rep(0.15,p-1),mean(phi))
log_like=function(theta,y,x)
{
   m = X%*%theta[1:p]
   mu = 1/(1+exp(m))
   phi=theta[p+1]
    -sum(dbeta(y,mu*phi,(1-mu)*phi,log=T))
}
theta=optim(theta,log_like,y=y,x=x)$par;theta


#####################
### Compute universal residuals
#####################
m = X%*%theta[1:p]
mu = 1/(1+exp(m))
phi=theta[p+1]
  
u1 = pbeta(y,mu*phi,(1-mu)*phi)
 
u1[u1==1] = 0.999
u1[u1==0] = 0.001
u = qnorm(u1)
  

#####################
### Transform continous predictors
### to increase atoms' flexilibity
#####################
## Transform continous predictors
aux = cbind(rep(1,n))
for(j in 2:(p1+1)) aux = cbind(aux,bs(X[,j]))
Xbs = cbind(aux,X[,(p1+2):p])


#####################
## Create dataframe -- named column containing the universal residuals as "u"
#####################
Data = data.frame(u=u,X=Xbs)


#####################
### Set mcmc parameters --  H: DPM truncation level
#####################
MCMC = list(nburn = 5000, nskip = 10, nsave = 1000, H = 15)


#####################
### Compute Bayes factor
#####################
BF = f.MCMC(Data,MCMC)
BF # BF should be smaller as n increases
