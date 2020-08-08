f.MCMC = function(Data,MCMC)
{
  # Load data
  X = as.matrix(Data[,-which(colnames(Data)=="u")])
  p = ncol(X)
  n = nrow(X)
  u1 = Data[,"u"]
  
  # MCMC parameters
  nburn = MCMC$nburn
  nskip = MCMC$nskip
  nsave = MCMC$nsave
  H = MCMC$H # DPM truncation level
  
  # Set hyperparameters and initial values
  # Stick-breaking
  a.alpha=0.25
  b.alpha=0.25
  alpha=rgamma(1,a.alpha,b.alpha)
  V=c(rbeta(H-1,1,10),1)
  w=V*c(1,exp(cumsum(log(1-V[1:(H-1)]))))
  
  
  # Atoms
  a0 = 1.5
  b0 = 0.5
  
  sigma=rinvgamma(H,a0,b0) 
  sigma.hat=1
  B = matrix(0,nrow=H,ncol=p,byrow=T)
  B[1,] = as.vector(lm(u~.-1,data=Data)$coef) 
  f = function(j)
  {
    a.alpha=0.25
    b.alpha=0.25
    alpha=rgamma(1,a.alpha,b.alpha)
    V=c(rbeta(H-1,1,10),1)
    V*c(1,exp(cumsum(log(1-V[1:(H-1)]))))
  }
  ng=round(apply(sapply(1:1e4,f)*n,1,mean))+1
  g = rep(1,H)
  Gamma0 = t(X)%*%X
  
  # Labels
  s=rep(1,n)
  
  # Aux objects
  sigma0 = rep(1,H) 
  B0 = matrix(0,nrow=H,ncol=p,byrow=T)
  
  # Store MCMC samples
  B.output=matrix(-11111,nrow=nsave,ncol=p*H)
  sigma.output=matrix(-11111,nrow=nsave,ncol=H)
  w.output=matrix(-11111,nrow=nsave,ncol=H)
  
  # Run MCMC
  Jsave=0
  for(J in 1:(nburn+nskip*nsave))
  {
    
    tracker = rep(0,5)
    #############################################
    # Step 1. Update V
    nh=table(c(1:H,s))-1
    nl=sapply(1:(H-1),function(j)sum(nh[(j+1):H]))
    V=c(sapply(1:(H-1),function(j)rbeta(1,1+nh[j],alpha+nl[j])),1)
    V[-H][V[-H]==1]=0.9999
    w=V*c(1,exp(cumsum(log(1-V[1:(H-1)]))))
    if(sum(V[-H]==1)>0) break()
    
    #############################################
    # Step 2. Update alpha
    alpha = rgamma(1,a.alpha+H-1,b.alpha-sum(log(1-V[1:(H-1)])))
    if(alpha==0) break()
    tracker[1]=1      
    
    
    #############################################
    # Step 3. Update s
    m=X%*%t(B)
    f.s = function(i)
    {
      probs=log(w)+dnorm(u1[i],m[i,],sqrt(sigma),log=T)
      probs=exp(probs-max(probs))
      sample(1:H,1,prob=probs)
    }
    s=sapply(1:n,f.s)
    tracker[2]=1     
    
    #Step 6. Update sigma and beta
    p1 = 0
    for(j in 1:H)
    {
      if(sum(s==j)>0)
      {
        Gamma0 = (t(X)%*%X)/g[j]
        mu0 = rep(0,p)
        Gamman = (t(X[s==j,,drop=F])%*%X[s==j,,drop=F]+Gamma0)
        mun = ginv(Gamman)%*%t(X[s==j,,drop=F])%*%u1[s==j]
        an = a0 + sum(s==j)/2
        bn = b0 +0.5*(t(u1[s==j])%*%u1[s==j]-t(mun)%*%Gamman%*%mun)
        
        p1 = p1 + (dmvnorm(B0[j,],mu0,as.matrix(forceSymmetric(sigma0[j]*ginv(Gamma0))),log=T)+
                     dinvgamma(sigma0[j],a0,b0,log=T))-
          (dmvnorm(B0[j,],mun,as.matrix(forceSymmetric(sigma0[j]*ginv(Gamman))),log=T)+
             dinvgamma(sigma0[j],an,bn,log=T))
      }
    }
    
    z = rbinom(1,1,0.5/(0.5+0.5*exp(p1)))
    if(z==1)
    {
      B = B0
      sigma = sigma0
    }
    if(z==0)
    {
      for(j in 1:H)
      {
        if(sum(s==j)>0)
        {
          Gamma0 = (t(X)%*%X)/g[j]
          mu0 = rep(0,p)
          Gamman = (t(X[s==j,,drop=F])%*%X[s==j,,drop=F]+Gamma0)
          mun = ginv(Gamman)%*%t(X[s==j,,drop=F])%*%u1[s==j]
          an = a0 + sum(s==j)/2
          bn = b0 +0.5*(t(u1[s==j])%*%u1[s==j]-t(mun)%*%Gamman%*%mun)
          sigma[j] = rinvgamma(1,an,bn)
          B[j,] = rmvnorm(1,mun,as.matrix(forceSymmetric(sigma[j]*ginv(Gamman))))
        }
        if(sum(s==j)==0)
        {
          Gamma0 = (t(X)%*%X)/g[j]
          mu0 = rep(0,p)
          sigma[j] = rinvgamma(1,a0,b0)  
          B[j,] =rmvnorm(1,mu0,as.matrix(forceSymmetric(sigma[j]*ginv(Gamma0))))
        }
      }
    }
    
    # Update g
    if(B[1,1]!=0)
    {
      for(j in 1:H)
      {
        if(sum(s==j)>0) g[j]= 
            rinvgamma(1,(p+1)/2,
                      B[j,]%*%t(X[s==j,])%*%X[s==j,]%*%B[j,]/(2*sigma[j])+(ng[j]/2))
        if(sum(s==j)==0) g[j]=rinvgamma(1,1/2,ng[j]/2)
      }
    }
    for(j in 1:H)
    {
      if(B[1,1]==0) g[j]=rinvgamma(1,(p+1)/2,ng[j]/2)
    }      

    if(J==1)
      print("Burning ...")
    
    if(J<=nburn)
      if(J%%1000==0)
        print(paste("Burning",J,"of",nburn))
#      if(J%%1000==0)
#        print(paste("sigma=",sigma))
    
    if(sum(J==nburn+nskip*(1:nsave))==1)
    {
      Jsave=Jsave+1
      B.output[Jsave,]=as.vector(B)
      sigma.output[Jsave,]=sigma
      w.output[Jsave,]=w
      if(Jsave%%100==0)
        print(paste("Saving",Jsave,"of",nsave))
#      if(Jsave%%100==0)
#        print(paste("sigma=",sigma))
    }
  }
  if(mean(sigma.output[,1]==1)==1)
  {
    output = list(posterior.prob=mean(sigma.output[,1]==1), BF = 999, 
                  message="posterior.prob =1 implies BF > 999")
  }
  BF = 999
  if(mean(sigma.output[,1]==1)<1)
  {
    BF = mean(sigma.output[,1]==1)/(1-mean(sigma.output[,1]==1))
    output = list(posterior.prob=mean(sigma.output[,1]==1), BF=BF)
  }
  output
}
