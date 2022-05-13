#### To compute empirical power i.e. when beta_2=0####

#### For ZINB  Distribution
#####################################################

library(MASS)
library(VGAM)

t=c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50)
for(p in 1:11)
{
  n=20
  k=100
  lrt=c()
  
  lcount=0 
  
  beta0=1
  beta1=0.5
  #beta2=t[p]
  beta2=0
  
  for(i in 1:k){
    
    y=c();x1=c();x2=c();m=c();qf=c()
    
    for(j in 1:n){
      x1[j]=rnorm(1,0,1)
      x2[j]=rnorm(1,0,1)
      m[j]=exp(beta0+beta1*x1[j]+beta2*x2[j])
      y[j] <- rzinegbin(1,pstr0=0.3,mu=m[j], size=30)
    }
    
    ############Likelihood Ratio Statistics#####
    fit=vglm(y~x1+x2, zinegbinomial)
    fit0=vglm(y~x1, zinegbinomial)
    lrt[i]=2*(logLik(fit)-logLik(fit0))
    
    if(lrt[i]>3.84) {lcount=lcount+1}
    
  }
  write.table(t(c(lcount/100)),file="ZINB Empirical Level_n=20.txt",
              append=TRUE,row.names=F,col.names=F)  
}


