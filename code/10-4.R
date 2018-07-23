n<-length(y)
p<-dim(X)[2]


#b) inizializziamo la prior:

pmn.beta<-c(1.1, 0)
psd.beta<-c(2.75, 0.22)

metropolis<-function(tuning, nsimul) {	
  pmn.beta<-c(1.1, 0)
  psd.beta<-c(2.75, 0.22)
  var.prop<-tuning
  beta<-rep(0, p)
  S<-nsimul
  BETA<-matrix(0, nrow=S, ncol=p)
  ac<-0
  set.seed(1)
  library(coda)
  for(s in 1:S){
    beta.p<-t(rmvnorm(1, beta, var.prop))
    lhr<- sum(log(dbinom(y,1,exp(X%*%beta.p)/(1+exp(X%*%beta.p))))) + dnorm(beta.p[1],pmn.beta[1],psd.beta[1],log=TRUE) + dnorm(beta.p[2],pmn.beta[2],psd.beta[2],log=TRUE)-sum(log(dbinom(y,1,exp(X%*%beta)/(1+exp(X%*%beta)))))-dnorm(beta[1],pmn.beta[1],psd.beta[1],log=TRUE)-dnorm(beta[2],pmn.beta[2],psd.beta[2],log=TRUE)
  }
  if(log(runif(1))<lhr){beta<-beta.p; ac<-ac+1}
  BETA[s,]<-beta 
  Ef<-effectiveSize(BETA)
  
  cat("acceptance rate=", ac/S, "\n")
  cat("effective sample size=", Ef, "\n")
  return (BETA)
}

nsimul <- 1000
var.prop <- var(log(y + 1 / 2)) * solve(t(X) %*% X)
var.prop