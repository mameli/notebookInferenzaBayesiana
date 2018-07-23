#Funzione per campionare da una distribuzione normale multivariata
rmvnorm <- function(n, mu, Sigma)
{
  #samples form the multivariate normal distribution
  E<-matrix(rnorm(n*length(mu)), n, length(mu))
  t( t(E%*%chol(Sigma)) +c(mu))
}
#Lettura dei dati :
dati<-as.matrix(dati<-read.table(
  "/home/mameli/git/notebookInferenzaBayesiana/code/msparrownest.dat", 
  col.names=c("Y" , "X")))
head(dati)

#Matrice del modello e vettore delle osservazioni :
X = cbind(1, dati[,2])
head(X)

y = dati[,1]
head(y)

n<-length(y)
p<-dim(X)[2]


#b) inizializziamo la prior:

pmn.beta<-c(1.1, 0)
psd.beta<-c(2.75, 0.22)

library(coda)
metropolis<-function(tuning, nsimul) {	
  pmn.beta<-c(1.1, 0)
  psd.beta<-c(2.75, 0.22)
  var.prop<-tuning
  beta<-rep(0, p)
  S<-nsimul
  BETA<-matrix(0, nrow=S, ncol=p)
  ac<-0
  set.seed(1)
  for(s in 1:S){
    beta.p<-t(rmvnorm(1, beta, var.prop))
    lhr<- sum(log(dbinom(y,1,exp(X%*%beta.p)/(1+exp(X%*%beta.p))))) + dnorm(beta.p[1],pmn.beta[1],psd.beta[1],log=TRUE) + dnorm(beta.p[2],pmn.beta[2],psd.beta[2],log=TRUE)-sum(log(dbinom(y,1,exp(X%*%beta)/(1+exp(X%*%beta)))))-dnorm(beta[1],pmn.beta[1],psd.beta[1],log=TRUE)-dnorm(beta[2],pmn.beta[2],psd.beta[2],log=TRUE)
    if(log(runif(1))<lhr) beta<-beta.p; ac<-ac+1
    BETA[s,]<-beta 
  }
  Ef<-effectiveSize(BETA)
  
  cat("acceptance rate=", ac/S, "\n")
  cat("effective sample size=", Ef, "\n")
  return (BETA)
}

nsimul <- 1000
var.prop <- var(log(y + 1 / 2))*solve(t(X) %*% X)
var.prop

beta.post <- metropolis(var.prop, nsimul)

var.prop <- var.prop*9
nsimul <- 5000
beta.post <- metropolis( var.prop, nsimul)

skips<-seq(5,nsimul,by=10)
plot(skips ,beta.post[skips ,1], type="l", xlab="iteration",ylab=expression(alpha), col = "deeppink")
plot(skips, beta.post[skips ,2] ,type="l",xlab="iteration", ylab=expression (beta), col = "deepskyblue")

par(mfrow=c(1,2))
x<-seq(-10, 10, by=0.1)
plot(x, dnorm(x, pmn.beta[1], psd.beta[1]) ,ylim=c(0,0.25), type="l",lwd=1, xlab=expression(alpha), ylab="density", col = "deepskyblue")
lines(density(beta.post[ ,1]) , col = "deeppink")
legend("topright",legend=c("Prior","Posterior"),lty = c(1,1),cex=0.7, bty="n",seg.len=1.5,  col = c("deepskyblue","deeppink"))
y<-seq(-2,2,by=0.01)
plot(y,dnorm(y,pmn.beta[2],psd.beta[2]), ylim=c(0,3), type="l", lwd=1, xlab=expression(beta), ylab="density", col = "deepskyblue")
lines(density(beta.post[ ,2]) , col="deeppink")
legend("topright",legend=c("Prior","Posterior"),lty = c(1,1),cex=0.7,bty="n", seg.len=1.5, col = c("deepskyblue","deeppink"))

f<-exp(t(X%*%t(beta.post)))/(1+exp(t(X%*%t(beta.post))))
qE<-apply(f,2,quantile,probs=c(0.025,0.975))
par(mfrow=c(1,1))
plot(c(10,15),range(c(0,qE)),type="n",xlab="wingspan",ylab="f")
lines(qE[1,],col="deepskyblue",lwd=2)
lines(qE[2,],col="deeppink",lwd=2)