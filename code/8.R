y.barA<-75.2
y.barB<-77.5
sA<-7.3
sB<-8.1
nA<-nB<-n_<-16

delta0<-c(-4 ,-2 ,0 ,2 ,4)
tau20<-c(10, 50, 100, 500)
mu0<-75
gamma20<-100
v0<-2
sigma20<-100

delta <-(y.barA-y.barB)/2
mu<-(y.barA+y.barB)/2

nsimul<-1000


gibbs<-array(NA,c(nsimul,3,length(delta0)*length(tau20)))

#Ciclo Gibbs:

v<-1
for(j in delta0){
  for(k in tau20){
    for(i in 1:nsimul) {
      
      #Aggiorniamo sigma:
      
      vn<-v0+nA+nB
      vnsigma2n<-v0*sigma20+(n_-1)*(sA^2+sB^2)+n_*(y.barA^2+y.barB^2)+2*n_*(mu^2+delta^2-mu*(y.barA+y.barB)+delta*(y.barB-y.barA))
      sigma2<-1/rgamma(1,vn /2,vnsigma2n/2)
      
      #Aggiorniamo mu:
      
      gamma2n<-gamma20*sigma2/(sigma2+(nA+nB)*gamma20)
      mun<-gamma2n*(mu0/gamma20+n_*(y.barA+y.barB)/sigma2)
      mu<-rnorm(1,mun,sqrt(gamma2n))
      
      #Aggiorniamo delta:
      
      tau2n<-k*sigma2 /(sigma2+(nA+nB)*k)
      deltan<-tau2n*(j/k+n_*(y.barA-y.barB)/sigma2)
      delta<-rnorm(1,deltan,sqrt(tau2n))
      
      gibbs[i,,v] <- c(mu,delta,sigma2)
    }
    v<-v+1
  }
}

colnames(gibbs)<-c("mu","delta","sigma2")

#a )

#Probabilita' a posteriori che la semi-differenza tra le medie sia negativa
#+ intervallo di credibilita' a posteriori per la semi-differenza tra le
#medie + correlazione a priori e a posteriori tra la media del primo gruppo
#e quella del secondo (tutto per ognuna delle possibili prior definite su
#delta):
correlazione.post<-matrix (NA,5,4)
probabilita.post<-corrrelazione.post
quantili.post<-NULL
v=1
for(i in 1:5){
  for(j in 1:4){
    probabilita.post[i,j]<-mean(gibbs[,2,v]<0)
    quantili.post<-rbind(quantili.post,quantile(gibbs[,2,v],c(0.025,0.975)))
    correlazione.post[i,j]<-cor(gibbs[,1,v]+gibbs[,2,v],
                                gibbs[,1,v]-gibbs[,2,v])
    v<-v+1
  }
}

rownames(probabilita.post)<-rownames(correlazione.post)<-
  c("delta0=-4","delta0=-2","delta0=0","delta0=2","delta0=4")
colnames(probabilita.post)<-colnames(correlazione.post)<-
  c("tau20=10"," tau20=50"," tau20=100"," tau20=500")
rownames(quantili.post)<-c("delta0=-4 tau20=10"," delta0=-4 tau20=50",
                           "delta0=-4 tau20=100","delta0=-4 tau20=500","delta0=-2 tau20=10",
                           "delta0=-2 tau20=50","delta0=-2 tau20=100","delta0=-2 tau20=500",
                           "delta0=0 tau20=10","delta0=0 tau20=50","delta 0=0 tau20=100",
                           "delta0=0 tau20=500","delta0=2 tau20=10","delta0=2 tau20=50",
                           "delta0=2 tau20=100","delta0=2 tau20=500","delta0=4 tau20=10",
                           "delta0=4 tau20=50","delta0=4 tau20=100","delta0=4 tau20=500")

#i)
probabilita.post

par(mfrow=c(3,2))
for(i in 1:5){
  plot(tau20,probabilita.post[i,],pch=20,
       xlab="tau^2[0]",
       ylab="p(delta<0 given Y)", type="l")
       legend("right", paste("delta0=",delta0[i]), bty="n")
}

quantili.post

correlazione.prior<-c(0.81,0.33,0,-0.67)

correlazione.post

win.graph()
par(mfrow=c(2,2))
x<-seq(-50,50,by=0.01)
v<-0
for(j in 1:4){
  plot(x,dnorm(x,delta0[1],sqrt(tau20[j])), xlim=c(-10,10),ylim=c(0,0.3),
       xlab=expression(delta),ylab="density",type="l",col="grey")
  lines(density(gibbs[,2,j]))
  legend("topright",legend=c("posterior","prior"),lwd=c(2,2),col=c("black",
                                                                   "gray"),bty="n")
  text(5.5,0.15,paste("tau20=",tau20[j]))
}

win.graph()
par(mfrow=c(2,2))
x<-seq(-50,50,by=0.01)
v<-0
for(j in 1:4){
  plot(x,dnorm(x,delta0[5],sqrt(tau20[j])),xlim=c(-10 ,10),ylim=c(0,0.3),
       xlab=expression(delta),ylab="density",type="l",col="grey")
  lines(density(gibbs[,2,j+16]))
  legend("topright",legend=c("posterior","prior"),lwd=c(2,2),col=c("black",
                                                                   "gray"),bty="n")
  text(5.5,0.15, paste("tau20=",tau20[j]))
}

library(coda)
effectiveSize(gibbs[,,1])