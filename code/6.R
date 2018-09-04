A<-c(1, 0, 0, 1, 2, 2, 1, 5, 2, 0, 0, 0, 0, 0,
     0, 1, 1, 1, 0, 0, 0, 1, 1, 2, 1, 3, 2, 0, 0,
     3, 0, 0, 0, 2, 1, 0, 2, 1, 0, 0, 1, 3, 0, 1,
     1, 0, 2, 0, 0, 2, 2, 1, 3, 0, 0, 0, 1, 1)

B<-c(2, 2, 1, 1, 2, 2, 1, 2, 1, 0, 2, 1, 1, 2,
     0, 2, 2, 0, 2, 1, 0, 0, 3, 6, 1, 6, 4, 0, 3, 2, 0, 1,
     0, 0, 0, 3, 0, 0, 0, 0, 0, 1, 0, 4, 2, 1, 0, 0, 1, 0,
     3, 2, 5, 0, 1, 1, 2, 1, 2, 1, 2, 0, 0, 0, 2, 1,
     0, 2, 0, 2, 4, 1, 1, 1, 2, 0, 1, 1, 1, 1, 0, 2, 3, 2,
     0, 2, 1, 3, 1, 3, 2, 2, 3, 2, 0, 0, 0, 1, 0, 0,
     0, 1, 2, 0, 3, 3, 0, 1, 2, 2, 2, 0, 6, 0, 0, 0, 2, 0,
     1, 1, 1, 3, 3, 2, 1, 1, 0, 1, 0, 0, 2, 0, 2, 0,
     1, 0, 2, 0, 0, 2, 2, 4, 1, 2, 3, 2, 0, 0, 0, 1, 0, 0, 1,
     5, 2, 1, 3, 2, 0, 2, 1, 1, 3, 0, 5, 0, 0, 2,
     4, 3, 4, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 2, 0, 0, 1, 1, 0,
     2, 1, 3, 3, 2, 2, 0, 0, 2, 3, 2, 4, 3, 3, 4,
     0, 3, 0, 1, 0, 1, 2, 3, 4, 1, 2, 6, 2, 1, 2, 2)

a_theta<-2
b_theta<-1
v_gamma<-c(8, 16, 32, 64, 128)

n<-5000
set.seed(150)

valori_medie<-NULL

gibbs<-array(NA, c(n, 2, 5))

nA<-length(A)
nB<-length(B)

ytotA<-sum(A)
ytotB<-sum(B)
ytot<-ytotA+ytotB

gamma<-1
theta<-1
k<-1

for(i in v_gamma){
  a_gamma<-b_gamma<-i
  for(j in 1:n){
    theta<-rgamma(1, a_theta+ytot, b_theta+nA+nB*gamma)
    gamma<-rgamma(1, a_gamma+ytotB, b_gamma+nB*theta)
    gibbs[j,,k]<-cbind(theta, gamma)
  }
  k<-k+1
}

for (i in 1:5){
  valori_medie<-c(valori_medie, mean(gibbs[,1,i]*gibbs[,2,i]
                                     -gibbs[,1,i]))
}
valori_medie

leg<-v_gamma
leg.txt<-rep("ab=bb",5)

par(mfrow=c(3,2))
for(i in 1:5){
  hist(gibbs[,1,i],prob=T,col="deeppink",ylim=c(0,5),xlim=c(0,2.5),
       ylab="densita' a posteriori",xlab="numero medio di figli",main="")
lines(density(gibbs[,1,i]))
hist(gibbs[,1,i]*gibbs[,2,i],prob=T,col="deepskyblue",add=T)
lines(density(gibbs[,1,i]*gibbs[,2,i]))
}
for(i in 1:5){
hist(gibbs[,1,i],prob=T,col="deeppink",ylim=c(0,5),xlim=c(0,2.5),
ylab="densita' a posteriori",xlab="numero medio di figli",main="")
lines(density(gibbs[,1,i]*gibbs[,2,i]))
hist(gibbs[,1,i]*gibbs[,2,i],prob=T,col="deepskyblue",add=T)
lines(density(gibbs[,1,i]*gibbs[,2,i]))
text(2, 4.5, paste("a_gamma=",leg[i]))
text(2, 3.5, paste("b_gamma=",leg[i]))
}

x<-seq(0, 10, by=0.01)
par(mfrow=c(1,1))
plot(x, dgamma(x,8,8),type="l",xlim=c(0,5),ylim=c(0,5),
     xlab=expression(gamma),ylab=expression(p(gamma)),
     main="gamma a priori",col=1)


for(j in 2:5){
  curve(dgamma(x,v_gamma[j],v_gamma[j]),add=T,col=j)
}
legend(2,4,c("a_gamma=b_gamma=8","a_gamma=b_gamma=16",
             "a_gamma=b_gamma=32", "a_gamma=b_gamma=64",
             "a_gamma=b_gamma=128"),col=c(1,2,3,4,5),lty=1)