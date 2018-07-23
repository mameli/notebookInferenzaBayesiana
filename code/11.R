rmvnorm<-function(n,mu, Sigma ) {
  p<-length(mu)
  res<-matrix ( 0, nrow=n, ncol=p)
  if ( n>0 & p>0 ) {
    E<-matrix( rnorm(n*p),n, p)
    res<-t ( t (E%*%chol ( Sigma ) ) +c (mu) )
    #R <- chol(A)
  }
  res
}


rinvwish <- function(n, nu0, iS0)
{
  sL0 <- chol(iS0)
  S <- array(dim = c (dim(iS0), n))
  for (i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(iS0)[1]), nu0, dim(iS0)[1]) %*% sL0
    S[, , i] <- solve(t(Z) %*% Z)
  }
  S[, , 1:n]
}
rwish <- function(n, nu0, S0)
{
  sS0 <- chol(S0)
  S <- array(dim = c(dim(S0), n))
  for (i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[, , i] <- t(Z) %*% Z
  }
  S[, , 1:n]
}

dati<-read.table("~/git/notebookInferenzaBayesiana/code/pdensity.dat", header=TRUE)
head (dati)

ids<-unique (dati$plot )
m<-length ( ids )
Y<-list() ; X<-list() ; N<-NULL
for ( j in 1:m){
  Y[[ j ]]<-dati [ dati [,1]== ids [ j ], 3]
  N[j]<- sum( dati$plot==ids [ j ])
  xj<-dati [ dati [,1]== ids [ j ], 2]
  X[[ j ]]<-cbind ( rep (1,N[ j ]), xj, xj^2 )
}
p<-dim(X[[1]])[2]
N


Y[[1]]

S2.OLS<-BETA.OLS<-NULL
for( j in 1:m){
  fit<-lm(Y[[j]]~-1+X[[j]])
  BETA.OLS<-rbind (BETA.OLS, c(fit$coef ) )
  S2.OLS<-c ( S2.OLS, summary( fit )$sigma ^2)
}

colnames (BETA.OLS)<-c ("1"," x "," x^2")
BETA.OLS

par (mfrow=c (1,2) )
plot ( range ( dati [,2]), range ( dati [,3]), type="n", xlab="planting density ", ylab="Expected yield ",main="Regression lines OLS")
for (j in 1:m) {
  curve (BETA.OLS[j, 1] + BETA.OLS[j, 2] * x + BETA.OLS[j, 3] * x ^ 2, col = "deepskyblue", lty = 2, add = T)
}
BETA.OLS.MEAN <- apply (BETA.OLS, 2, mean)
curve (BETA.OLS.MEAN[1] + BETA.OLS.MEAN[2] * x + BETA.OLS.MEAN[3] * x ^ 2,lwd = 2, col = "deeppink", add = T)

BETA.OLS.MEAN


theta<-BETA.OLS.MEAN; Sigma<-cov(BETA.OLS) ; sigma2<-mean(S2.OLS)
eta0 <-4; S0<-Sigma
mu0<-theta ; L0<-Sigma
v0<-2; sigma20<-sigma2

#c)

nsimul=10000

beta<-BETA.OLS
Sigma<-cov(BETA.OLS)
sigma2<-mean(S2.OLS)

Sigma.post<-matrix (0,p,p)
BETA.pp<-THETA.POST<-S2.POST<-NULL
BETA.POST<-BETA.OLS*0
SIGMA.POST<-array (0, c(p,p, nsimul ) )
set.seed (1)

for ( s in 1: nsimul ) {
  for ( j in 1:m){
    Vj<-solve ( solve (Sigma) + t(X[[j]])%*%X[[j]]/sigma2 )
    Ej<-Vj%*%( solve (Sigma)%*%theta + t(X[[j]])%*%Y[[j]]/sigma2 )
    beta [ j,]<-rmvnorm(1,Ej, Vj)
  }
  
  #theta 
  Lm<- solve ( solve (L0) + m* solve (Sigma) )
  mum<- Lm%*%( solve (L0)%*%mu0 + solve (Sigma)%*%apply ( beta,2,sum) )
  theta<-t (rmvnorm(1,mum,Lm) )
  
  #Sigma
  mtheta<-matrix ( theta,m,p, byrow=TRUE)
  Sigma<-solve ( rwish (1, eta0+m, solve ( S0+t ( beta-mtheta)%*%(beta-mtheta) ) ) )
  #Campioniamo la varianza residua:
  
  RSS<-0
  for ( j in 1:m) {
    RSS<-RSS+sum((Y[[j]]-X[[j]]%*%beta[j,] )^2 )
  }
  sigma2<-1/rgamma(1,( v0+sum(N) ) /2, (v0*sigma20+RSS)/2 )
  #Immagazziziniamo i valori appena campionati:
  S2.POST<-c(S2.POST, sigma2 ) ;THETA.POST<-rbind (THETA.POST, t ( theta ) )
  Sigma.post<-Sigma.post+Sigma ; BETA.POST<-BETA.POST+beta
  SIGMA.POST[,, s]<-Sigma
  #Campioniamo dalla posterior pdeeppinkictive dei coefficienti di regressione che ci servira' per il punto e).
  BETA.pp<-rbind (BETA.pp,rmvnorm(1, theta, Sigma) )
}

colnames (THETA.POST)<-c(" theta1 "," theta2 "," theta3 ")
colnames (BETA.POST)<-colnames (BETA.pp)<-c(" beta1 "," beta2 "," beta3 ")
BETA.PM<-BETA.POST/nsimul
plot (range (c (0, 10)), range (c (0, 10)), type = "n", xlab = "planting density ", ylab = "Yield ", main = "Bayesian regression lines ")
for (j in 1:m) {
  curve (BETA.PM[j, 1] + BETA.PM[j, 2] * x + BETA.PM[j, 3] * x ^ 2, lty = 2, col = "deepskyblue ", add = T)
}
curve (mean(THETA.POST[, 1]) + mean(THETA.POST[, 2]) * x + mean(THETA.POST[, 3]) * x ^ 2, lwd = 2, add = T, col = "deeppink")

effectiveSize (THETA.POST)

n<-10000

THETA.PRIOR<-rmvnorm(n,mu0,L0) ; colnames (THETA.PRIOR)<-colnames (THETA.POST)
SIGMA.PRIOR<-rinvwish(n, eta0, solve( S0 ) )
head(THETA.PRIOR)

head(SIGMA.PRIOR)

#A priori e a posteriori di theta a confronto ( plot ):

par (mfrow = c(2, 2))
for (i in 1:3) {
  plot (density (THETA.POST[, i]), xlab = paste (" theta ", i), main = "", col = "deeppink")
  lines (density (THETA.PRIOR[, i]), col = "deepskyblue ")
  legend ("topright", legend = c(" Prior ", " Posterior "), col = c(" deepskyblue ", " deeppink "), lty = 1, bty = "n", cex = 0.7)
}

par (mfrow = c(2, 2))
for (i in 1:3) {
  plot (density (log(SIGMA.POST[i, i,])), type = "l", xlab = paste ("Sigma", i), col = "deeppink", ylab = "density ", main = "")
  lines (density (log(SIGMA.PRIOR[i, i,])), col = "deepskyblue ", lwd = 2)
  legend("topright", legend = c(" Prior ", " Posterior "), col = c(" deepskyblue ", " deeppink "), bty = "n", lty = 1)
}

x0<-c (2,4,6,8)
raccolto<-NULL
for ( i in x0){
  x<-c (1, i, i ^2)
  raccolto<-cbind ( raccolto,BETA.pp%*%x)
}
colnames ( raccolto )<-x0

#Confrontiamo adesso le quattro distribuzioni dei valori attesi,
#graficamente e con un indice sintetico, la media:

par (mfrow = c (1, 1))
hist (raccolto [, 1], prob = T, main = "Raccolto atteso a posteriori ", xlab = "Raccolto ", xlim = c (2, 13), ylim = c (0, 1.5), col = "deepskyblue")
lines (density (raccolto [, 1]), col = "deepskyblue")
hist (raccolto [, 2], prob = T, col = "deeppink ", add = T)
lines (density (raccolto [, 2]), col = "deeppink")
hist (raccolto [, 3], prob = T, col = " chartreuse ", add = T)
lines (density (raccolto [, 3]), col = " chartreuse ")
hist (raccolto [, 4], prob = T, col = " darkorange ", add = T)
lines (density (raccolto [, 4]), col = "darkorange ")
legend ("topright",legend = c("x=2", "x=4", "x=6", "x=8"),col = c("deepskyblue", " deeppink ", "chartreuse", " lightsalmon "), lty = 1)
apply (raccolto, 2, mean)

apply ( raccolto,2,mean)

y.pred<-rnorm( nsimul,BETA.pp%*%x, S2.POST)

quantile(y.pred, c (0.025,0.975)) 
          
hist(y.pred, prob=T, main="Predittiva a posteriori ", xlab="raccolto ")
lines (y.pred )
