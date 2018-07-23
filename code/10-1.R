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