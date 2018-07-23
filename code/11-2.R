rmvnorm <- function(n, mu, Sigma)
{
  E <- matrix(rnorm(n * length(mu)), n, length(mu))
  t(t(E %*% chol(Sigma)) + c(mu))
}

plot.hdr2d <- function(x,
                       prob = c(.025, .25, .5, .75, .975),
                       bw = c(5, 5),
                       cols = gray(((length(prob) - 1):1) / length(prob)),
                       xlim = range(x[, 1]),
                       ylim = range(x[, 2]),
                       ...)
{
  plot(c(0, 0),
       xlim = xlim,
       ylim = ylim,
       type = "n",
       ...)
  add.hdr2d(x, prob, bw, cols)
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

dati <-read.table("~/git/notebookInferenzaBayesiana/code/pdensity.dat", header = TRUE, sep = "")
head(dati)
gruppo <- unique(dati$plot)
m <- length(gruppo)
Y <- list()
X <- list()
N <- NULL
for (j in 1:m)
{
  Y[[j]] <-
    dati[dati[, 1] == gruppo[j], 3]
  N[j] <- sum(dati$plot == gruppo[j])
  xj <- dati[dati[, 1] == gruppo[j], 2]
  X[[j]] <- cbind(rep(1, N[j]), xj, xj ^ 2)
}
###stima dei singoli coefficienti ols####
S2.LS <- BETA.LS <- NULL
for (j in 1:m) {
  fit <- lm(Y[[j]] ~ -1 + X[[j]])
  BETA.LS <- rbind(BETA.LS, c(fit$coef))
  S2.LS <- c(S2.LS, summary(fit)$sigma ^ 2)
}
colnames(BETA.LS) <- c("1", "x", "x^2")
BETA.LS

####
par(mfrow = c(1, 2))
plot(
  range(dati [, 2]) ,
  range (dati [, 3]) ,
  type = "n" ,
  xlab = "planting density ",
  ylab = "yield "
)
for (j in 1:m) {
  curve(BETA.LS[j, 1] + BETA.LS[j, 2] * x + BETA.LS[j, 3] * x ^ 2, col = "gray",add = T)
}
BETA.LS.MEAN <- apply(BETA.LS, 2 , mean)
curve (BETA.LS.MEAN[1] + BETA.LS.MEAN[2] * x + BETA.LS.MEAN[3] * x ^ 2, lwd = 2, add = T)

##### hierarchical regression model
p <- dim(X[[1]])[2]
theta <-mu0 <- apply(BETA.LS, 2, mean) #come valori iniziali di theta e muo utilizziamo le stime ols
nu0 <- 1
#parametro che entra nella a priori delle sigma^2j
s20 <- mean(S2.LS) # parametro che entra nella a priori delle sigma^2j
eta0 <- p + 2
#parametro per inverse Wishart per Sigma
Sigma <-S0 <- L0 <- cov(BETA.LS)
# Prior variance of the Betas (\Lambda_0) and
#S0 set to the empirical covariance of the beta_OLS
#(matrice di covarianza tra i coefficienti di regressione totale)
BETA <- BETA.LS #assegna ai beta le stime ols
S2 <- S2.LS #è il vettore delle sigma^2 j
THETA.ps <- NULL
NU0 <- NULL
#matrici che fa comodo calcolare fuori dal ciclo.
iL0 <- solve(L0)
#matrice di var e cov dei theta
iSigma <- solve(Sigma)
#valore iniziale per la varianza dei beta
ACR <- NULL
Sigma.ps <- matrix(0, p, p)
SIGMA.PS <- NULL
BETA.ps <- NULL
S20.ps <- NULL
NU0.ps <- NULL
S2.ps <- NULL
set.seed(1)
acs <- 0
#L0 è la matrice che entra nella a priori di theta
#s0 è la matrice che entra nella a priori di SIGMA
#SIGMA è la matrice di var e cov che entra nella a priori dei beta\\

for (s in 1:10000) {
  ##update beta_j
  for (j in 1:m)
  {
    Vj <- solve(iSigma + t(X[[j]]) %*% X[[j]] / S2[j])
    Ej <- Vj %*% (iSigma %*% theta + t(X[[j]]) %*% Y[[j]] / S2[j])
    BETA[j, ] <- rmvnorm(1, Ej, Vj)
  }
  ##
  ##update theta
  Lm <- solve(iL0 + m * iSigma)
  mum <- Lm %*% (iL0 %*% mu0 + iSigma %*% apply(BETA, 2, sum))
  theta <- t(rmvnorm(1, mum, Lm))
  ##
  ##update Sigma
  mtheta <- matrix(theta, m, p, byrow = TRUE)
  iSigma <-rwish(1, eta0 + m, solve(S0 + t(BETA - mtheta) %*% (BETA - mtheta)))
  ##
  ##update s2j
  RSS <- 0
  for (j in 1:m) {
    RSS <- sum((Y[[j]] - X[[j]] %*% BETA[j, ]) ^ 2)
    S2[j] <- 1 / rgamma(1, (nu0 + N[[j]]) / 2, (nu0 * s20 + RSS) / 2)
  }
  ##
  #update sigma0
  sumsig = 0
  for (j in 1:m) {
    sumsig <- sumsig + (1 / S2[j])
  }
  s20 <- rgamma(1, (m * nu0 / 2) + 2, (nu0 / 2 * sumsig) + 2)
  #update di nu0 bisogna fare un passo di metropolis
  nu0.star <- sample(max(1, nu0 - 2):min(nu0 + 2, 100), 1)
  #questa è una proposal non simmetrica, che dipende dal valore campionato precedentemente
  log.star = 0
  log.s = 0
  for (j in 1:m) {
    log.star = log.star + dgamma(1 / S2[j], nu0.star, s20, log = TRUE)
    log.s = log.s + dgamma(1 / S2[j], nu0, s20, log = TRUE)
  }
  log.r = log.star - log.s + log(1 / (min(nu0.star + 2, 100) + 1
                                      - max(1, nu0.star - 2)) - log(1 / (min(nu0 +
                                                                               2, 100) + 1 - max(1, nu0 - 2))))
  if (log(runif(1)) < log.r) {
    nu0 <- nu0.star
    acs <- acs + 1
  }
  NU0 <- c(NU0, nu0)
  S20.ps<-c(S20.ps,s20);
  THETA.ps<-rbind(THETA.ps,t(theta))
  #Sigma.ps<-Sigma.ps+solve(iSigma) ;
  BETA.ps<-rbind(BETA.ps,t(BETA))
  SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
  S2.ps<-rbind(S2.ps,t(S2))
}

library(coda)
effectiveSize(S20.ps) #\sigma0^2
