#Numero di osservazioni e numero dei parametri si cui fare inferenza:
n<-length(y)
p<-dim(X)[2]
#b)
#Si veda il setting per la formulazione della prior:
pmn.beta<-c(1.1, 0)
psd.beta<-c(2.75, 0.22)
#c)
#Creiamo una funzione che approssima la distribuzione a posteriori dei
#coefficienti di regressione mediante algoritmo di Metropolis-Hastings . E'
#richiesto di aggiustare la distribuzione proposta in modo che il tasso di
#accettazione sia ragionevole e di considerare un numero di iterazioni che
#portano ad una effective sample size di circa 1000: pertanto la funzione
#creata prende in ingresso questi due parametri 
#Osservazione :
#abbiamo appena detto che la matrice di varianza e covarianza della
#distribuzione proposal e ' inserita come input della funzione in modo da
#vedere come cambiano i risultati al variare di questa e poterla cosi '
#scegliere in maniera adeguata a rispondere alle rischieste . Ricordiamo
#pero ' che essa e ' fissa per ogni catena , ovvero ogni catena ha la sua . In
#verita ' in alcuni casi si puo ' anche aggiustare durante l ' algoritmo a
#patto che siano soddisfatte certe condizioni : in teoria si puo' quindi
#estrarre la varianza ogni volta ma a patto che la distribuzione da cui
#viene estratta non dipenda dai valori dei parametri estratti nella catena
#(a parte quelli dell 'iterazione precedente).
#Questo succede in situazioni piu' complicate in cui a volte vorremmo poter > #fare piccoli passi e a volte grandi : quindi e ' possibile fare un certo
#numero di iterazioni con varianza piccola e un altro con varianza grande , > #sempre sotto la condizione che queste due varianze siano prespecificate o > #comunque prese random e non dipendere dai valori estratti durante la >#catena). In questo modo la distribuzione proposal e' piu' flessibile e
#riesce ad esplorare piu ' facilmente la distribuzione a posteriori . Come >#esempio concreto si puo' considerare il caso in cui la distribuzione a
#posteriori di interesse e ' bimodale e quindi per approssimarla campionando
#secondo una algoritmo di Metropolis c'e' bisogno a volte di grandi apssi
#per passare da una moda all'altra, altre di piccoli passi in modo che il
#tasso di accettazione non sia basso. Si aumenta cosi' la cosiddetta
#capacita' di mixing dell'algoritmo.
metropolis<-function(tuning, nsimul) {
	#setting della distribuzione a priori:
	pmn.beta<-c(1.1, 0)
	psd.beta<-c(2.75, 0.22)
	#setting della distribuzione proposal: si sceglie come spesso si fa una
	#distribuzione normale multivariata con media vettore nullo e matrice
	#di varianza e covarianza quella specificata in ingresso nella funzione.
	var.prop<-tuning
	#valore inizial del vettore dei coefficenti:
	beta<-rep(0, p)
	#numero di simulazioni
	S<-nsimul
	#Vettore in cui immagazzino i valori campionati della distribuzione a posteriori:
	BETA<-matrix(0, nrow=S, ncol=p)
	#contatore del numero di accetazioni:
	ac<-0
	set.seed(1)
	library(coda)
	#algoritmo metropolis (dal momento che la proposal e' simmetrica siamo in
	#questo caso pasrticolare di Metropolis-Hastings):
	 for(s in 1:S){
	 	#proposta dei coefficenti di regressione:
	 	beta.p<-t(rmvnorm(1, beta, var.prop))
	 	#Rapporto di metropolis:
	 	#puo' essere calcolato in piu' modi tra cui usare la forma funzionale
	 	#della verosomiglianza trovata al punto a (e' scritto per completezza in
	 	#commento) o usando la funzione dbinom di R. Facciamo relativamente a
	 	#questa alcune considerazioni:
	 	#-) la funzione prende in ingresso tre elementi: il vettore y delle
	 	#osservazioni, la dimensione n e il vettore della probabilita' p. In
	 	#questo caso, dal momento che la dimensione specificata e' 1, la funzione
	 	#restituisce un vettore della stessa dimensione del vettore risposta il cui
	 	#i-mo elemento e' il valore della densita' Bernoulliana calcolata
	 	#nell'elemento i-mo di y con probabiita' pari all'i-mo elemento del
	 	#vettore delle probabilita' p. La somma dei logaritmi di tali densita' e' la
	 	#log- verosomiglianza.
	 	#-) Il vettore di probabilita' specificato nella funzione dbinom e'
	 	#calcolato secondo la formulazione del modello logistico come xpit del
	 	#predittore lineare. Il rapporto di metropolis contiene due
	 	#log-verosomiglianze: una calcolata considerando come vettore di 
	 	#probababilita' quello calcolato con il vettore dei coefficenti
	 	#all'iterazione corrente e una considerando quello calcolato
	 	#all'iterazione precedente.
	 	#-) Il rapporto di metropolis, oltre al rapporto tra verosomiglianze,
	 	#contiena nche il rapporto tra la densita' della prior specificata per
	 	#i coefficenti di regressione valutata all'iterazione corrente e quella
	 	#valutata all'iterazione precedente. Usiamo l'opzione log= TRUE perche'
	 	#lavoriamo con log-verosomiglizne.
	 	#-) Nulla viete di lavorare con le verosogiglianze, ma usiamo le
	 	#log-verosogmilianze perche' migliori da pun punto di vista
	 	#computazionale.
	 	lhr<- sum(log(dbinom(y,1,exp(X%*%beta.p)/(1+exp(X%*%beta.p))))) + dnorm(beta.p[1],pmn.beta[1],psd.beta[1],log=TRUE) + dnorm(beta.p[2],pmn.beta[2],psd.beta[2],log=TRUE)-sum(log(dbinom(y,1,exp(X%*%beta)/(1+exp(X%*%beta)))))-dnorm(beta[1],pmn.beta[1],psd.beta[1],log=TRUE)-dnorm(beta[2],pmn.beta[2],psd.beta[2],log=TRUE)
		}
		if(log(runif(1))<lhr){beta<-beta.p; ac<-ac+1}BETA[s,]<-beta 
		Ef<-effectiveSize(BETA)
		#Tasso di accettazione ed effective sample size:
		cat("acceptance rate=", ac/S, "\n")
		cat("effective sample size=", Ef, "\n")
		return (BETA)
	}
	#cerchiamo adesso un numero di iterazioni una matrice di varianze e
	#covarianza dalla distribuzione proposal in modo da avere un buon tasso di
	#accettazione ed una effective sample size di circa 1000 come richiesto
	#dall'esercizio
	#cominciamo con 1000 iterazioni (sicuramente un numero troppo ottimistico)
	#e con una matrice di varianza e covarianza simile a quella della g-prior:
	nsimul<-1000
	var.prop <-var(log(y+1/2))*solve(t(X)%*%X)
	var.prop

	 
}