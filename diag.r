############################################################################################################################################
##	"Shared Clustering" diagnostic 
##	All denotations are the same as in data_gen.r once involved
############################################################################################################################################
library(sna)	#potscalered.mcmc()
#library(mclust)	#adjustedRandIndex() code is copied
#-------------------------------------------------------------------------------------------------------------------------------------------

adjustedRandIndex <- function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

pseudo_ari<-function
(C,C1,C2)
{
	N<-length(C)
	K<-max(C)
	Ctemp<-rep(0,N)
	tab<-table(C1,C2)
	for (k in 1:(K-1)){
		wheremax<-which(tab==max(tab),arr.ind = TRUE)[1,]
		wheremax<-as.numeric(c(rownames(tab)[wheremax[1]],colnames(tab)[wheremax[2]]))
		Ctemp[C2==wheremax[2]]<-wheremax[1]
		tab<-tab[-wheremax[1],-wheremax[2]]
	}
	for (i in 1:K)
		if (length(which(Ctemp==i))==0)
			break
	Ctemp[Ctemp==0]<-i
	cont<-(table(C,C1)+table(C,Ctemp))/2
	b<-colSums(cont)
	a<-rowSums(cont)
	nomi<- sum(cont*(cont-1)/2) - sum(a*(a-1)/2) * sum(b*(b-1)/2)/( N*(N-1)/2 )
	denom<- 0.5*( sum(a*(a-1)/2) + sum(b*(b-1)/2) ) - sum(a*(a-1)/2) * sum(b*(b-1)/2)/( N*(N-1)/2 )
	ARI<-nomi/denom
	ARI
}

tra<-function
(dat,chain=3,cnt=100)
{
	#Extract data information and initialize hyper-parameters
	K<-length(dat$Clus$cn)
	q<-ncol(dat$X)
	u0<-rep(0,q)
	A<-matrix(0.01,1,1)
	v0<-q+3
	T<-v0*diag(q)
	a<-rep(5,K)
	px<-c(rep(1,K))/K
	C<-rgenc(K,nrow(dat$X))
	E<-array(0,c(q,q,K))
	u<-matrix(0,ncol=q,nrow=K)
	beta<-log10(ncol(dat$Y)/K)
	psi_vec<-runif(K*(K+1)/2)
	plotlik<-matrix(0,nrow=cnt,ncol=chain)
	pdis<-matrix(0,nrow=cnt,ncol=chain)
	udis<-matrix(0,nrow=cnt,ncol=chain)
	Edis<-matrix(0,nrow=cnt,ncol=chain)
	#pxdis<-matrix(0,nrow=cnt,ncol=chain)
	cmat<-matrix(0,nrow=chain,ncol=nrow(dat$X))
	for(rep in 1:cnt){
		one<-Gibbsone(dat$X,u0,A,v0,T,a,px,C,E,u,dat$Y,psi_vec,beta)
		E<-one$E
		Edis[rep,1]<-sqrt(sum((as.vector(E)-as.vector(dat$PHI$E))^2))
		u<-one$u
		udis[rep,1]<-sqrt(sum((as.vector(u)-as.vector(dat$PHI$u))^2))
		psi_vec<-one$psi
		pdis[rep,1]<-sqrt(sum((psi_vec-as.vector(dat$PSI[upper.tri(dat$PSI,diag=TRUE)]))^2))
		px<-one$px
		C<-one$C
		plotlik[rep,1]<-jopolik(dat,u,E,v0,T,u0,A,psi_vec,beta,C)
		cmat[1,]<-C
	}
	cat(1,"sample finished.\n")
	par(mfrow=c(2,2))
	plot(plotlik[,1],type="l",main="Joint Posterior Log-Likelihood",xlab="cnt",ylab="value")	
	for (ch in 2:chain){
		px<-c(rep(1,K))/K
		C<-rgenc(K,nrow(dat$X))
		E<-array(0,c(q,q,K))
		u<-matrix(0,ncol=q,nrow=K)
		psi_vec<-runif(K*(K+1)/2)
		for(rep in 1:cnt){
			one<-Gibbsone(dat$X,u0,A,v0,T,a,px,C,E,u,dat$Y,psi_vec,beta)
			E<-one$E
			Edis[rep,ch]<-sqrt(sum((as.vector(E)-as.vector(dat$PHI$E))^2))
			u<-one$u
			udis[rep,ch]<-sqrt(sum(( c(sort(u[,1]),sort(u[,2])) - c(sort(dat$PHI$u[,1]),sort(dat$PHI$u[,2])) )^2))
			psi_vec<-one$psi
			pdis[rep,ch]<-sqrt(sum((psi_vec-as.vector(dat$PSI[upper.tri(dat$PSI,diag=TRUE)]))^2))
			px<-one$px
			C<-one$C
			plotlik[rep,ch]<-jopolik(dat,u,E,v0,T,u0,A,psi_vec,beta,C)
			cmat[ch,]<-C
		}
		cat(ch,"sample finished.\n")
		lines(plotlik[,ch],col=rainbow(chain)[ch])
	}	
	plot(Edis[,1],type="l",main="parameter: E",xlab="cnt",ylab="distance")
	for (ch in 2:chain)
		lines(Edis[,ch],col=rainbow(chain)[ch])
	plot(udis[,1],type="l",main="parameter: u",xlab="cnt",ylab="distance")
	for (ch in 2:chain)
		lines(udis[,ch],col=rainbow(chain)[ch])
	plot(pdis[,1],type="l",main="parameter: psi",xlab="cnt",ylab="distance")
	for (ch in 2:chain)
		lines(pdis[,ch],col=rainbow(chain)[ch])
	C_true<-dat$Clus$C
	cat("Results of clustering:\n")
	#prmatrix(cmat)
	par(mfrow=c(1,1))
	ari<-rep(0,chain)
	for (i in 1:chain)
		ari[i]<-adjustedRandIndex(C_true,cmat[i,])
	return(list=c("rhat"=potscalered.mcmc(plotlik),"ari"=ari))
}