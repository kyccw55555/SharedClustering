#########################################################################################################################
#"Shared Clustering" data generating
#DENOTATION:
#Total number of nodes:N
#Total number of clusters K
#Cluster label vector C(length N)
#pp is a vector containing partition points for C
#Vectorial data X, with each component X[i,] and the dimension q
#Network data Y, adjacency matrix(N*N) with diagonal elements all 0
#Parameters:
#Gaussian: PHI={u,E}, where u is a K*q matrix and E is an array(with length K) of covariance matrix(q*q)
#Bernoulli: PSI=P, a K*K symmetric matrix containing Bernoulli probability 
#########################################################################################################################
#set.seed(4002)
library(MASS)	#mvrnorm()
#------------------------------------------------------------------------------------------------------------------------

#Generate parameter set PSI
genpsi<-function
(K,iimin=0.8,iimax=0.8,ijmin=0.2,ijmax=0.2)
{ 
    p<-matrix(ncol=K,nrow=K)
    for (i in 1:K){
        for (j in 1:K){
            if (i==j) p[i,j]=runif(1,min=iimin,max=iimax)
            else p[i,j]=runif(1,min=ijmin,max=ijmax)
            p[j,i]=p[i,j]
        }
    }    
    p
}

#Generate parameter set PHI
genphi<-function
(K,q=2,umin=1,umax=10,Emin=1,Emax=1,rhomin=0.5,rhomax=0.5)
{
	PHI<-NULL
	E<-array(0,c(q,q,K))
	u<-matrix(0,ncol=q,nrow=K)
	rownames(u)<-1:K
	colnames(u)<-1:q
	u[,1]<-sample(umin:umax,K,replace=F)
	#u[,1]<-(1:K)
	#u[,1]<-c(-10,0,10)
	for (i in 2:q)
		u[,i]<-u[,1]
	for (k in 1:K){
		for (i in 1:q)
			for (j in 1:i){
				if (i==j)
					E[i,j,k]<-runif(1,min=Emin,max=Emax)
				else
					E[i,j,k]<-runif(1,min=rhomin,max=rhomax)
			}
		ind<-upper.tri(E[,,k])
		E[,,k][ind]<-t(E[,,k])[ind]
	}
	PHI$u<-u; PHI$E<-E; PHI$q<-q;
	PHI
}

#Generate clustering structure, the vector C
genc<-function
(K, N, pp=c((N/K)*(1:(K-1)),N) )
{
#pp:"partiton points"
	#cat("pp:",pp,"\n")
	Clus<-NULL
	C<-rep(0,N)
	C[1:pp[1]]<-1
	for (i in 2:K)
		C[(pp[i-1]+1):pp[i]]<-i
    cn<-rep(0,K) #cluster numbers
	cn[1]<-pp[1]
	for (i in 2:K)
		cn[i]<-pp[i]-pp[i-1]
	Clus$C<-C; Clus$cn<-cn; Clus$pp<-pp
	Clus
}

#Generate adjacency matrix Y
geny<-function
(N,K,P,ccn) #ccn is the object returned by genc()
{
	C<-ccn$C
	cn<-ccn$cn
	Y<-matrix(rbinom(cn[1]*cn[1],1,P[1,1]),ncol=cn[1],nrow=cn[1])
	for (index2 in 2:K)
		Y<-cbind(Y, matrix(rbinom(cn[index2]*cn[1],1,P[1,index2]),ncol=cn[index2],nrow=cn[1]) )
	for (index1 in 2:K){
		tempY<-matrix(rbinom(cn[index1]*cn[1],1,P[index1,1]),ncol=cn[index1],nrow=cn[1])
		for (index2 in 2:K){
			tempY<-cbind(tempY, matrix(rbinom(cn[index2]*cn[index1],1,P[index1,index2]),ncol=cn[index2],nrow=cn[index1]) )
		}
		Y<-rbind(Y,tempY)
	}
	ind<-lower.tri(Y) #Based on the generated data for lower triangle
	Y[ind]<-t(Y)[ind]
	#Forbid self loop
	for (i in 1:N)
		Y[i,i]<-0
	#For easy detection
	rownames(Y)<-(1:N)
	colnames(Y)<-C
	Y
}

#Generate vectorial data X
genx<-function
(N,K,phi,ccn)
{
	u<-phi$u
	E<-phi$E
	q<-phi$q
	C<-ccn$C
	pp<-ccn$pp
	cn<-ccn$cn
	X<-mvrnorm(pp[1],u[1,],E[,,1])
	for (k in 2:K)
		X<-rbind(X,mvrnorm(pp[k]-pp[k-1],u[k,],E[,,k]))
	rownames(X)<-C
	colnames(X)<-1:q
	X
}

#Generate the composite data
gendata<-function
(N,K,q=2)
{
	dat<-NULL
	PSI<-genpsi(K)
	PHI<-genphi(K,q)
	ccn<-genc(K,N)
	Y<-geny(N,K,PSI,ccn)
	X<-genx(N,K,PHI,ccn)
	dat$Clus<-ccn; dat$X<-X; dat$PHI<-PHI; dat$Y<-Y; dat$PSI<-PSI
	dat
}