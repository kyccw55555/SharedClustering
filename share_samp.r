############################################################################################################################################
##	"Shared Clustering" parameter sampling -- speed improved version
##	All denotations are the same as in data_gen.r once involved
############################################################################################################################################
#set.seed(4002)
############################################################################################################################################

rgenc<-function
(K,N)
{ 
#Generate valid C randomly
    mark=0
    while(mark!=1){
        C<-sample(1:K,replace=TRUE,N)
        if (length(which(table(C)<=1))==0)
			mark=1			
    }
    C
}

validc<-function
(K,C)
{	
    mark=0
	z=0
	if (length(which(table(C)<=1))==0)
		mark=1		
	else 
		z<-min(which(table(C)<=1))
    return(list("ind"=mark,"miss"=z))
}

extrb<-function
(Y,C,x,y)
{
	N<-ncol(Y)
	na<-(1+N)*N/2-N
	b<-rep(NA,na)
	diag(Y)<-rep(NA,N)
	b<-Y[C==x,C==y]
	b<-as.vector(b[!is.na(b)])
	if (x==y)
		b<-c(rep(1,sum(b)/2),rep(0,(length(b)-sum(b))/2))
	b
}

vectomat<-function
(c,pr)
{
#Probability Vector To Matrix
	prmat<-matrix(NA,ncol=c,nrow=c)
	prmat[upper.tri(prmat, diag=TRUE)]<-pr
	prmat[lower.tri(prmat)]<-t(prmat[upper.tri(prmat)])
	prmat
}

lndMvn<-function	#Exactly as in {bayesm}, this function can only deal with one single observation
(x,mu,rooti)
{
#
# changed 12/05 by Rossi to include normalizing constant
#
# function to evaluate log of MV NOrmal density with  mean mu, var Sigma
# Sigma=t(root)%*%root   (root is upper tri cholesky root)
# Sigma^-1=rooti%*%t(rooti)   
# rooti is in the inverse of upper triangular chol root of sigma
#          note: this is the UL decomp of sigmai not LU!
#                Sigma=root'root   root=inv(rooti)
#
z=as.vector(t(rooti)%*%(x-mu))
return(  -(length(x)/2)*log(2*pi) -.5*(z%*%z) + sum(log(diag(rooti))))
}

lndIWishart<-function	#Exactly as in {bayesm}
(nu,V,IW)
{
# 
# P. Rossi 12/04
#
# purpose: evaluate log-density of inverted Wishart
#    includes normalizing constant
#
# arguments:
#   nu is d. f. parm
#   V is location matrix
#   IW is the value at which the density should be evaluated
#
# output:
#   value of log density
#
# note: in this parameterization, E[IW]=V/(nu-k-1)
#
k=ncol(V)
Uiw=chol(IW)
lndetVd2=sum(log(diag(chol(V))))
lndetIWd2=sum(log(diag(Uiw)))
#
# first evaluate constant
#
const=((nu*k)/2)*log(2)+((k*(k-1))/4)*log(pi)
arg=(nu+1-c(1:k))/2
const=const+sum(lgamma(arg))
return(-const+nu*lndetVd2-(nu+k+1)*lndetIWd2-.5*sum(diag(V%*%chol2inv(Uiw))))
}

rwishart<-function	#Exactly as in {bayesm}
(nu,V)
{
#
# function to draw from Wishart (nu,V) and IW
# 
# W ~ W(nu,V)
# E[W]=nuV
#
# WI=W^-1
# E[WI]=V^-1/(nu-m-1)
# 
#
m=nrow(V)
df=(nu+nu-m+1)-(nu-m+1):nu
if(m >1) {
T=diag(sqrt(rchisq(c(rep(1,m)),df)))	#rchisq(c(rep(1,m)),df) is equivalent to rchisq(m,df)
T[lower.tri(T)]=rnorm((m*(m+1)/2-m))}
else
{T=sqrt(rchisq(1,df))}
U=chol(V)
C=t(T)%*%U
CI=backsolve(C,diag(m))	#Inverse of C, same to solve(C), this way faster?
#
#   C is the upper triangular root of Wishart
#      therefore, W=C'C  this is the LU decomposition 
#      Inv(W) = CICI'  Note:  this is the UL decomp not LU!
#
return(list(W=crossprod(C),IW=crossprod(t(CI)),C=C,CI=CI))	#W=t(C)%*%C
#  W is Wishart draw,  IW is W^-1
}

rdirichlet<-function(n,alpha)
## generate n random deviates from the Dirichlet function with shape
## parameters alpha
{
    l<-length(alpha);
    x<-matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE);
    sm<-x%*%rep(1,l);
    x/as.vector(sm);
}

ddirichlet<-function(x,alpha)
## probability density for the Dirichlet function, where x=vector of
## probabilities
## and (alpha-1)=vector of observed samples of each type
## ddirichlet(c(p,1-p),c(x1,x2)) == dbeta(p,x1,x2)
{

  dirichlet1 <- function(x, alpha)
    {
      logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
      s<-sum((alpha-1)*log(x))
      exp(sum(s)-logD)

    }

  # make sure x is a matrix
  if(!is.matrix(x))
    if(is.data.frame(x))
      x <- as.matrix(x)
    else
      x <- t(x)

  if(!is.matrix(alpha))
    alpha <- matrix( alpha, ncol=length(alpha), nrow=nrow(x), byrow=TRUE)

  if( any(dim(x) != dim(alpha)) )
    stop("Mismatch between dimensions of 'x' and 'alpha'.")

  pd <- vector(length=nrow(x))
  for(i in 1:nrow(x))
    pd[i] <- dirichlet1(x[i,],alpha[i,])

  # Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
  pd[ apply( x, 1, function(z) any( z <0 | z > 1)) ] <- 0
  pd[ apply( x, 1, function(z) all.equal(sum( z ),1) !=TRUE) ] <- 0
  pd
}

rmultireg<-function	#Exactly as in {bayesm}
(Y,X,Bbar,A,nu,V)	#X here is a nobincomp*1 vector with all elements one
{
#
# revision history:
#    changed 1/11/05 by P. Rossi to fix sum of squares error
#
# purpose:
#    draw from posterior for Multivariate Regression Model with
#    natural conjugate prior
# arguments:
#    Y is n x m matrix
#    X is n x k
#    Bbar is the prior mean of regression coefficients  (k x m)
#    A is prior precision matrix
#    nu, V are parameters for prior on Sigma
# output:
#    list of B, Sigma draws of matrix of coefficients and Sigma matrix
# model:
#    Y=XB+U  cov(u_i) = Sigma
#    B is k x m matrix of coefficients
# priors:  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
#                   betabar=vec(Bbar)
#                   beta = vec(B) 
#          Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)
n=nrow(Y)
m=ncol(Y)
k=ncol(X)	#k=1
#
# first draw Sigma
#
RA=chol(A)	#Just sqrt(A[1,1])
W=rbind(X,RA)
Z=rbind(Y,RA%*%Bbar)
#   note:  Y,X,A,Bbar must be matrices!
IR=backsolve(chol(crossprod(W)),diag(k))
#                      W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
Btilde=crossprod(t(IR))%*%crossprod(W,Z)   
#                      IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
S=crossprod(Z-W%*%Btilde)
#                      E'E
rwout=rwishart(nu+n,chol2inv(chol(V+S)))
#
# now draw B given Sigma
#   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
#       Cov=(X'X + A)^-1  = IR t(IR)  
#       Sigma=CICI'    
#       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
#	so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
#			Z_mk is m x k matrix of N(0,1)
#	since vec(ABC) = (C' (x) A)vec(B), we have 
#		B = Btilde + IR Z_mk CI'
#
B = Btilde + IR%*%matrix(rnorm(m*k),ncol=m)%*%t(rwout$CI)
return(list(B=B,Sigma=rwout$IW))
}

rcompsC<-function
(index,xi,px,u,E,C,Y,N,psi)
{
#xi is q*1 vector - one observation
#px is prior probability 
	K<-max(unique(C)) 
	q<-ncol(u)
	ll<-rep(0,K)
	poc<-rep(0,K)
	hc<-rep(0,K)
	for (i in 1:K){
		rooti<-backsolve(chol(E[,,i]),diag(q))
		ll[i]<-lndMvn(xi,u[i,],rooti)
		Yindex<-Y[index,-index]
		Ci<-C[-index]
		hc[i]<-sum( Yindex[Ci>=i]*log(psi[Ci[Ci>=i]*(Ci[Ci>=i]-1)/2+i])
			+(1-Yindex[Ci>=i])*log(1-psi[Ci[Ci>=i]*(Ci[Ci>=i]-1)/2+i])
			)+
			sum(
			Yindex[Ci<i]*log(psi[i*(i-1)/2+Ci[Ci<i]])
			+(1-Yindex[Ci<i])*log(1-psi[i*(i-1)/2+Ci[Ci<i]])
			)
	}
	max<-max(ll)
	poc<-exp(ll-rep(max,K)+hc)*px
	poc<-poc/sum(poc)
	rand<-runif(1)
	label<-min(which(cumsum(poc)>=rand))
	return(list("label"=label,"poc"=poc))
}

Gibbsone<-function
(X,u0,A,v0,T,a,px,C,E,u,Y,psi,beta)	
{
	#Sample phi (E, u) for X and psi for Y
	N<-nrow(X)
	K<-max(unique(C))	
	u0<-t(as.matrix(u0))	
	for(i in 1:K){
		nobincomp<-sum(C==i)         # get number of observations "in" component i
		if (nobincomp>0) {             # if more than one obs in this component, draw from posterior
			Xi<-X[C==i,]
			dim(Xi)<-c(nobincomp,ncol(X))
			temp<-rmultireg(Xi,matrix(rep(1,nobincomp),ncol=1),u0,A,v0,T)
			u[i,]<-as.vector(temp$B)
			E[,,i]<-temp$Sigma	
		} 
		else { # else draw from the prior
			rw<-rwishart(v0,chol2inv(chol(T)))	
			#T is the prior of the Scale Matrix Sigma, thus chol2inv(chol(T)) is the "prior" of Precision 
			u[i,]<-as.vector(t(u0) + (rw$CI %*% rnorm(ncol(X)))/sqrt(A[1,1]))
			E[,,i]<-rw$IW
		}
		for (j in 1:i){
			b<-extrb(Y,C,i,j)
			psi[i*(i-1)/2+j]<-rbeta(1,sum(b)+beta,length(b)-sum(b)+beta)
		}
	}
	#Sample C -- the most time consuming part
	dc<-rep(0,N)
	for (index in 1:N){
		sampc<-rcompsC(index,X[index,],px,u,E,C,Y,N,psi)
		C[index]<-sampc$label
		judge<-validc(K,C)
		if (judge$ind!=1)
			C[index]<-judge$miss
		dc[index]<-sampc$poc[C[index]]
		#cat(index,"dc_i",dc[index],"\n")
	}		
	#Sample px
	a<-a+table(C)
	px<-rdirichlet(1,a)
	dpx<-ddirichlet(px,a)	#Density but not Log-density
	
	out<-NULL
	out$px<-px; out$dpx<-dpx; out$C<-C; out$dc<-dc; out$E<-E; out$u<-u; out$psi<-psi
	out
}

jopolik<-function	#Calculate joint posterior Log-likelihood at one stage
(dat,u,E,v0,T,u0,A,psi,beta,C) 
{	
	K<-length(dat$Clus$cn)
	jplik<-0
	#For vectorial data one by one i.e. 1:N
	for (i in 1:nrow(dat$X))
			jplik<-(jplik
				+lndMvn(dat$X[i,],u[C[i],],backsolve(chol(E[,,C[i]]),diag(ncol(dat$X))))
				+lndIWishart(v0,T,E[,,C[i]])
				+lndMvn(u[C[i],],u0,backsolve(chol(E[,,C[i]]*(A[1,1]^-1)),diag(ncol(dat$X)))))
	#For network data
	for (i in 1:K){
		for (j in 1:i){
			b<-extrb(dat$Y,C,i,j)
			jplik<-(jplik
				+sum( log(psi[i*(i-1)/2+j])*(b==1) )+sum( log(1-psi[i*(i-1)/2+j])*(b==0) )
				+log(gamma(2*beta)/(gamma(beta)^2))+log(psi[i*(i-1)/2+j])*(beta-1)+log(1-psi[i*(i-1)/2+j])*(beta-1))
		}#for j
	}#for i
	jplik
}

run<-function
(dat,cnt,burnin=1000)
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
	#px<-c(0.1,0.1,0.8)
	C<-rgenc(K,nrow(dat$X))
	E<-array(0,c(q,q,K))
	u<-matrix(0,ncol=q,nrow=K)
	beta<-log10(ncol(dat$Y)/K)
	psi_vec<-runif(K*(K+1)/2)
	plotlik<-rep(0,cnt)
	maxjp<--Inf
	E_pmode<-E
	u_pmode<-u
	psi_vec_pmode<-psi_vec
	px_pmode<-px
	C_pmode<-C
	#cat("Initial labels:",C,"\n")
	for(rep in 1:cnt){
		one<-Gibbsone(dat$X,u0,A,v0,T,a,px,C,E,u,dat$Y,psi_vec,beta)
		E<-one$E
		u<-one$u
		psi_vec<-one$psi
		px<-one$px
		C<-one$C
		#cat("number",rep,"label:	",C,"\n")
		#cat("#",rep,class(one$dc),"\n")
		plotlik[rep]<-jopolik(dat,u,E,v0,T,u0,A,psi_vec,beta,C)+log(one$dpx)+sum(log(one$dc))
		if (plotlik[rep]>maxjp && rep>burnin){
			maxjp<-plotlik[rep]
			E_pmode<-E
			u_pmode<-u
			psi_vec_pmode<-psi_vec
			px_pmode<-px
			C_pmode<-C
		}
	}
	#plot(plotlik,type="l",main="Joint Posterior Log-Likelihood",xlab="cnt",ylab="value")
	PSI<-vectomat(K,psi_vec_pmode)
	res<-NULL
	res$px<-px_pmode; res$C<-C_pmode; res$E<-E_pmode; res$u<-u_pmode; res$psi<-PSI
	#cat("sample finished.\n")
	res
}
