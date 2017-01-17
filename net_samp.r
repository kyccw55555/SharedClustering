############################################################################################################################################
##	"Shared Clustering" parameter sampling -- speed improved version for network data only
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

rcompsCy<-function
(index,px,C,Y,N,psi)
{
#xi is q*1 vector - one observation
#px is prior probability 
	K<-max(unique(C)) 
	#q<-ncol(u)
	poc<-rep(0,K)
	hc<-rep(0,K)
	for (i in 1:K){
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
	poc<-exp(hc)*px
	poc<-poc/sum(poc)
	rand<-runif(1)
	label<-min(which(cumsum(poc)>=rand))
	return(list("label"=label,"poc"=poc))
}

Gibbsoney<-function
(a,px,C,Y,psi,beta)	
{
	#Sample phi (E, u) for X and psi for Y
	N<-ncol(Y)
	K<-max(unique(C))		
	for(i in 1:K){
		for (j in 1:i){
			b<-extrb(Y,C,i,j)
			psi[i*(i-1)/2+j]<-rbeta(1,sum(b)+beta,length(b)-sum(b)+beta)
		}
	}
	#Sample C -- the most time consuming part
	dc<-rep(0,N)
	for (index in 1:N){
		sampc<-rcompsCy(index,px,C,Y,N,psi)
		C[index]<-sampc$label
		judge<-validc(K,C)
		if (judge$ind!=1)
			C[index]<-judge$miss
		dc[index]<-sampc$poc[C[index]]
	}		
	#Sample px
	a<-a+table(C)
	px<-rdirichlet(1,a)
	dpx<-ddirichlet(px,a)
	
	out<-NULL
	out$px<-px; out$dpx<-dpx; out$C<-C; out$dc<-dc; out$psi<-psi
	out
}

jopoliky<-function	#Calculate joint posterior likelihood at one stage
(dat,psi,beta,C) 
{	
	K<-length(dat$Clus$cn)
	jplik<-0
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

runy<-function
(dat,cnt,burnin=1000)
{
	#Extract data information and initialize hyper-parameters
	K<-length(dat$Clus$cn)
	q<-ncol(dat$X)
	a<-rep(5,K)
	px<-c(rep(1,K))/K
	C<-rgenc(K,nrow(dat$X))
	beta<-log10(ncol(dat$Y)/K)
	psi_vec<-runif(K*(K+1)/2)
	plotlik<-rep(0,cnt)
	maxjp<--Inf
	psi_vec_pmode<-psi_vec
	px_pmode<-px
	C_pmode<-C
	#cat("Initial labels:",C,"\n")
	for(rep in 1:cnt){
		one<-Gibbsoney(a,px,C,dat$Y,psi_vec,beta)
		psi_vec<-one$psi
		px<-one$px
		C<-one$C
		#cat("number",rep,"label:	",C,"\n")
		#cat("#",rep,"Success!\n")		
		plotlik[rep]<-jopoliky(dat,psi_vec,beta,C)+log(one$dpx)+sum(log(one$dc))
		if (plotlik[rep]>maxjp && rep>burnin){
			maxjp<-plotlik[rep]
			psi_vec_pmode<-psi_vec
			px_pmode<-px
			C_pmode<-C
		}
	}
	#plot(plotlik,type="l",main="Joint Posterior Log-Likelihood -- network data",xlab="cnt",ylab="value")
	PSI<-vectomat(K,psi_vec_pmode)
	res<-NULL
	res$px<-px_pmode; res$C<-C_pmode; res$psi<-PSI
	#cat("sample finished.\n")
	res
}
