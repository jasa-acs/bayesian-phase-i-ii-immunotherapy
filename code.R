## R version 3.0.0 was used. The following R packages were used: MCMCpack, dlm, msm, mvtnorm and cubature.
## The code takes roughly 20 hours to complete 100 simulated trials. Suggest to use computer cluster if a large number of simulations are needed.


## sensitivity analyses in the manuscript

## replace the original utility 
## Uti <- array(0,c(2,3,2)) # order: tox, eff, immuno
## Uti[,,1] <- matrix(c(0,0,50,10,80,35),nrow=2)
## Uti[,,2] <- matrix(c(5,0,70,20,100,45),nrow=2)
## with the following
## Uti <- array(0,c(2,3,2)) # order: tox, eff, immuno
## Uti[,,1] <- matrix(c(5,0,50,30,80,55),nrow=2)
## Uti[,,2] <- matrix(c(10,5,70,40,100,60),nrow=2)

## sensitivity on the sample size 
## replace N.max <- 60
## with N.max <- 42

## sensitivity on prior estimates of alphas
## replace theta <- c(1,5,.5,2)
## with theta <- c(2,4,.4,3) or 
## theta <- c(1.5,6,.6,2.5)


library(MCMCpack);library(dlm);library(msm);library(mvtnorm);library(cubature);#library(R2Cuba)

rm(list = ls())

set.seed(seed=1)



N.post <- 2000
N.burnin <- 500

doses <- c(.1,.3,.5,.7,.9)
n.dose <- length(doses)

#Utility
Uti <- array(0,c(2,3,2)) # order: tox, eff, immuno
Uti[,,1] <- matrix(c(0,0,50,10,80,35),nrow=2)
Uti[,,2] <- matrix(c(5,0,70,20,100,45),nrow=2)


cut.Y.I <- 4


dens <- function(x,dose,rho,alpha0,alpha1,alpha2,alpha3,beta0,beta1,beta2,beta3,gamma0,gamma1,gamma2,sigma2,m.Y.I,sd.Y.I,m.d,sd.d) {
	x1 <- x[1] 
	x2 <- x[2] 
	y <- x[3] 
	y.s <- (y-m.Y.I)/(2*sd.Y.I)
	if (sd.d==0) dose.s <- dose-m.d
	if (sd.d!=0) dose.s <- (dose-m.d)/(2*sd.d)
	mu1 <- beta0+beta1*dose.s+ifelse(y>beta3,1,0)*beta2*y.s
	mu2 <- gamma0+gamma1*y.s+gamma2*y.s^2
	mu.y <- alpha0+alpha1*(dose^alpha3)/(alpha2^alpha3+dose^alpha3)
	z <- (x1-mu1)^2+(x2-mu2)^2-2*rho*(x1-mu1)*(x2-mu2)
	return(exp(-(y-mu.y)^2/(2*sigma2)-z/(2*(1-rho^2)))/(2*sqrt(2*pi*sigma2)*pi*sqrt(1-rho^2)))
}


tox.prob <- function(v,dose) {
	rho <- v[1]; beta0 <- v[2]; beta1 <- v[3]; beta2 <- v[4]; beta3 <- v[5];
	gamma0 <- v[6]; gamma1 <- v[7]; gamma2 <- v[8]; sigma2 <- v[9]; 
	m.Y.I <- v[10]; sd.Y.I <- v[11]; m.d <- v[12]; sd.d <- v[13];
	alpha0 <- v[14]; alpha1 <- v[15]; alpha2 <- v[16]; alpha3 <- v[17]; xi <- v[18]
	Y_I_mean <- alpha0+alpha1*(dose^alpha3)/(alpha2^alpha3+dose^alpha3)

	y.sample <- rnorm(4000,Y_I_mean,sd=sqrt(sigma2))
	return(mean(pnorm(0,mean=(beta0+beta1*(dose-m.d)/(2*ifelse(sd.d==0,.1,sd.d))+ifelse(y.sample>beta3,1,0)*beta2*(y.sample-m.Y.I)/(2*sd.Y.I)),sd=1,lower=F)))
}
eff.prob <- function(v,dose) {
	rho <- v[1]; beta0 <- v[2]; beta1 <- v[3]; beta2 <- v[4]; beta3 <- v[5];
	gamma0 <- v[6]; gamma1 <- v[7]; gamma2 <- v[8]; sigma2 <- v[9]; 
	m.Y.I <- v[10]; sd.Y.I <- v[11]; m.d <- v[12]; sd.d <- v[13];
	alpha0 <- v[14]; alpha1 <- v[15]; alpha2 <- v[16]; alpha3 <- v[17]; xi <- v[18]
	Y_I_mean <- alpha0+alpha1*(dose^alpha3)/(alpha2^alpha3+dose^alpha3)

	y.sample <- rnorm(4000,Y_I_mean,sd=sqrt(sigma2))
	return(mean(pnorm(0,mean=(gamma0+gamma1*(y.sample-m.Y.I)/(2*sd.Y.I)+gamma2*((y.sample-m.Y.I)/(2*sd.Y.I))^2),sd=1,lower=F)))
}


norm2d <- function(v,rho)
	rmvnorm(1,mean=v,sigma=matrix(c(1,rho,rho,1),nrow=2))

n.y <- 10000
get.utility <- function(v,dose) {

	rho <- v[1]; beta0 <- v[2]; beta1 <- v[3]; beta2 <- v[4]; beta3 <- v[5];
	gamma0 <- v[6]; gamma1 <- v[7]; gamma2 <- v[8]; sigma2 <- v[9]; 
	m.Y.I <- v[10]; sd.Y.I <- v[11]; m.d <- v[12]; sd.d <- v[13];
	alpha0 <- v[14]; alpha1 <- v[15]; alpha2 <- v[16]; alpha3 <- v[17]; xi <- v[18]

	prob <- array(0,c(2,3,2))
	Y.I.mean <- alpha0+alpha1*(dose^alpha3)/(alpha2^alpha3+dose^alpha3)
	y.sample <- rnorm(n.y,Y.I.mean,sd=sqrt(sigma2))
	

	mu.T <- beta0+beta1*(dose-m.d)/(2*ifelse(sd.d==0,.1,sd.d))+ifelse(y.sample>beta3,1,0)*beta2*(y.sample-m.Y.I)/(2*sd.Y.I)
	mu.E <- gamma0+gamma1*(y.sample-m.Y.I)/(2*sd.Y.I)+gamma2*((y.sample-m.Y.I)/(2*sd.Y.I))^2
	mu.TE <- cbind(mu.T,mu.E)

	sample2d <- t(apply(mu.TE,1,norm2d,rho=rho))
	sample3d <- cbind(y.sample,sample2d)

	ll.y.I <- c(-Inf,cut.Y.I)
	uu.y.I <- c(cut.Y.I,Inf)

	prob <- array(0,c(2,3,2))
	for (i in 1:2) {
		prob[1,1,i] <- sum((sample3d[,1]<=uu.y.I[i]) & (sample3d[,1]>ll.y.I[i]) & (sample3d[,3]<=0) & (sample3d[,2]<=0))/n.y
		prob[1,2,i] <- sum((sample3d[,1]<=uu.y.I[i]) & (sample3d[,1]>ll.y.I[i]) & (sample3d[,3]>0) & (sample3d[,3]<xi) & (sample3d[,2]<=0))/n.y
		prob[1,3,i] <- sum((sample3d[,1]<=uu.y.I[i]) & (sample3d[,1]>ll.y.I[i]) & (sample3d[,3]>xi) & (sample3d[,2]<=0))/n.y
		prob[2,1,i] <- sum((sample3d[,1]<=uu.y.I[i]) & (sample3d[,1]>ll.y.I[i]) & (sample3d[,3]<=0) & (sample3d[,2]>0))/n.y
		prob[2,2,i] <- sum((sample3d[,1]<=uu.y.I[i]) & (sample3d[,1]>ll.y.I[i]) & (sample3d[,3]>0) & (sample3d[,3]<xi) & (sample3d[,2]>0))/n.y
		prob[2,3,i] <- sum((sample3d[,1]<=uu.y.I[i]) & (sample3d[,1]>ll.y.I[i]) & (sample3d[,3]>xi) & (sample3d[,2]>0))/n.y
	}
	uti <- sum(Uti[,,1]*prob[,,1]+Uti[,,2]*prob[,,2])
	return(uti)
}		


summary.mcmc <- function(matrix.uti) {

	mean.matrix.uti <- apply(matrix.uti,2,mean)

	
	uti <- rep(0,n.dose)
	tox.mcmc <- rep(0,n.dose)
	eff.mcmc <- rep(0,n.dose)
	for (i in 1:n.dose) {
		tox.prob.mcmc <- apply(matrix.uti,1,tox.prob,dose=doses[i])
		eff.prob.mcmc <- apply(matrix.uti,1,eff.prob,dose=doses[i])
		tox.mcmc[i] <- sum(tox.prob.mcmc<phi.T)/nrow(matrix.uti)
		eff.mcmc[i] <- sum(eff.prob.mcmc>phi.E)/nrow(matrix.uti)
		uti[i] <- get.utility(mean.matrix.uti,dose=doses[i])
	}
	list(tox=tox.mcmc,eff=eff.mcmc,uti=uti)
}


phi.T <- .3
phi.E <- .3



a1 <- .1 
b1 <- .1
a2 <- 1 
b2 <- 2
theta <- c(1,5,.5,2)
tau2 <- (4*theta)^2 
b.alpha <- theta/tau2
a.alpha <- theta^2/tau2
mu0 <- 0 

sigma0.2 <- (1.3)^2 
U_xi <- 6
c1 <- 0 
c2 <- 9



outcome <- function(dose.ind,coh.size,Y_I_mean.true,sigma2.true,beta0.true,beta1.true,beta2.true,beta3.true,gamma0.true,gamma1.true,gamma2.true,rho.true,xi.true) {
	Y_I <- rnorm(coh.size,Y_I_mean.true[dose.ind],sqrt(sigma2.true))

	Z_T_mean <- beta0.true+beta1.true*doses[dose.ind]+ifelse(Y_I>beta3.true,1,0)*beta2.true*Y_I
	Z_E_mean <- gamma0.true+gamma1.true*Y_I+gamma2.true*Y_I^2

	Z <- matrix(0,coh.size,2)
	for (i in 1:coh.size) 
		Z[i,] <- rmvnorm(1,mean=c(Z_T_mean[i],Z_E_mean[i]),sigma=matrix(c(1,rho.true,rho.true,1),nrow=2)) 

	Y_T <- ifelse(Z[,1]>0,1,0)
	Y_E <- rep(0,coh.size)
	for (i in 1:coh.size) {
		if (Z[i,2]>xi.true) Y_E[i] <- 3
		if (Z[i,2]<0) Y_E[i] <- 1
		if ((Z[i,2]>0) & (Z[i,2]<xi.true)) Y_E[i] <- 2
	}
	
	list(Y_I=Y_I,Y_T=Y_T,Y_E=Y_E)
}

log.likelihood.Y_I <- function(Y_I,d,alpha0,alpha1,alpha2,alpha3,sigma2) 
	-sum((Y_I-alpha0-alpha1*d^alpha3/(alpha2^alpha3+d^alpha3))^2)/(2*sigma2) 

log.likelihood.beta3 <- function(Y_I.uns,Y_I,Z_T,Z_E,d,beta0,beta1,beta2,beta3,gamma0,gamma1,gamma2,rho) {
	top.part1 <- sum((Z_T-beta0-beta1*d-ifelse(Y_I.uns>beta3,1,0)*beta2*Y_I)^2)
	top.part2 <- sum((Z_E-gamma0-gamma1*Y_I-gamma2*Y_I^2)^2)
	top.part3 <- -2*rho*sum((Z_T-beta0-beta1*d-ifelse(Y_I.uns>beta3,1,0)*beta2*Y_I)*(Z_E-gamma0-gamma1*Y_I-gamma2*Y_I^2))
	return(-(top.part1+top.part2+top.part3)/(2*(1-rho^2)))
}

log.likelihood.rho <- function(rho,Y_I.uns,Y_I,Z_T,Z_E,d,beta0,beta1,beta2,beta3,gamma0,gamma1,gamma2,n) {
	top.part1 <- sum((Z_T-beta0-beta1*d-ifelse(Y_I.uns>beta3,1,0)*beta2*Y_I)^2)
	top.part2 <- sum((Z_E-gamma0-gamma1*Y_I-gamma2*Y_I^2)^2)
	top.part3 <- -2*rho*sum((Z_T-beta0-beta1*d-ifelse(Y_I.uns>beta3,1,0)*beta2*Y_I)*(Z_E-gamma0-gamma1*Y_I-gamma2*Y_I^2))
	return(-(n/2)*log(1-rho^2)-(top.part1+top.part2+top.part3)/(2*(1-rho^2)))
}


mcmc.arms <- function(Y_I,Y_T,Y_E,d,n) {

	m.Y.I <- mean(Y_I); sd.Y.I <- sd(Y_I)
	m.d <- mean(d); sd.d <- sd(d)
	Y_I.s <- (Y_I-m.Y.I)/(2*sd.Y.I)
	if (sd.d==0) d.s <- d-m.d+rnorm(n,sd=.05)
	if (sd.d!=0) d.s <- (d-m.d)/(2*sd.d)

	sigma2_t <- rep(0,N.post)
	alpha0_t <- rep(0,N.post)
	alpha1_t <- rep(0,N.post)
	alpha2_t <- rep(0,N.post)
	alpha3_t <- rep(0,N.post)
	beta0_t <- rep(0,N.post)
	beta1_t <- rep(0,N.post)
	beta2_t <- rep(0,N.post)
	beta3_t <- rep(0,N.post)
	gamma0_t <- rep(0,N.post)
	gamma1_t <- rep(0,N.post)
	gamma2_t <- rep(0,N.post)
	Z_T_t <- matrix(0,N.post,n)
	Z_E_t <- matrix(0,N.post,n)
	xi_t <- rep(0,N.post)
	rho_t <- rep(0,N.post)
	
	# set initial values
	alpha0 <- 2
	alpha1 <- 5
	alpha2 <- .6
	alpha3 <- 2
	beta0 <- 0
	beta1 <- 0
	beta2 <- 0
	beta3 <- 6
	gamma0 <- 0
	gamma1 <- 0
	gamma2 <- 0
	xi <- 2
	rho <- .4

	Z_E_1 <- ifelse(Y_E==1,-.5,Y_E)
	Z_E_2 <- ifelse(Z_E_1==2,1,Z_E_1)
	Z_E <- ifelse(Z_E_2==3,9,Z_E_2)

	for (ite in 1:N.post) {

		sigma2 <- rinvgamma(1,n/2+a1,b1+.5*sum((Y_I-alpha0-alpha1*d^alpha3/(alpha2^alpha3+d^alpha3))^2))
		sigma2_t[ite] <- sigma2

		#alpha0

		logden <- function(x)
			log.likelihood.Y_I(Y_I,d,x,alpha1,alpha2,alpha3,sigma2)-(a.alpha[1]-1)*log(x)-b.alpha[1]*x
		alpha0 <- arms(alpha0,logden,function(x) ((x>0)*(x<5)),1)
		alpha0_t[ite] <- alpha0


		#alpha1
		logden <- function(x)
			log.likelihood.Y_I(Y_I,d,alpha0,x,alpha2,alpha3,sigma2)-(a.alpha[2]-1)*log(x)-b.alpha[2]*x
		alpha1 <- arms(alpha1,logden,function(x) ((x>3)*(x<9)),1)
		alpha1_t[ite] <- alpha1


		#alpha2
		logden <- function(x)
			log.likelihood.Y_I(Y_I,d,alpha0,alpha1,x,alpha3,sigma2)-(a.alpha[3]-1)*log(x)-b.alpha[3]*x
		alpha2 <- arms(alpha2,logden,function(x) ((x>0)*(x<2.5)),1)
		alpha2_t[ite] <- alpha2

		#alpha3
		logden <- function(x)
			log.likelihood.Y_I(Y_I,d,alpha0,alpha1,alpha2,x,sigma2)-(a.alpha[4]-1)*log(x)-b.alpha[4]*x
		alpha3 <- arms(alpha3,logden,function(x) ((x>0)*(x<5)),1)
		alpha3_t[ite] <- alpha3


		#xi, Z_T and Z_E
		
		min_xi <- ifelse(sum(Y_E==2)==0,0,max(max(subset(Z_E,Y_E==2)),0))
		max_xi <- ifelse(sum(Y_E==3)==0,U_xi,min(min(subset(Z_E,Y_E==3)),U_xi))
		xi <- runif(1,min_xi,max_xi)
		xi_t[ite] <- xi
		xi_cut <- c(-Inf,0,xi,Inf)
		zeta_cut <- c(-Inf,0,Inf)

		var.Z <- 1-rho^2
		mean.Z_T <- beta0+beta1*d.s+ifelse(Y_I>beta3,1,0)*beta2*Y_I.s+rho*(Z_E-gamma0-gamma1*Y_I.s-gamma2*Y_I.s^2)
		Z_T <- rtnorm(n,mean.Z_T,sqrt(var.Z),lower=zeta_cut[Y_T+1],upper=zeta_cut[Y_T+2])

		mean.Z_E <- gamma0+gamma1*Y_I.s+gamma2*Y_I.s^2+rho*(Z_T-beta0-beta1*d.s-ifelse(Y_I>beta3,1,0)*beta2*Y_I.s)

		Z_E <- rtnorm(n,mean.Z_E,sqrt(var.Z),lower=xi_cut[Y_E],upper=xi_cut[Y_E+1])

		Z_T_t[ite,] <- Z_T
		Z_E_t[ite,] <- Z_E

		#beta0
		mu.beta0.1 <- (sum(Z_T-beta1*d.s-ifelse(Y_I>beta3,1,0)*beta2*Y_I.s)-rho*sum(Z_E-gamma0-gamma1*Y_I.s-gamma2*Y_I.s^2))/n
		sigma2.beta0.1 <- (1-rho^2)/n
	
		sigma2.beta0 <- 1/(1/sigma2.beta0.1+1/sigma0.2)
		mu.beta0 <- sigma2.beta0*(mu.beta0.1/sigma2.beta0.1+mu0/sigma0.2)

		beta0 <- rnorm(1,mu.beta0,sqrt(sigma2.beta0))
		beta0_t[ite] <- beta0

		#beta1
		var.bottom <- sum(d.s^2)
		mu.beta1.1 <- (sum(d.s*(Z_T-beta0-ifelse(Y_I>beta3,1,0)*beta2*Y_I.s))-rho*sum(d.s*(Z_E-gamma0-gamma1*Y_I.s-gamma2*Y_I.s^2)))/var.bottom
		sigma2.beta1.1 <- (1-rho^2)/var.bottom
		
		sigma2.beta1 <- 1/(1/sigma2.beta1.1+1/sigma0.2)
		mu.beta1 <- sigma2.beta1*(mu.beta1.1/sigma2.beta1.1+mu0/sigma0.2)

		beta1 <- rnorm(1,mu.beta1,sqrt(sigma2.beta1))
		beta1_t[ite] <- beta1

		#beta2
		if (sum(Y_I>beta3)!=0) {
			Y_I.beta2 <- Y_I.s[Y_I>beta3]
			Z_T.beta2 <- Z_T[Y_I>beta3]
			Z_E.beta2 <- Z_E[Y_I>beta3]
			d.beta2 <- d.s[Y_I>beta3]
		
			var.bottom <- sum(Y_I.beta2^2)
			sigma2.beta2.1 <- (1-rho^2)/var.bottom
			mu.beta2.1 <- ( sum(Y_I.beta2*(Z_T.beta2-beta0-beta1*d.beta2))-rho*sum(Y_I.beta2*(Z_E.beta2-gamma0-gamma1*Y_I.beta2-gamma2*Y_I.beta2^2)) )/var.bottom

			sigma2.beta2 <- 1/(1/sigma2.beta2.1+1/sigma0.2)
			mu.beta2 <- sigma2.beta2*(mu.beta2.1/sigma2.beta2.1+mu0/sigma0.2)

			beta2 <- rnorm(1,mu.beta2,sqrt(sigma2.beta2))
		}
		beta2_t[ite] <- beta2

		#beta3

		logden <- function(x)
			log.likelihood.beta3(Y_I,Y_I.s,Z_T,Z_E,d.s,beta0,beta1,beta2,x,gamma0,gamma1,gamma2,rho)
		beta3 <- arms(beta3,logden,function(x) ((x>c1)*(x<c2)),1)
		beta3_t[ite] <- beta3

		#gamma0
		mu.gamma0.1 <- (sum(Z_E-gamma1*Y_I.s-gamma2*Y_I.s^2)-rho*sum(Z_T-beta0-beta1*d.s-ifelse(Y_I>beta3,1,0)*beta2*Y_I.s))/n
		sigma2.gamma0.1 <- (1-rho^2)/n

		sigma2.gamma0 <- 1/(1/(sigma2.gamma0.1)+1/(sigma0.2))
		mu.gamma0 <- sigma2.gamma0*(mu.gamma0.1/sigma2.gamma0.1+mu0/sigma0.2)
	
		gamma0 <- rnorm(1,mu.gamma0,sqrt(sigma2.gamma0))
		gamma0_t[ite] <- gamma0

		#gamma1
		mu.gamma1.1 <- (sum(Y_I.s*(Z_E-gamma0-gamma2*Y_I.s^2))-rho*sum(Y_I.s*(Z_T-beta0-beta1*d.s-ifelse(Y_I>beta3,1,0)*beta2*Y_I.s)))/sum(Y_I.s^2)
		sigma2.gamma1.1 <- (1-rho^2)/sum(Y_I.s^2)

		sigma2.gamma1 <- 1/(1/(sigma2.gamma1.1)+1/(sigma0.2))
		mu.gamma1 <- sigma2.gamma1*(mu.gamma1.1/sigma2.gamma1.1+mu0/sigma0.2)

		gamma1 <- rnorm(1,mu.gamma1,sqrt(sigma2.gamma1))
		gamma1_t[ite] <- gamma1

		#gamma2
		mu.gamma2.1 <- (sum(Y_I.s^2*(Z_E-gamma0-gamma1*Y_I.s))-rho*sum(Y_I.s^2*(Z_T-beta0-beta1*d.s-ifelse(Y_I>beta3,1,0)*beta2*Y_I.s)))/sum(Y_I.s^4)
		sigma2.gamma2.1 <- (1-rho^2)/sum(Y_I.s^4)

		sigma2.gamma2 <- 1/(1/(sigma2.gamma2.1)+1/(sigma0.2))
		mu.gamma2 <- sigma2.gamma2*(mu.gamma2.1/sigma2.gamma2.1+mu0/sigma0.2)
		
		gamma2 <- rnorm(1,mu.gamma2,sqrt(sigma2.gamma2))
		gamma2_t[ite] <- gamma2

		#rho
		logden <- function(x)
			log.likelihood.rho(x,Y_I,Y_I.s,Z_T,Z_E,d.s,beta0,beta1,beta2,beta3,gamma0,gamma1,gamma2,n)
		rho <- arms(rho,logden,function(x) ((x>0)*(x<1)),1)
		rho_t[ite] <- rho

		#print(ite)
	}
	
	ind <- seq((N.burnin+1),N.post)
	matrix.post <- cbind(alpha0_t,alpha1_t,alpha2_t,alpha3_t,beta0_t,beta1_t,beta2_t,beta3_t,gamma0_t,gamma1_t,gamma2_t,sigma2_t,xi_t)[ind,]
	matrix.uti <- cbind(rho_t,beta0_t,beta1_t,beta2_t,beta3_t,gamma0_t,gamma1_t,gamma2_t,sigma2_t,m.Y.I,sd.Y.I,m.d,sd.d,alpha0_t,alpha1_t,alpha2_t,alpha3_t,xi_t)[ind,]
	
	return(matrix.uti)
}

main <- function(coh.size,C_T,C_E) {

	
	
		alpha0.true <- 2
		alpha1.true <- 5
		alpha2.true <- .5
		alpha3.true <- 3
		sigma2.true <- 1
		beta0.true <- -2
		beta1.true <- 1.5
		beta2.true <- .1
		beta3.true <- 4
		gamma0.true <- -.9
		gamma1.true <- .5
		gamma2.true <- -.065
		rho.true <- .5
		xi.true <- .9
	
	Y_I_mean.true <- rep(0,n.dose)
	for (i in 1:n.dose) 
		Y_I_mean.true[i] <- alpha0.true+alpha1.true*(doses[i])^alpha3.true/(alpha2.true^alpha3.true+(doses[i])^alpha3.true)



	Y_I <- NULL; Y_T <- NULL; Y_E <- NULL; d <- NULL; n <- 0
	dose.ind <- 1

		while (n <= (N.max-coh.size)) {
			out <- outcome(dose.ind,coh.size,Y_I_mean.true,sigma2.true,beta0.true,beta1.true,beta2.true,beta3.true,gamma0.true,gamma1.true,gamma2.true,rho.true,xi.true)
			Y_I <- c(Y_I,out$Y_I)
			Y_E <- c(Y_E,out$Y_E)
			Y_T <- c(Y_T,out$Y_T)
			d <- c(d,rep(doses[dose.ind],coh.size))
	
			n <- n + coh.size
			h.dose <- 5*(max(unique(d))+.1)

			mcmc.tmp <- mcmc.arms(Y_I,Y_T,Y_E,d,n)
			summ.tmp <- summary.mcmc(mcmc.tmp)
			tox.summ <- summ.tmp$tox
			if ((tox.summ[h.dose]>.5) & (h.dose!=n.dose)) dose.ind <- h.dose+1 else {

	
				eff.summ <- summ.tmp$eff
				capA.ind <- which((tox.summ>C_T) & (eff.summ>C_E))
				if (length(capA.ind)==0) {
					dose.ind <- 0
					break
				} else {
					uti <- summ.tmp$uti[capA.ind]
	
					if (n==N.max) dose.ind <- capA.ind[which.max(uti)]
					if (n!=N.max) {
						uti.normalize <- uti/sum(uti)
						cumu.uti <- uti
						for (i in 1:length(uti)) 
						cumu.uti[i] <- sum(uti.normalize[1:i])
						r <- runif(1,0,1)
						dose.ind <- capA.ind[min(which(r<cumu.uti))]
							if (dose.ind>h.dose) dose.ind <- h.dose+1
					}
				}
			}
		}
	
	list(dose.ind=dose.ind,d=d)

}
N.max <- 60
n.sim <- 50
dose.recomm <- rep(0,n.sim)
alloc <- matrix(0,n.sim,N.max)

for (i in 1:n.sim) {
	sim.tmp <- main(coh.size=3,C_T=0.05,C_E=0.05)
	dose.recomm[i] <- sim.tmp$dose.ind
	if (length(sim.tmp$d)==N.max)
		alloc[i,] <- sim.tmp$d else
		alloc[i,] <- c(sim.tmp$d,rep(0,(N.max-length(sim.tmp$d))))
	print(sim.tmp$dose.ind)
	write.table(dose.recomm,"recomm.txt",quote=F,row.names=F)
	write.table(alloc,"alloc.txt",quote=F,row.names=F)
}
		
