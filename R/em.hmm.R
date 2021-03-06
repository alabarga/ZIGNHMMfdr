
em.hmm = function(x, alttype="mixnormal", L=2, maxiter=1000, nulltype=2, symmetric= FALSE, seed = 100, burn = 200, ptol=1e-2, core = 2, v = F)
{
	NUM = length(x)
	ptol = 1e-2
	difference = 1
	niter = 0
	tmpdir = paste(paste("fdrtmp",format(Sys.time(), "%H_%M_%S"),sep="_"),abs(rnorm(1))*10e15,sep="_")
	dir.create(tmpdir)
	EMvar = tryCatch({ em.hmm.runseed(x, alttype, L, seed, burn, nulltype, maxiter, symmetric, ptol, core, v = v, tmpdir=tmpdir) }, error = function(e) {unlink(tmpdir,r=T)})
	best_EMvar = list()
	while(length(best_EMvar) <= 1)
	{
		niter = niter + 1
		best_EMvar = em.hmm.loadseed(EMvar[[niter]], tmpdir)
		if(v) print(paste("best seed : ",best_EMvar$logL))
		best_EMvar = try(em.hmm.EM(x, best_EMvar, alttype, L, maxiter, nulltype, symmetric, ptol, E=T, v = v))
		if(length(EMvar) <= 1)
		{
			if(v) print("error: Trying next best seed")
			best_EMvar = em.hmm.loadseed(EMvar[[niter]], tmpdir)
		}
	}
	lfdr = best_EMvar$gamma[, 1]
	unlink(tmpdir,r=T)
	if(length(best_EMvar) > 1)
	{
		logL  =  best_EMvar$logL
		if(v) print(logL)
		if (nulltype > 0) {
			BIC = logL-(3*L+3)*log(NUM)/2 
		} else {
			BIC = logL-(3*L+2)*log(NUM)/2 
		}
		em.var = list(ptheta = best_EMvar$ptheta, pii = best_EMvar$pii, A = best_EMvar$A, pc = best_EMvar$pc, f0 = best_EMvar$f0, f1 = best_EMvar$f1, LIS=lfdr, logL=logL, BIC=BIC, gamma = best_EMvar$gamma, dgamma = best_EMvar$dgamma)
		return (em.var)
	} else {
		return(-1)
	}
}

em.hmm.loadseed = function(EMvar, tmpdir)
{
	load(paste(paste(tmpdir, "seed_",sep="/"),EMvar$file,sep=""))
	return(seed_EMvar)
}


em.hmm.runseed = function(x, alttype, L, seed, burn, nulltype, maxiter, symmetric, ptol, core, v, tmpdir=tmpdir)
{
	EMvar = mclapply(1:seed, FUN = function(x, zval, alttype, L, burn, nulltype, symmetric, ptol, seed, v, tmpdir){
		seed_EMvar = list(logL=-Inf)
		while(length(seed_EMvar)<=2)
		{
			if(v) print(paste(paste(paste("seed : ",x),"/"),seed))
			seed_Evar     = em.hmm.init(zval, 0.95, 0.2, L, symmetric, v = v)
			seed_Mvar     = try(em.hmm.M(zval, seed_Evar, alttype, L, nulltype, symmetric, v = v))
			if( length(seed_Mvar) != 1)
			{
				rm(seed_Evar)
				seed_EMvar    = try(em.hmm.EM(zval, seed_Mvar, alttype, L, burn, nulltype, symmetric, ptol, E=F, v = v))
				rm(seed_Mvar)
			}
			else
			{
				if(v) print("Error in M")
			}
		}
		return( tryCatch({
			seed_file = paste(paste(tmpdir, "seed_",sep="/"),x,sep="")
			save(seed_EMvar, file=seed_file)
			return(list(logL = seed_EMvar$logL, file = x))
		}, error = function(e) {
			if(v) print(paste("error: in seed",x))
			return(list(logL = -Inf, file = NA))
		}))
	}, mc.cores = core, zval = x, alttype=alttype, L=L, burn=burn, nulltype=nulltype, symmetric=symmetric, ptol=ptol, seed=seed, v = v, tmpdir = tmpdir)
	
	logL = c()
	for(i in 1:seed)
	{
		logL[i] = ifelse(length(EMvar[[i]]) > 1, (EMvar[[i]])$logL, -Inf)
	}
	EMvar_ordered = list()
	j = 1
	for(i in order(logL, decreasing = T))
	{
		EMvar_ordered [[j]] = EMvar[[i]]
		j = j + 1
	}
	
	return(EMvar_ordered)
}



em.hmm.init = function(x, A11, A22, L = 2, symmetric, v)
{
#	if(v) print("init")
	NUM       = length(x)
	A         = array(0,c(2,2, NUM-1))
	
	
	tmp = NA
	while(sum(is.na(tmp))>0 & length(tmp) == 1)
	{
		A_11 = 1
		A_22 = 1
		while(A_11 + 1e-4 >= 1 | A_11 - 1e-4 <= 0)
		{
			A_11 = rbeta(1, A11, 1-A11)
		}
		while(A_22 + 1e-4 >= 1 | A_22 - 1e-4 <= 0)
		{
			A_22 = rbeta(1, A22, 1-A22)
		}
		A[1,1,] = A_11
		A[1,2,] = 1 - A[1,1,1]
		A[2,2,] = A_22
		A[2,1,] = 1 - A[2,2,1]
		tmp = try(inverse.rle( list(values=rep(1:2,NUM) ,lengths=1+rgeom( 2*NUM, rep( c( A[1,2,1], A[2,1,1] ), NUM) )))[1:NUM] - 1)
	}
	gamma = matrix(rep(NA, NUM*2), NUM, 2, byrow=TRUE)
	gamma[,1] = tmp
	gamma[,2] = 1-gamma[,1]
	
	gamma[(gamma[,1] == 1),1] = 0.9
	gamma[(gamma[,1] == 0),1] = 0.1
	gamma[(gamma[,2] == 1),2] = 0.9
	gamma[(gamma[,2] == 0),2] = 0.1
	
	if(symmetric)
	{
		L = L*2
	}
	omega     = t(rmultinom(NUM, L, rep(1/L, L)))
	omega     = omega[,]/L
	return(list(gamma = gamma, dgamma = c(A[1,1,1], A[2,2,1]), omega = omega))
}

em.hmm.EM = function(x, Mvar, alttype, L, maxiter, nulltype, symmetric, ptol, E, v)
{
	NUM = length(x)
	difference = 1
	converged=TRUE
	niter = 0
	logL.old = -Inf
	logL = 0
	
	Evar = list()
	Mvar.old = Mvar
	while(difference > ptol && niter <= maxiter)
	{
		niter = niter + 1
		Evar     = try(em.hmm.E(x, Mvar, alttype, L, symmetric, v = v))
		if( length(Evar) == 1)
		{
#			if(v) print("Error in E")
			converged=FALSE;
			break
		}
		
		logL  =  -sum(log(Evar$c0))
		
		if( logL < logL.old & abs(logL - logL.old) > 0.1)
		{
			if(v) print(paste(paste(paste("Error in EM : logL increasing by", logL.old - logL ), "iteration :"), niter))
			converged=FALSE;
			break
		}
		
		Mvar     = try(em.hmm.M(x, Evar, alttype, L, nulltype, symmetric, v = v))
		if( length(Mvar) == 1)
		{
#			if(v) print("Error in M")
			converged=FALSE;
			break
		}
		
		df1 = abs(Mvar.old$A - Mvar$A)
		df2 = abs(Mvar.old$f1 - Mvar$f1)
		df3  =  ifelse(abs(logL - logL.old) > 0.1, abs(logL - logL.old), 0)
		difference = max(df1, df2, df3)
		
		if( is.na(difference) )
		{
#			if(v) print("Error in EM : NA result")
			converged=FALSE;
			break
		}
		
		logL.old = logL
		Mvar.old = Mvar
		
	}
	if(!converged)
	{
		return(-1)
	}
	if(v) print(-sum(log(Evar$c0)))
	if(E)
	{
		return(list(pii=Mvar$pii, ptheta = Mvar$ptheta, pc=Mvar$pc, A=Mvar$A, f0=Mvar$f0, f1=Mvar$f1, logL = logL, gamma = Evar$gamma))
	}
	else
	{
		rm(Evar)
		return(list(pii=Mvar$pii, ptheta = Mvar$ptheta, pc=Mvar$pc, A=Mvar$A, f0=Mvar$f0, f1=Mvar$f1, logL = logL))
	}
}

em.hmm.E = function(x, Mvar, alttype, L, symmetric, v)
{
#	if(v) print("E step")
	res = list()
	NUM = length(x)
	delta = length(x[x==0])/length(x)
	f0x = c()
	f1x = c()
	omega = c()
	if(symmetric)
	{
		L = L*2
	}
	
	# f1x
	if(alttype == "kernel")
	{
		f1x = Mvar$f1
	}
	else
	{
		f1x = rep(0, NUM)
		if(L == 1)
		{
			f1x = dnorm(x, Mvar$f1[1], Mvar$f1[2])
		}
		else
		{
			for (c in 1:L)
			{
				f1x = f1x+Mvar$pc[c]*dnorm(x, Mvar$f1[c, 1], Mvar$f1[c, 2])
			}
		}
	}
	
	# f2x
	if(length(Mvar$f0) < 3)
	{
		f0x<-dnorm(x, Mvar$f0[1], Mvar$f0[2])
	}
	else
	{
		f0x = delta * (x==0) + (1-delta)*dnorm(x, Mvar$f0[1], Mvar$f0[2])# * (x!=0)
		f1x[x==0] = 0
	}
	
	# forward
	alpha<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	c0<-rep(0, NUM)
	alpha[1, 1]<-Mvar$pii[1]*f0x[1]
	alpha[1, 2]<-Mvar$pii[2]*f1x[1]
	c0[1]<-1/sum(alpha[1, ])
	alpha[1, ]<-c0[1]*alpha[1, ]
	alpha.tmp <- .C('calAlpha',alpha=as.numeric(alpha),c0=as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))
	alpha <- alpha.tmp$alpha
	dim(alpha) <- c(NUM,2)
	c0 <- alpha.tmp$c0
	
	# backward
	beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	beta[NUM, 1]<-c0[NUM]
	beta[NUM, 2]<-c0[NUM]
	beta.tmp <- .C('calBeta',beta=as.numeric(beta),as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))
	beta <- beta.tmp$beta
	dim(beta) <- c(NUM,2)
	
	# lfdr
	lfdr<-rep(0, NUM)
	lfdr.tmp <- .C('calLfdr',as.numeric(alpha),as.numeric(beta),lfdr=as.numeric(lfdr),as.integer(NUM))
	lfdr <- lfdr.tmp$lfdr
	
	# gamma & dgamma
	gamma<-matrix(1:(NUM*2), NUM, 2, byrow=TRUE)
	gamma[NUM, ]<-c(lfdr[NUM], 1-lfdr[NUM])
	dgamma<-array(rep(0, (NUM-1)*4), c(2, 2, (NUM-1)))
	gamma.tmp <- .C('calGamma',as.numeric(alpha),as.numeric(beta),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),gamma=as.numeric(gamma),dgamma=as.numeric(dgamma),as.integer(NUM))
	gamma <- gamma.tmp$gamma
	dgamma <- gamma.tmp$dgamma
	dim(gamma) <- c(NUM,2)
	dim(dgamma) <- c(2, 2, (NUM-1))
	
	if(alttype != "kernel" & L > 1)
	{
		# omega
		omega = matrix(rep(0, NUM*L), NUM, L, byrow=TRUE)
		for (c in 1:L)
		{ 
			f1c = dnorm(x, Mvar$f1[c, 1], Mvar$f1[c, 2])
			omega[, c] = gamma[, 2] * Mvar$pc[c]*f1c/f1x
			if(length(Mvar$f0) >= 3)
			{
				omega[x==0, c] = 0
			}
		}
	}
#	plot(res$pr[,1],col=color.scale(abs(x*2),c(1,1,0),0,c(0,1,1), alpha=1), pch="."); abline(h=c(0.95, median(Mvar$A[1,1,]), median(Mvar$A[2,2,1])));
	return(list(gamma = gamma, dgamma = dgamma, omega = omega, c0 = c0))
}

em.hmm.M = function(x, Evar, alttype, L, nulltype, symmetric, v)
{
	NUM = length(x)
	pii  =  c(0.95, 0.05)
	A  =  array(c(0.95, 0.05, 0.05, 0.95),c(2,2, NUM-1))
	f0 = c(0, 1)
	pc  =  rep(1, L)/L
	mus  =  rep(0, L)
	sds  =  rep(1, L)
	f1  =  cbind(mus, sds)
	ptheta  =  apply(Evar$gamma,2,sum)/NUM
	
	# transitions
	for (i in 1:2)
	{
		pii[i]  =  Evar$gamma[1, i]
	}
	if(is.vector(Evar$dgamma)) # if we are at the first step (i.e the init phase)
	{
		A[1,1,] = Evar$dgamma[1]
		A[1,2,] = 1 - A[1,1,1]
		A[2,2,] = Evar$dgamma[2]
		A[2,1,] = 1 - A[2,2,1]
	}
	else
	{
		for (i in 1:2)
		{
			for (j in 1:2)
			{ 
				q1  =  sum(Evar$dgamma[i, j, ])
				q2  =  sum(Evar$gamma[1:(NUM-1), i])
				A[i, j,]  =  q1/q2
			}
		}
	}
	
	
	# f0
	q5  =  sum(Evar$gamma[, 1]*x)
	mu0  =  q5/sum(Evar$gamma[, 1])
	q6  =  sum(Evar$gamma[, 1]*(x-mu0)*(x-mu0))
	sd0  =  sqrt(q6/sum(Evar$gamma[, 1]))
	f0 = c(mu0, sd0)
	
	if(nulltype == 0){
		f0  =  c(0,1)
	}
	if(nulltype == 1)
	{
		f0  =  c(0,1,-1)
	}
	
	# f1
	if(alttype == "kernel")
	{
		if(symmetric)
		{
			kern.f1  =  density(c(x,2*f0[1]-x),weights=c(Evar$gamma[,2],Evar$gamma[,2])/sum(c(Evar$gamma[,2],Evar$gamma[,2])))
			f1  =  approx(kern.f1$x, kern.f1$y, x, rule = 2, ties="ordered")$y
		}
		else
		{
			kern.f1  =  density(x,weights=Evar$gamma[,2]/sum(Evar$gamma[,2]))
			f1  =  approx(kern.f1$x, kern.f1$y, x, rule = 2, ties="ordered")$y
		}
		
		if(nulltype == 1)
		{
			f1[x == 0] = 0
		}
	}
	else
	{
		if(L == 1)
		{
			q1  =  sum(Evar$gamma[, 2])
			q2  =  sum(Evar$gamma[, 2]*x)
			mu1  =  q2/q1
			q3  =  sum(Evar$gamma[, 2]*(x-mu1)*(x-mu1))
			sd1  =  sqrt(q3/q1)
			f1  =  c(mu1, sd1)
		}
		else
		{
			mus  =  1:L
			sds  =  1:L
			for (c in 1:L)
			{
				q1  =  sum(Evar$omega[, c])
				q2  =  sum(Evar$gamma[, 2])
				pc[c]  =  q1/q2
				q3  =  sum(Evar$omega[, c]*x)
				mus[c]  =  q3/q1
				q4  =  sum(Evar$omega[, c]*(x-mus[c])*(x-mus[c]))
				sds[c]  =  sqrt(q4/q1)
			}
			f1  =  cbind(mus, sds)
		}
		if(symmetric & (L == 1 | L%%2 == 0))
		{
			f1_bis = f1
			f1_bis[,1] = -f1_bis[,1]
			f1 = rbind(f1, f1_bis)
			pc = c(pc,pc)
		}
	}
	return(list(pii = pii, ptheta = ptheta, pc = pc, A = A, f0 = f0, f1 = f1))
}
