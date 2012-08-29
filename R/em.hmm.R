
em.hmm = function(x, alttype="mixnormal", L=2, maxiter=1000, nulltype=2, symmetric= FALSE, seed = 100)
{
	NUM = length(x)
	ptol = 1e-2
	diff = 1
	niter = 0
	converged = TRUE
	
	Mvar = em.hmm.runseed(x, alttype, L, seed, nulltype, symmetric)
	logL.iter = Mvar$logL.iter
	Mvar.old = Mvar
	
	while(diff>ptol && niter<maxiter)
	{
		Evar     = em.hmm.E(x, Mvar, alttype, L)
		Mvar     = em.hmm.M(x, Evar, alttype, L, nulltype, symmetric)
		
		logL.iter  =  c(logL.iter,-sum(log(Evar$c0)))
		
		df1 = abs(Mvar.old$A - Mvar$A)
		df2 = abs(Mvar.old$f1 - Mvar$f1)
		diff = max(df1, df2)
		
		if (is.na(diff)) {
			print("Error in EM")
			converged=FALSE;
			break;
		}
		niter = niter + 1
		Mvar.old = Mvar
	}
	
	lfdr = Evar$gamma[, 1]
	if(converged)
	{
		logL  =  -sum(log(Evar$c0))
		if (nulltype > 0) {
			BIC = logL-(3*L+3)*log(NUM)/2 
		} else {
			BIC = logL-(3*L+2)*log(NUM)/2 
		}
		em.var  =  list(ptheta = Mvar$ptheta, pii=Mvar$pii, A=Mvar$A, pc=Mvar$pc, f0=Mvar$f0, f1=Mvar$f1, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged, gamma = Evar$gamma, dgamma = Evar$dgamma, logL.iter = logL.iter) 
	} else {
		logL  =  (-Inf)
		BIC  =  logL  =  (-Inf)
		em.var  =  list(ptheta = Mvar.old$ptheta, pii=Mvar.old$pii, A=Mvar.old$A, pc=Mvar.old$pc, f0=Mvar.old$f0, f1=Mvar.old$f1, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged, gamma = Evar$gamma, dgamma = Evar$dgamma, logL.iter = logL.iter)
	}
	return (em.var)
}

em.hmm.runseed = function(x, alttype="mixnormal", L=2, seed=100, nulltype=2, symmetric= FALSE)
{
	niter = 0
	seed_result = list()
	seed_logL = c()
	Mvar = list(logL = -Inf)
	
	while(niter <= seed)
	{
		tmp  =  try(em.hmm.seed(x, alttype, L, nulltype, symmetric))
		if(length(tmp)==1)
		{
			print("Numerical ERROR: Rerunning...")
		}
		else
		{
			if(!is.na(tmp$logL))
			{
				if(tmp$logL > Mvar$logL)
				{
					Mvar = tmp
				}
				niter = niter + 1
			}
		}
	}
	return(list(pii=Mvar$pii, ptheta = Mvar$ptheta, pc=Mvar$pc, A=Mvar$A, f0=Mvar$f0, f1=Mvar$f1, logL.iter = Mvar$logL.iter))
}

em.hmm.seed = function(x, alttype="mixnormal", L=2, nulltype=2, symmetric= FALSE)
{
	NUM = length(x)
	ptol = 1e-2
	diff = 1
	converged=TRUE
	niter = 0
	logL.iter = c()
	logL = NA
	
	Evar = em.hmm.init(NUM, 0.95, 0.2, L)
	Mvar = em.hmm.M(x, Evar, alttype, L, nulltype, symmetric)
	Mvar.old = Mvar
	
	while(diff>ptol && niter <= 10)
	{
		Evar     = em.hmm.E(x, Mvar, alttype, L)
		Mvar     = em.hmm.M(x, Evar, alttype, L, nulltype, symmetric)
		
		logL.iter  =  c(logL.iter,-sum(log(Evar$c0)))
		
		df1 = abs(Mvar.old$A - Mvar$A)
		df2 = abs(Mvar.old$f1 - Mvar$f1)
		diff = max(df1, df2)
		
		if (is.na(diff)) {
			print("Error in EM")
			converged=FALSE;
			break;
		}
		niter = niter + 1
		Mvar.old = Mvar
	}
	if(converged)
	{
		logL = -sum(log(Evar$c0))
	}
	else
	{
		logL = NA
	}
	return(list(pii=Mvar$pii, ptheta = Mvar$ptheta, pc=Mvar$pc, A=Mvar$A, f0=Mvar$f0, f1=Mvar$f1, logL.iter = logL.iter, logL = logL))
}

em.hmm.init = function(NUM, a11, a22, L = 2)
{
	a  =  array(0,c(2,2, NUM-1))
	a_11 = 1
	a_22 = 1
	while(a_11 == 1 | a_11 == 0 | a_22 == 1 | a_22 == 0)
	{
		a_11 = rbeta(1, a11, 1-a11)
		a_22 = rbeta(1, a22, 1-a22)
	}
	a[1,1,] = a_11
	a[1,2,] = 1 - a[1,1,1]
	a[2,2,] = a_22
	a[2,1,] = 1 - a[2,2,1]
	
	gamma = matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	omega = matrix(rep(0, NUM*L), NUM, L, byrow=TRUE)
	
	while(sum(is.na(gamma[,1]))>0)
	{
		gamma[,1] = inverse.rle( list(values=rep(1:2,NUM) ,lengths=1+rgeom( 2*NUM, rep( c( a[1,2,1], a[2,1,1] ), NUM) )))[1:NUM] - 1
	}
	gamma[,2] = 1-gamma[,1]
	
	gamma[(gamma[,1] == 1),1] = 0.9
	gamma[(gamma[,1] == 0),1] = 0.1
	gamma[(gamma[,2] == 1),2] = 0.9
	gamma[(gamma[,2] == 0),2] = 0.1
	
	omega     = t(rmultinom(NUM, L, rep(1/L, L)))
	omega     = omega[,]/L
	return(list(gamma = gamma, dgamma = c(a[1,1,1], a[2,2,1]), omega = omega))
}

em.hmm.E = function(x, Mvar, alttype, L)
{
	res = list()
	if(alttype == "kernel")
	{
		res = forwardbackward1.kernel(x, Mvar$pii, Mvar$A, Mvar$f0, Mvar$f1)
		res$wt = NA
	}
	else
	{
		if(L == 1)
		{
			res = forwardbackward1(x, Mvar$pii, Mvar$A, Mvar$f0, Mvar$f1)
			res$wt = NA
		}
		else
		{
			res = forwardbackward(x, Mvar$pii, Mvar$A, Mvar$pc, Mvar$f0, Mvar$f1)
		}
	}
	
	return(list(gamma = res$pr, dgamma = res$ts, c0 = res$rescale, omega = res$wt))
}

em.hmm.M = function(x, Evar, alttype, L, nulltype, symmetric)
{
	NUM = length(x)
	pii  =  c(0.95, 0.05)
	ptheta  =  c(0.95,0.05)
	A  =  array(c(0.95, 0.05, 0.05, 0.95),c(2,2, NUM-1))
	f0 = c(0, 1)
	pc  =  rep(1, L)/L
	mus  =  rep(0, L)
	sds  =  rep(1, L)
	f1  =  cbind(mus, sds)
	
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
	
	ptheta  =  apply(Evar$gamma,2,sum)/NUM
	
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
			if(symmetric & L%%2 == 0)
			{
				tmp      = order(mus)
				mus      = mus[tmp]
				sds      = sds[tmp]
				if(abs(mus[1]) < abs(mus[L]))
				{
					mus      = c(-mus[(L/2+1):L], mus[(L/2+1):L])
				}
				else
				{
					mus      = c(mus[1:(L/2)], abs(mus[1:(L/2)]))
				}
			}
			f1  =  cbind(mus, sds)
		}
	}
	return(list(pii = pii, ptheta = ptheta, pc = pc, A = A, f0 = f0, f1 = f1))
}
