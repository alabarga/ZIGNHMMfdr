
em.nhmm = function(x, Z, dist, dist.included=TRUE, alttype='mixnormal', L=2, maxiter=1000, nulltype=1, symmetric=FALSE, seed = 1000, burn = 20, iter.CG = 1000, ptol=1e-3, core = 2, v)
{
	NUM = length(x)
	
	if(length(Z)>0){
		if(is.vector(Z)==TRUE)
		{
			Z = matrix(Z,ncol=1)
			Z = scale(Z)
		}
		else
		{
			Z = apply(Z,2,scale)
		}
		if(dim(Z)[1] != NUM)
		{
			if(v) print("'Error: x and Z are not compatible")
			return(-1)
		}
		if(length(dist) > 0 & dim(Z)[1] != length(dist))
		{
			if(v) print("Error: dist and Z are not compatible")
			return(-1)
		}
	}
	Z = cbind(dist,Z)
	
	ptol = 1e-2
	difference = 1
	niter = 0
	EMvar = em.nhmm.runseed(x, Z, dist.included, alttype, L, seed, burn, nulltype, maxiter, symmetric, iter.CG, ptol, core, v = v)
	best_EMvar = list()
	while(length(best_EMvar) <= 1)
	{
		niter = niter + 1
		best_EMvar = EMvar[[niter]]
		if(v) print(paste("best seed : ",best_EMvar$logL))
		best_EMvar = try(em.nhmm.EM(x, Z, dist.included, best_EMvar, alttype, L, maxiter, nulltype, symmetric, iter.CG, ptol, E = T, v = v))
		if(length(EMvar) <= 1)
		{
			if(v) print("error: Trying next best seed")
			best_EMvar = EMvar[[niter]]
		}
	}
	
	lfdr = best_EMvar$gamma[, 1]
	if(length(best_EMvar) > 1)
	{
		logL  =  best_EMvar$logL
		if(v) print(logL)
		if (nulltype > 0) {
			BIC  =  logL - (3*L + dim(Z)[2] + 2 + 1)*log(NUM)/2 
		} else {
			BIC  =  logL - (3*L + dim(Z)[2] + 2)*log(NUM)/2 
		}
		em.var = list(ptheta = best_EMvar$ptheta, pii = best_EMvar$pii, A = best_EMvar$A, pc = best_EMvar$pc, f0 = best_EMvar$f0, f1 = best_EMvar$f1, LIS=lfdr, logL=logL, BIC=BIC, trans.par = best_EMvar$trans.par, gamma = best_EMvar$gamma, dgamma = best_EMvar$dgamma)
		return (em.var)
	} else {
		return(-1)
	}
}

em.nhmm.runseed = function(x, Z, dist.included, alttype, L, seed, burn, nulltype, maxiter, symmetric, iter.CG, ptol, core, v)
{
	EMvar = mclapply(1:seed, FUN = function(x, zval, alttype, L, burn, nulltype, symmetric, ptol, seed, v){
		seed_EMvar = list()
		while(length(seed_EMvar)<=1)
		{
			if(v) print(paste(paste(paste("seed : ",x),"/"),seed))
			seed_Evar     = em.nhmm.init(zval, Z, dist.included, 0.95, 0.2, L, symmetric, v = v)
			seed_Mvar     = try(em.nhmm.M(zval, Z, dist.included, seed_Evar, alttype, L, nulltype, symmetric, iter.CG, ptol, v = v))
			if( length(seed_Mvar) == 1)
			{
				if(v) print("Error in M")
				converged=FALSE;
				break
			}
			rm(seed_Evar)
			seed_EMvar    = try(em.nhmm.EM(zval, Z, dist.included, seed_Mvar, alttype, L, burn, nulltype, symmetric, iter.CG, ptol, E = F, v = v))
			rm(seed_Mvar)
		}
		return(seed_EMvar)
	}, mc.cores = core, zval = x, alttype=alttype, L=L, burn=burn, nulltype=nulltype, symmetric=symmetric, ptol=ptol, seed=seed, v = v)

	
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

em.nhmm.init = function(x, Z, dist.included, A11, A22, L, symmetric, v)
{
	NUM = length(x)
	
	A         = array(0,c(2,2, NUM-1))
	trans.par = array(0,c(2,3+dim(Z)[2]))
	
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
	
	dgamma    = array(0,c(2,2, NUM-1))
	dgamma[1,1,] = as.numeric(gamma[-dim(gamma)[1],2] == 0 & gamma[-1,1] == 0)
	dgamma[2,1,] = as.numeric(gamma[-dim(gamma)[1],2] == 1 & gamma[-1,1] == 0)
	dgamma[1,2,] = as.numeric(gamma[-dim(gamma)[1],2] == 0 & gamma[-1,2] == 1)
	dgamma[2,2,] = as.numeric(gamma[-dim(gamma)[1],2] == 1 & gamma[-1,2] == 1)
	
	for(i in 1:2)
	{
		tmp = glm(dgamma[i,2,] ~ Z[-1,], family=binomial("logit"))
		trans.par[2,1+i]     = tmp$coefficient[1]
		trans.par[2,-c(1:3)] = tmp$coefficient[-1]
	}
	if(dist.included)
	{
		trans.par[2,4] = abs(trans.par[2,4])
	}
	
	if(symmetric)
	{
		L = L*2
	}
	omega     = t(rmultinom(NUM, L, rep(1/L, L)))
	omega     = omega[,]/L
	
	return(list(gamma = gamma, dgamma = dgamma, omega = omega, trans.par = trans.par))
}

em.nhmm.EM = function(x, Z, dist.included, Mvar, alttype, L, maxiter, nulltype, symmetric, iter.CG, ptol, E, v)
{
	NUM = length(x)
	difference = 1
	converged=TRUE
	niter = 0
	logL.old = -Inf
	logL = 0
	T = 1
	Evar = list()
	Mvar.old = Mvar
	
	while(difference > ptol && niter <= maxiter)
	{
		niter = niter + 1
		Evar     = try(em.nhmm.E(x, Mvar.old, alttype, L, symmetric, v = v))
		if( length(Evar) == 1)
		{
			if(v) print("Error in E")
			converged=FALSE;
			break
		}
		
		logL  =  -sum(log(Evar$c0))
		
		if( logL < logL.old & abs(logL - logL.old) > 0.1)
		{
			if(v) print(paste(paste(paste("Error in EM : logL increasing by", logL.old - logL ), "iteration :"), niter))
			converged=FALSE
			break
		}
		
		Mvar     = try(em.nhmm.M(x, Z, dist.included, Evar, alttype, L, nulltype, symmetric, iter.CG, ptol, v = v))
		if( length(Mvar) == 1)
		{
			if(v) print("Error in M")
			converged=FALSE
			break
		}
		
		df1 = abs(Mvar.old$trans.par[2,-1] - Mvar$trans.par[2,-1])
		df2 = abs(Mvar.old$f1 - Mvar$f1)
		df3  =  ifelse(abs(logL - logL.old) > 0.1, abs(logL - logL.old), 0)
		difference = max(df1, df2, df3)
		
		if( is.na(difference) )
		{
			if(v) print("Error in EM : NA result")
			converged=FALSE
			break
		}
		logL.old = logL
		Mvar.old = Mvar
	}
	if(!converged)
	{
		return(-1)
	}
	if(v) print(logL.old)
	if(E)
	{
		return(list(pii=Mvar.old$pii, ptheta = Mvar.old$ptheta, pc=Mvar.old$pc, A=Mvar.old$A, trans.par = Mvar.old$trans.par, f0=Mvar.old$f0, f1=Mvar.old$f1, logL = logL.old, gamma = Evar$gamma))
	}
	else
	{
		rm(Evar)
		return(list(pii=Mvar.old$pii, ptheta = Mvar.old$ptheta, pc=Mvar.old$pc, A=Mvar.old$A, trans.par = Mvar.old$trans.par, f0=Mvar.old$f0, f1=Mvar.old$f1, logL = logL.old))
	}
}


em.nhmm.E = function(x, Mvar, alttype, L, symmetric, v)
{
	if(symmetric)
	{
		L = L*2
	}
	
	res = list()
	if(alttype == "kernel")
	{
		res  =  forwardbackward1.kernel(x, Mvar$pii, Mvar$A, Mvar$f0, Mvar$f1)
		res$wt = NA
	}
	else
	{
		if(L == 1)
		{
			forwardbackward.res  =  forwardbackward1(x, Mvar$pii, Mvar$A, Mvar$f0, Mvar$f1)
			res$wt = NA
		}
		else
		{
			res  =  forwardbackward(x, Mvar$pii, Mvar$A, Mvar$pc, Mvar$f0, Mvar$f1)
		}
		
	}
	return(list(gamma = res$pr, dgamma = res$ts, omega = res$wt, c0 = res$rescale, trans.par = Mvar$trans.par))
}

#em.nhmm.E = function(x, Mvar, alttype, L, symmetric, v)
#{
#	res = list()
#	NUM = length(x)
#	delta = length(x[x==0])/length(x)
#	f0x = c()
#	f1x = c()
#	omega = c()
#	if(symmetric)
#	{
#		L = L*2
#	}
#	
#	# f1x
#	if(alttype == "kernel")
#	{
#		f1x = Mvar$f1
#	}
#	else
#	{
#		f1x = rep(0, NUM)
#		if(L == 1)
#		{
#			f1x = dnorm(x, Mvar$f1[1], Mvar$f1[2])
#		}
#		else
#		{
#			for (c in 1:L)
#			{
#				f1x = f1x+Mvar$pc[c]*dnorm(x, Mvar$f1[c, 1], Mvar$f1[c, 2])
#			}
#		}
#	}
#	
#	# f0x
#	if(length(Mvar$f0) < 3)
#	{
#		f0x=dnorm(x, Mvar$f0[1], Mvar$f0[2])
#	}
#	else
#	{
#		f0x = delta * (x==0) + dnorm(x, Mvar$f0[1], Mvar$f0[2]) * (x!=0)
#		f1x[x==0] = 0
#	}
#	
#	# forward
#	alpha=matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
#	c0=rep(0, NUM)
#	alpha[1, 1]=Mvar$pii[1]*f0x[1]
#	alpha[1, 2]=Mvar$pii[2]*f1x[1]
#	c0[1]=1/sum(alpha[1, ])
#	alpha[1, ]=c0[1]*alpha[1, ]
#	alpha.tmp = .C('calAlpha',alpha=as.numeric(alpha),c0=as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))
#	alpha = alpha.tmp$alpha
#	dim(alpha) = c(NUM,2)
#	c0 = alpha.tmp$c0
#	
#	# backward
#	beta=matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
#	beta[NUM, 1]=c0[NUM]
#	beta[NUM, 2]=c0[NUM]
#	beta.tmp = .C('calBeta',beta=as.numeric(beta),as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))
#	beta = beta.tmp$beta
#	dim(beta) = c(NUM,2)
#	
#	# lfdr
#	lfdr=rep(0, NUM)
#	lfdr.tmp = .C('calLfdr',as.numeric(alpha),as.numeric(beta),lfdr=as.numeric(lfdr),as.integer(NUM))
#	lfdr = lfdr.tmp$lfdr
#	
#	# gamma & dgamma
#	gamma=matrix(1:(NUM*2), NUM, 2, byrow=TRUE)
#	gamma[NUM, ]=c(lfdr[NUM], 1-lfdr[NUM])
#	dgamma=array(rep(0, (NUM-1)*4), c(2, 2, (NUM-1)))
#	gamma.tmp = .C('calGamma',as.numeric(alpha),as.numeric(beta),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),gamma=as.numeric(gamma),dgamma=as.numeric(dgamma),as.integer(NUM))
#	gamma = gamma.tmp$gamma
#	dgamma = gamma.tmp$dgamma
#	dim(gamma) = c(NUM,2)
#	dim(dgamma) = c(2, 2, (NUM-1))
#	
#	if(alttype != "kernel" & L > 1)
#	{
#		# omega
#		omega = matrix(rep(0, NUM*L), NUM, L, byrow=TRUE)
#		for (c in 1:L)
#		{ 
#			f1c = dnorm(x, Mvar$f1[c, 1], Mvar$f1[c, 2])
#			omega[, c] = gamma[, 2] * Mvar$pc[c]*f1c/f1x
#			if(length(Mvar$f0) >= 3)
#			{
#				omega[x==0, c] = 0
#			}
#		}
#	}
#	return(list(gamma = gamma, dgamma = dgamma, omega = omega, c0 = c0, trans.par = Mvar$trans.par))
#}

em.nhmm.M = function(x, Z, dist.included, Evar, alttype, L, nulltype, symmetric, iter.CG, ptol, v)
{
#	if(v) print("M step")
	NUM = length(x)
	f0 = c(0, 1)
	pc  =  rep(1, L)/L
	mus  =  rep(0, L)
	sds  =  rep(1, L)
	f1  =  cbind(mus, sds)
	ptheta  =  apply(Evar$gamma,2,sum)/NUM
	CG = list()
	
	# f0
	q5 = sum(Evar$gamma[, 1]*x)
	mu0 = q5/sum(Evar$gamma[, 1])
	
	q6 = sum(Evar$gamma[, 1]*(x-mu0)*(x-mu0))
	sd0 = sqrt(q6/sum(Evar$gamma[, 1]))
	
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
		if(symmetric == FALSE){
			kern.f1  =  density(x,weights=Evar$gamma[,2]/sum(Evar$gamma[,2]))
			f1  =  approx(kern.f1$x, kern.f1$y, x, rule = 2, ties="ordered")$y
		}
		
		if(symmetric == TRUE){
			kern.f1  =  density(c(x,2*f0[1]-x),weights=c(Evar$gamma[,2],Evar$gamma[,2])/sum(c(Evar$gamma[,2],Evar$gamma[,2])))
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
#	CG = try(em.nhmm.compute.CG(Z, dist.included, Evar$dgamma, Evar$gamma, Evar$trans.par, iter.CG = 10, ptol, v = v))
	CG = try(update.trans.prob.nhmm(Z, Evar$dgamma, Evar$gamma, Evar$trans.par[1,], Evar$trans.par[2,], iter.conj.grad = 10, dist.included))
	if(length(CG) == 1)
	{
		if(v) print("Error in CG")
		return(-1)
	}
#	return(list(pii = CG$pii, ptheta = ptheta, pc=pc, A = CG$A , trans.par = CG$trans.par, f0 = f0, f1 = f1))
	return(list(pii = CG$pii, ptheta = ptheta, pc=pc, A = CG$A , trans.par = rbind(CG$trans.par1, CG$trans.par2), f0 = f0, f1 = f1))
}

em.nhmm.compute.CG = function(Z, dist.included, dgamma, gamma, trans.par, iter.CG, ptol, v)
{
	gradient.old = em.nhmm.compute.gradient(Z, dist.included, dgamma, gamma, trans.par, v = v)
	trans.par.old = trans.par
	phi = -gradient.old
	difference = 1
	niter = 0
	while(difference > ptol & niter < iter.CG)
	{
		trans.par.old = trans.par
		niter = niter + 1
		tmp = try( em.nhmm.line.search(Z, dist.included, dgamma, gamma, trans.par, phi, iter.CG, ptol, v = v) )
		if(is.na(tmp$nu))
		{
			break
		}
		trans.par = tmp$trans.par
		gradient.new = em.nhmm.compute.gradient(Z, dist.included, dgamma, gamma, trans.par, v = v)
		PR = sum((gradient.new - gradient.old)*gradient.new) / sum(gradient.old^2)
		if(is.nan(PR) | PR < 0)
		{
			PR = 0
		}
		phi = -gradient.new + PR * phi
		gradient.old = gradient.new
		if(dist.included & trans.par[2,4]<0 )
		{
			trans.par[2,4] = abs(trans.par[2,4])
		}
		difference = max(abs(trans.par[2,-1] - trans.par.old[2,-1]))
	}
	if(!is.na(tmp$nu))
	{
		tmp = em.nhmm.pii.A(Z, dist.included, trans.par, v = v)
		return(list(pii = tmp$pii, A = tmp$A, trans.par = trans.par))
	}else
	{
		return(-1)
	}
}

em.nhmm.line.search = function(Z, dist.included, dgamma, gamma, trans.par, phi, iter.CG, ptol, v)
{
	N = dim(Z)[1] - 1
	nu = 0
	nu.old = 1
	niter = 0
	difference = 1
	phi[1,] = 0
	
	trans.par.new = trans.par
	
	fixe_1  = phi[,1] + sum( phi[,-c(1:3)]  * Z[1,] )
	fixe_2  = array(0,c(2,2,N))
	for(i in 1:2)
	{
		for(j in 1:2)
		{
			if(dist.included & i == 2)
				fixe_2[i,j,] = ( phi[j,i+1] + rowSums( t(c(-phi[j,4],phi[j,-c(1:4)]) * t(Z[-1,])) ) )
			else
				fixe_2[i,j,] = ( phi[j,i+1] + rowSums( t(phi[j,-c(1:3)] * t(Z[-1,])) ) )
		}
	}
	
	while(difference > ptol & niter < iter.CG)
	{
		niter = niter + 1
		trans.par.new = trans.par + nu * phi
		trans.prob = em.nhmm.pii.A(Z, dist.included, trans.par.new, v = v)
		dQ  =  sum( fixe_1   * (gamma[1,] - trans.prob$pii) )
		dQ2 = -sum( fixe_1^2 * trans.prob$pii * (1 - trans.prob$pii) )
		
		for(i in 1:2)
		{
			for(j in 1:2)
			{
				dQ  = dQ  + sum( fixe_2[i,j,]   * ( dgamma[i,j,] - ( gamma[-dim(gamma)[1],i] * trans.prob$A[i,j,] ) ) )
				dQ2 = dQ2 - sum( fixe_2[i,j,]^2 * trans.prob$A[i,j,] * ( 1 - trans.prob$A[i,j,] ) * gamma[-dim(gamma)[1], i] )
			}
		}
		
		nu = nu - dQ / dQ2
		difference = abs(dQ/ dQ2)
		if(is.na(difference) | niter > 100){
			nu = NaN
			break
		}
	}
	if(is.nan(dQ2))
	{
		if(v) print("Error in line search")
	}
	trans.par.new = trans.par + nu * phi
	return(list(nu = nu, trans.par = trans.par.new))
}

em.nhmm.compute.gradient = function(Z, dist.included, dgamma, gamma, trans.par, v)
{
	N = dim(Z)[1] - 1
	
	gradient      = array(0,c(2,3+dim(Z)[2]))
	
	trans.prob    = em.nhmm.pii.A(Z, dist.included, trans.par, v = v)
	
	gradient[2,1] = gamma[1,2] - trans.prob$pii[2]
	gradient[2,2] = sum( dgamma[1,2,] - (gamma[-dim(gamma)[1],1]*trans.prob$A[1,2,]) )
	gradient[2,3] = sum( dgamma[2,2,] - (gamma[-dim(gamma)[1],2]*trans.prob$A[2,2,]) )
	
	gradient[2, -c(1:3)] = (gamma[1,2] - trans.prob$pii[1]) * Z[1,]
	gradient[2, -c(1:3)] = (gradient[2, -c(1:3)] + colSums(( dgamma[1,2,] - (gamma[-dim(gamma)[1],1] * trans.prob$A[1,1,] ) ) * Z[-1,])) + (gradient[2, -c(1:3)] + colSums(( dgamma[2,2,] - (gamma[-dim(gamma)[1],2] * trans.prob$A[2,1,] ) ) * Z[-1,]))
	
	return(-gradient)
}

em.nhmm.pii.A = function(Z, dist.included, trans.par, v)
{
	N = dim(Z)[1] - 1
	pii = rep(0,2)
	A = array(0,c(2,2,N))
	
	pii[1] = exp( trans.par[1,1] + sum(trans.par[1,-c(1:3)] * Z[1,]) ) / ( exp( trans.par[1,1] + sum(trans.par[1,-c(1:3)] * Z[1,]) ) + exp( trans.par[2,1] + sum(trans.par[2,-c(1:3)] * Z[1,]) ) )
	pii[2] = 1 - pii[1]
	
	rho_Z_2    = rowSums(t(trans.par[2,-c(1:3)]*t(Z[-1,])))
	rho_Z_2bis = rowSums(t(trans.par[2,-c(1:3)]*t(Z[-1,])))
	if(dist.included)
	{
		rho_Z_2    = rowSums(t(c( trans.par[2,4], trans.par[2,-c(1:4)]) * t(Z[-1,])))
		rho_Z_2bis = rowSums(t(c(-trans.par[2,4], trans.par[2,-c(1:4)]) * t(Z[-1,])))
	}
	
	q1 = exp(trans.par[1,2]) + exp( trans.par[2,2] + rho_Z_2 )
	A[1,1,] = exp( trans.par[1,2] )              / q1
	A[1,2,] = exp( trans.par[2,2] + rho_Z_2 )    / q1
	q1 = exp(trans.par[1,3]) + exp( trans.par[2,3] + rho_Z_2bis )
	A[2,1,] = exp( trans.par[1,3] )              / q1
	A[2,2,] = exp( trans.par[2,3] + rho_Z_2bis ) / q1
	
	return(list(pii = pii, A = A))
}

