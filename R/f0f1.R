f0f1 = function (zval, LIS, k=F)
{
	# plot(density(zval[zval != 0]))
	hist(zval[zval != 0], nclass=sqrt(length(zval)), main="", freq=F, ylim=c(0,1), xlab="z-value", cex.main = 3, cex.lab = 2, cex.axis = 2.5)
	x_tmp = order(zval)
	zval_tmp = zval[x_tmp]
	
	delta = length(zval[zval==0])/length(zval)
	f0 = c()
	f0 = (zval_tmp==0)*delta + (1-delta) * LIS$ptheta[1] * dnorm(zval_tmp, LIS$f0[1], LIS$f0[2])
	lines(zval_tmp, f0, col="blue", lwd=3)
	
	if(k)
	{
		lines(zval_tmp,LIS$f1[x_tmp]*LIS$ptheta[2], col="red", lwd=3)
	}
	else
	{
		f1 = rep(0, length(zval))
		if(length(LIS$f1) == 2)
		{
			f1 = LIS$ptheta[2] * dnorm(zval_tmp, LIS$f1[1], LIS$f1[2])
		}
		else
		{
			for(ell in 1:dim(LIS$f1)[1])
			{
				f1l = LIS$ptheta[2] * LIS$pc[ell] * dnorm(zval_tmp, LIS$f1[ell,1], LIS$f1[ell,2])
				# abline(v=LIS$f1[ell,1], col="red")
				f1 = f1 + f1l
			}
		}
		lines(zval_tmp, f1, col="red", lwd=3)
		lines(zval_tmp, f1 + f0, col="green", lwd=2)
	}
	# lines(zval_tmp, LIS$f1[x_tmp]*LIS$ptheta[2] + f0, col="green", lwd=2)
}
