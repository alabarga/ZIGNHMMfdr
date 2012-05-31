forwardbackward1.indep.kernel <-
function(x, ptheta, pZ, f0, f1x)
{


## Initialize

NUM<-length(x)
delta = length(x[x==0])/length(x)
## Densities

f0x <- c()
if(length(f0) < 4)
{
	f0x<-dnorm(x, f0[1], f0[2])
}
else
{
	if(length(f0) == 4)
	{
		f0x<- delta * (x==0) + (1-delta)*dnorm(x, f0[1], f0[2]) * (x!=0)
	}
	if(length(f0) == 5)
	{
		f0x<- delta * (x==0) + (1-delta)* ( pZ[1] * dnorm(x, f0[1], f0[2]) + (1-pZ[1]) * dnorm(x, f0[1], f0[4])) * (x!=0)
	}
	f1x[x==0] = 0
}



gamma<-matrix(1:(NUM*2), NUM, 2, byrow=TRUE)

gamma[,1] <- ptheta[1]*f0x/(ptheta[1]*f0x + ptheta[2]*f1x)
gamma[,2] <- 1 - gamma[,1]

Z<-matrix(1:(NUM*2), NUM, 2, byrow=TRUE)

if(length(f0) == 5)
{
	Z[,1] <- pZ[1]*dnorm(x, f0[1], f0[2])/f0x * gamma[,1] * (1-delta)
	Z[,2] <- 1 - Z[,1]
}

lfdr <- gamma[,1]

c0 <- f0x*ptheta[1] + f1x*ptheta[2]

forwardbackward.var<-list(lf=lfdr, pr=gamma, pr2 = Z,rescale=c0)
return(forwardbackward.var)
  
}

