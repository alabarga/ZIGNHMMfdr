simdata.nhmm <-
function(NUM, pii, A, f0, pc, f1, alpha)
{


theta<-rep(0, NUM)
x<-rep(0, NUM)
nc<-length(pc)

if (nc==1)
{ 
  data <- simdata1.nhmm (NUM, pii, A, f0, f1, alpha)
}    

else
{
## generating the states
 # initial state
theta[1]<-rbinom(1, 1, pii[2])
 # other states
for (i in 2:NUM)
{
  if (theta[i-1]==0)
     theta[i]<-rbinom(1, 1, A[1, 2, i-1])
  else
     theta[i]<-rbinom(1, 1, A[2, 2, i-1])
}

x_f0      = c( rep(0, floor(NUM*alpha)), rnorm(floor(NUM*(1-alpha)), f0[1], f0[2]))
x_f1      = c()
for(ell in 1:nc)
{
	x_f1 = c(x_f1, rnorm(floor(NUM*pc[ell]), f1[ell,1], f1[ell,2]))
}

x_f0 = sample(x_f0, NUM, replace = T)
x_f1 = sample(x_f1, NUM, replace = T)

x = ifelse(theta == 0, x_f0, x_f1)

data<-list(s=theta, o=x)
}

return (data)

}

