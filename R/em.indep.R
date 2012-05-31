em.indep <-
function(x, alttype='mixnormal', L=2, maxiter=1000, nulltype = 2, symmetric = FALSE)
{

NUM<-length(x)

ptol<-1e-2
niter<-0

# Assuming it will converge 
converged=TRUE

### initializing model parameters
if(alttype == 'kernel'){

f0.new<-c(2, 1)

locfdr_p0 <- 0
if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1 | nulltype == 2)
{
	locfdr_p0 <- locfdr(x,plot=0)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

if(nulltype == 4)
{
	f0.new <- c(0,1,-1)
}
if(nulltype == 5)
{
	f0.new <- c(0,1,-1,0.5)
}

f1.new <- 0.5*dnorm(x,2,1)+0.5*dnorm(x,-2,1)

if(nulltype >= 4)
{
	f1.new[x == 0] = 0
}

ptheta.new <- c(0.95,0.05)
pZ.new <- c(0.5,0.5)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

ptheta.old <- ptheta.new
pZ.old <- pZ.new
f0.old <- f0.new
f1.old <- f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward1.indep.kernel(x, ptheta.old, pZ.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
Z <- forwardbackward.res$pr2
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

ptheta.new <- apply(gamma,2,sum)/NUM
pZ.new <- apply(Z,2,sum)/sum(gamma[, 1])

q5 <- sum(gamma[, 1]*x)
mu0 <- q5/sum(gamma[, 1])

q6 <- sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0 <- sqrt(q6/sum(gamma[, 1]))

f0.new<-c(mu0, sd0)

if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

if(nulltype == 4)
{
	f0.new <- c(0,1,-1)
}
if(nulltype == 5)
{
	q6bis <- sum( (Z[,1] - gamma[, 1])*x^2)
	sd0bis <- sqrt(q6/sum(Z[,1] - gamma[, 1]))
	f0.new <- c(0,1,-1,sd0bis)
}

if(symmetric == FALSE){
	kern.f1 <- density(x,weights=gamma[,2]/sum(gamma[,2]))
	f1.new <- approx(kern.f1$x, kern.f1$y, x, rule = 2, ties="ordered")$y
}

if(symmetric == TRUE){
	kern.f1 <- density(c(x,2*f0.new[1]-x),weights=c(gamma[,2],gamma[,2])/sum(c(gamma[,2],gamma[,2])))
	f1.new <- approx(kern.f1$x, kern.f1$y, x, rule = 2, ties="ordered")$y
}

if(nulltype >= 4)
{
	f1.new[x == 0] = 0
}

df2 <- abs(f1.old-f1.new)
diff <- max(df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

 lfdr <- gamma[, 1]
if (converged) {
 logL<- sum(log(c0))
if (nulltype > 0) {
	BIC <- logL-(3+2-2)*log(NUM)/2 
} else {
	BIC <- logL-(3-2)*log(NUM)/2 
}

 em.var<-list(ptheta=ptheta.new, f0=f0.new, f1=kern.f1, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged) 
} else {
 BIC<- logL<- (-Inf)
 em.var<-list(ptheta=ptheta.old, f0=f0.old, f1=kern.f1, LIS=lfdr, logL=logL, BIC=, ni=niter, converged=converged)
}

}

if(alttype == 'mixnormal'){

########
# L=1
########

if (L==1)
{

f0.new<-c(2, 1)

locfdr_p0 <- 0
if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1 | nulltype == 2)
{
	locfdr_p0 <- locfdr(x,plot=0)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

if(nulltype == 4)
{
	f0.new <- c(0,1,-1)
}
if(nulltype == 5)
{
	f0.new <- c(0,1,-1,0.5)
}

f1.new<-c(2, 1)
ptheta.new <- c(0.95,0.05)
pZ.new <- c(0.5, 0.5)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

ptheta.old <- ptheta.new
pZ.old <- pZ.new
f0.old <- f0.new
f1.old <- f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward1.indep(x, ptheta.old, pZ.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
Z <- forwardbackward.res$pr2
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

ptheta.new <- apply(gamma,2,sum)/NUM
pZ.new <- apply(Z,2,sum)/sum(gamma[, 1])

q5 <- sum(gamma[, 1]*x)
mu0 <- q5/sum(gamma[, 1])

q6 <- sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0 <- sqrt(q6/sum(gamma[, 1]))

f0.new<-c(mu0, sd0)

if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

if(nulltype == 4)
{
	f0.new <- c(0,1,-1)
}
if(nulltype == 5)
{
	q6bis <- sum( (Z[,1] - gamma[, 1])*x^2)
	sd0bis <- sqrt(q6/sum(Z[,1] - gamma[, 1]))
	f0.new <- c(0,1,-1,sd0bis)
}

q1 <- sum(gamma[, 2])
q2 <- sum(gamma[, 2]*x)
mu1 <- q2/q1
q3 <- sum(gamma[, 2]*(x-mu1)*(x-mu1))
sd1 <- sqrt(q3/q1)
f1.new <- c(mu1, sd1)

df2 <- abs(f1.old-f1.new)
diff <- max(df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

 lfdr <- gamma[, 1]
if (converged) {
 logL<- sum(log(c0))
if (nulltype > 0) {
	BIC<-logL-(3*L+2)*log(NUM)/2 
} else {
	BIC<-logL-(3*L)*log(NUM)/2 
}

 em.var<-list(ptheta=ptheta.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged) 
} else {
 BIC<- logL<- (-Inf)
 em.var<-list(ptheta=ptheta.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, BIC=, ni=niter, converged=converged)
}

}

#######
# L>1
#######

else if (L>1)
{

ptheta.new <- c(0.95, 0.05)
pZ.new <- c(0.5, 0.5)
pc.new <- rep(1, L)/L
mus <- seq(from=-1, by=1.5, length=L)
sds <- rep(1, L)
f0.new <- c(2, 1)

locfdr_p0 <- 0
if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1 | nulltype == 2)
{
	locfdr_p0 <- locfdr(x,plot=0)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

if(nulltype == 4)
{
	f0.new <- c(0,1,-1)
}
if(nulltype == 5)
{
	f0.new <- c(0,1,-1,0.5)
}

f1.new<-cbind(mus, sds)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter <- niter+1

ptheta.old <- ptheta.new
pZ.old <- pZ.new
pc.old <- pc.new
f0.old <- f0.new
f1.old <- f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward.indep(x, ptheta.old, pZ.old, pc.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
Z <- forwardbackward.res$pr2
omega <- forwardbackward.res$wt
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

ptheta.new <- apply(gamma,2,sum)/NUM
pZ.new <- apply(Z,2,sum)/sum(gamma[, 1])

q5 <- sum(gamma[, 1]*x)
mu0 <- q5/sum(gamma[, 1])

q6 <- sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0 <- sqrt(q6/sum(gamma[, 1]))

f0.new <- c(mu0, sd0)

if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

if(nulltype == 4)
{
	f0.new <- c(0,1,-1)
}
if(nulltype == 5)
{
	q6bis <- sum( (Z[,1] - gamma[, 1])*x^2)
	sd0bis <- sqrt(q6/sum(Z[,1] - gamma[, 1]))
	f0.new <- c(0,1,-1,sd0bis)
}


mus <- 1:L
sds <- 1:L

for (c in 1:L)
{

  q1 <- sum(omega[, c])
  q2 <- sum(gamma[, 2])
  pc.new[c] <- q1/q2
  
  q3 <- sum(omega[, c]*x)
  mus[c] <- q3/q1

  q4 <- sum(omega[, c]*(x-mus[c])*(x-mus[c]))
  sds[c] <- sqrt(q4/q1)

}

f1.new <- cbind(mus, sds)

df2 <- abs(f1.old-f1.new)
diff <- max(df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

lfdr <- gamma[, 1]
if (converged) {
 logL <- sum(log(c0))
if (nulltype > 0) {
	BIC <- logL-(3*L+2)*log(NUM)/2 
} else {
	BIC <- logL-(3*L)*log(NUM)/2 
}
 em.var<-list(ptheta=ptheta.new, pc=pc.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged) 
} else {
 logL <- (-Inf)
 BIC <- logL<- (-Inf)
 em.var <-list(ptheta=ptheta.old, pc=pc.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged)
}

}
}

return (em.var)
}

