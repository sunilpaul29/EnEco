######################################################################
##Function 1 lags
makelag<- function(y,nlag) 
{
  #arg y is a matrix
  # nlag number of lages
  nr=nrow(y)
  nc=ncol(y)
  y1 <- y[1:(nr-nlag),]
  ylag <- rbind(matrix(0,nrow=nlag,ncol=nc),y1)
  return(ylag)
} 
######################################################################
# Wishart distribution
wish<-function(h,n){
    m <- nrow(h)
    A <- t(chol(h))%*% matrix( rnorm(m*n,mean=0,sd=1), m, n) 
    A <- A%*%t(A)
}

######################################################################
# Function 2  Priors

ts_prior <- function(rawdat,tau,p,plag){
yt <- t(rawdat[((plag+1):(tau+plag)),])
m  <- p + (plag*(p^2)) # number of state vectors
Zt <-matrix(ncol=m) # Zt
for (i in (plag+1):(tau+plag))
{
    ztemp <- diag(1,p)
    for (j in 1:plag)
        {       
        xlag = rawdat[(i-j),(1:p)]
        xtemp = matrix(0,nrow=p, ncol=p*p)
        for(jj in 1:p){
         xtemp[jj,(((jj-1)*p)+1):(jj*p)] = xlag  
        }
    ztemp <- cbind(ztemp,xtemp)
    }
    Zt <-rbind(Zt,ztemp)
}

Zt <-Zt[-1,]

   
# see the main document for details
vbar <- matrix(0,nrow=m,ncol=m) 
xhy <- matrix(0,ncol=1,nrow=m) 

for (i in 1:tau){
    zhat1 <- Zt[(((i-1)*p+1):(i*p)),]
    vbar <- vbar + (t(zhat1)%*%zhat1)
    xhy <- xhy + t(zhat1)%*%yt[,i]
    }

vbar <- solve(vbar)
aols  <- vbar%*%xhy

sse2 <- matrix(0,nrow=p,ncol=p)

for (i in 1:tau){
    zhat1 <- Zt[((i-1)*p+1):(i*p),]
    sse2 = sse2 + (yt[,i] - (zhat1%*%aols))%*%t(yt[,i] - zhat1%*%aols)
}

hbar = sse2/tau # element wise right division

vbar <- matrix(0,ncol=m, nrow=m)

for (i in 1 : tau){
    zhat1 = Zt[((i-1)*p+1):(i*p),]
    vbar = vbar + t(zhat1)%*%solve(hbar)%*%zhat1
    }

library(Matrix) # for chol()

vbar <- solve(vbar)
achol <- t(chol(hbar))

ssig  <- matrix(0,nrow=p,ncol=p)


for( i in 1:p){
    ssig[i,i] = achol[i,i] 
    for (j in 1:p){
        achol[j,i] = achol[j,i]/ssig[i,i]
 }
}

achol  <- solve(achol)
numa <-  p*(p-1)/2
a0 <- matrix(0,nrow=numa,ncol=1)

ic = 1

for (i in 2:p){
    for (j in 1:(i-1)){
        a0[ic,1] <- achol[i,j]
        ic <- ic+1
    }
}

ssig1 <- matrix(0,nrow=p,ncol=1)

for (i in 1:p){
    ssig1[i,1] <-log(ssig[i,i]^2)

}

hbar1 <- solve(tau*hbar)

hdraw <- matrix(0,nrow=p,ncol=p)
a02mo <- matrix(0,nrow=numa,ncol=numa)
a0mean <-matrix(0,nrow=numa,ncol=1)

for(irep in 1:4000){
 hdraw <-  wish(hbar1,tau) 
  hdraw <- solve(hdraw)
achol <- t(chol(hdraw))
ssig <- matrix(0,nrow=p,ncol=p)
for (i in 1:p){
    ssig[i,i] <- achol[i,i]
    for (j in 1:p){
        achol[j,i] <- achol[j,i]/ssig[i,i]
    }
}
achol <- solve(achol)
a0draw <- matrix(0,nrow=numa,ncol=1)
ic=1
for (i in 2:p){
    for (j in 1:(i-1))
        {
        a0draw[ic,1] <- achol[i,j]
        ic <- ic+1
        
    }
}
a02mo <- a02mo + a0draw%*%t(a0draw)
    a0mean <- a0mean + a0draw
}

a02mo  <-  a02mo/4000
a0mean <- a0mean/4000;
a02mo <- a02mo - a0mean%*%t(a0mean)
results <- list(   B_OLS = aols,VB_OLS =vbar,A_OLS=a0,sigma_OLS =ssig1,VA_OLS =a02mo)
return(results)
}
#############################################################
# Carter and Kohn (1994), On Gibbs sampling for state space models.
carter_kohn <- function(y,Z,Ht,Qt,m,p,t,B0,V0){

# Kalman Filter
bp <- B0

Vp <- V0
bt <- matrix(0,nrow=t, ncol=m)
Vt <- matrix(0, nrow=m^2,ncol=t)
log_lik <- 0
R <- matrix(0,nrow=p,ncol=p)
H <- matrix(0,nrow=p,ncol=m)

for (i in 1:t){

R[1:p,] <- Ht[((i-1)*p+1):(i*p),]
H[1:p,] <- Z[((i-1)*p+1):(i*p),]
cfe <- y[,i] - H%*%bp   #conditional forecast error
f<- H%*%Vp%*%t(H)+R #variance of the conditional forecast error
inv_f<- solve(f)
#log_lik<-log(det(f))+t(cfe)%*%inv_f%*%cfe
btt <- bp+Vp%*%t(H)%*%inv_f%*%cfe
 Vtt <- Vp -Vp%*%t(H)%*%inv_f%*%H%*%Vp
if(i<t){
    bp <- btt
    Vp <- Vtt+Qt
    
}
bt[i,] <- btt
Vt[,i] <- c(Vtt)
}

# draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw <-matrix(0,nrow=t,ncol=m)

bdraw[t,] <- mvrnorm(mu=btt,Sigma=Vtt,n=1)

# Backward recurssions
for (i in 1:(t-1)){
    bf <- (bdraw[(t-i+1),])
    btt <- (bt[(t-i),])
    Vtt <- matrix(Vt[,t-i],nrow=m)
    f <- Vtt + Qt
    inv_f <- solve(f)
    cfe <- bf - btt
    bmean <- btt + Vtt%*%inv_f%*%cfe
    bvar <- Vtt - Vtt%*%inv_f%*%Vtt
    bdraw[t-i,] = mvrnorm(mu=bmean,Sigma=bvar,n=1); #bmean' + randn(1,m)*chol(bvar);
}
bdraw <- t(bdraw)
    return(bdraw)
    }
#############################################################
# Carter and Kohn (1994), On Gibbs sampling for state space models.
carter_kohn2 <- function(y,Z,Ht,Qt,m,p,t,B0,V0,kdraw){

# Kalman Filter

bp <- B0

Vp <- V0
bt <- matrix(0,nrow=t, ncol=m)
Vt <- matrix(0, nrow=m^2,ncol=t)
log_lik <- 0

# R <- matrix(0,nrow=p,ncol=p)
 H <- matrix(0,nrow=p,ncol=m)

for (i in 1:t){

# R[1:p,] <- Ht[((i-1)*p+1):(i*p),]
H[1:p,] <- Z[((i-1)*p+1):(i*p),]
# F <- diag(m)    
cfe <- y[,i] - H%*%bp   #conditional forecast error
f<- H%*%Vp%*%t(H)+Ht[,,i] #variance of the conditional forecast error
# inv_f<- solve(f)
inv_f<- t(H)/f    
#log_lik<-log(det(f))+t(cfe)%*%inv_f%*%cfe
btt <- bp+Vp%*%inv_f%*%cfe
 Vtt <- Vp -Vp%*%inv_f%*%H%*%Vp
if(i<t){
    bp <- btt
    Vp <- Vtt+kdraw[i,]%*%Qt
    
}
bt[i,] <- t(btt)# check the dim
Vt[,i] <- c(Vtt)
}

# draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw <-matrix(0,nrow=t,ncol=m)

bdraw[t,] <- mvrnorm(mu=btt,Sigma=Vtt,n=1)

# Backward recurssions
for (i in 1:(t-1)){
    bf <- (bdraw[(t-i+1),])
    btt <- (bt[(t-i),])
    Vtt <- matrix(Vt[,t-i],nrow=m)
    f <-  Vtt+kdraw[t-i,]%*%Qt
#     inv_f <- solve(f)
    inv_f <- Vtt/f
    cfe <- bf - btt
    bmean <- btt + inv_f%*%cfe
    bvar <- Vtt - inv_f%*%Vtt
    bdraw[t-i,] = mvrnorm(mu=bmean,Sigma=bvar,n=1); #bmean' + randn(1,m)*chol(bvar);
}
bdraw <- t(bdraw)
    return(bdraw)
   }

## function Repeats columns, and rows in a matrix
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

SVRW2 <- function(Ystar,h,sig,sigma_prmean,sigma_prvar,TVP_Sigma){
    

T <- length(h)

# normal mixture
pi  <- c(0.0073, .10556, .00002, .04395, .34001, .24566, .2575)
mi <-c(-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518, -1.08819) - 1.2704 
# means already adjusted!! %%
sigi <- c(5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261)
sqrtsigi <- sqrt(sigi)


# Sample S from a 7-point distrete distribution
temprand <- as.matrix(runif(T))




q <-rep.row(pi,T)*dnorm(rep.col(Ystar,7),(rep.col(h,7)+rep.row(mi,T)),rep.row(sqrtsigi,T)) 
# element wise multiplication
q <- q/rep.col(rowSums(q),7)

S <- 7-as.matrix(rowSums(1*(rep.col(temprand,7)<t(apply(q,1,cumsum)))))+1

vart <- array(0,c(1,1,T))
yss1 <- matrix(0,nrow=T,ncol=1)

for (i in 1:T){
    imix <-S[i]
    vart[1,1,i] <-sigi[imix]
    yss1[i,1] <- Ystar[i,1]-mi[imix]
}

h <- carter_kohn2(t(yss1),matrix(1,nrow=T,ncol=1),
                  vart,sig,m=1,p=1,t=T,sigma_prmean,sigma_prvar,kdraw =(TVP_Sigma*matrix(1,nrow=T,ncol=1)))


return(h)
}
#########################################
#Computes a correlation matrix from a
#               variance-covariance matrix.
corrvc <- function(vc){
if (is.complex(vc) == TRUE){       
    print("ERROR: Not implemented for complex arguments.")
} else {
std <- sqrt(as.matrix(diag(vc)))    
}

     output = vc/(std%*%t(std))

    return(output)
}
