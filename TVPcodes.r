#######################################################################################
Note: The R codes are written based on MATLAB codes prepared by Garry Koop and vailable 
at https://sites.google.com/site/garykoop/home/computer-code-2.
########################################################################################


rm(list=ls())

## Load the packaes required for the programme
#library("largeList")
library("MASS") # for mvrnorm
source("TVP_funs.R")

#################Load the data##################################

Y <- read.csv("data.csv",header=T)
Y <- as.matrix(Y[,c(-1,-4)]) # For Model1- with retruns
#Y <- as.matrix(Y[,c(-1,-5)]) # For Model2- with volatility.

#Number of observations and dimension of X and Y
t <- length(Y[,1])
M=ncol(Y)

#Number of factors & lags:
tau <- 40 # tau is the size of the training sample
p <- 2 # p is number of lags in the VAR part
numa <- M*(M-1)/2 # Number of lower triangular elements of A_t (other than 0's and 1's)

#Generate lagged Y matrix. This will be part of the X matrix
ylag1 <- makelag(Y,1)
ylag2 <- makelag(Y,2)

ylag  <- cbind(ylag1,ylag2) # Y is [T x M]. ylag is [T x (Mp)]
# Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
#ylag = ylag(p+tau+1:t,:)
ylag <- ylag[-(1:(tau+p)),]

K  <- M + p*(M^2) # K is the number of elements in the state vector

Z <- matrix(0,ncol=K,nrow=(t-tau-p)*M)#Create Z_t matrix

for (i in 1:(t-tau-p))
{
    ztemp <- diag(1,M)
    for (j in 1:p)
        {       
        xtemp = ylag[i,(((j-1)*M+1):(j*M))]
        xtemp = t(kronecker(diag(1,M),xtemp))
        ztemp = cbind(ztemp, xtemp)  
        }
    Z[(((i-1)*M+1):(i*M)),] <- ztemp
}


y <- t(Y[((tau+p+1):t),])
t <- ncol(y)
####################### Preliminaries#############################################
#Set some Gibbs - related preliminaries
nrep = 30000  # Number of replications
nburn = 10000 # Number of burn-in-draws
#it_print = 100;  %Print in the screen every "it_print"-th iteration
bn <- 3000


####PRIORS########
# the ts_prior() function is given seperately
priors <- ts_prior(Y,tau,M,p)
A_OLS <- priors$A_OLS
B_OLS <- priors$B_OLS
VA_OLS <- priors$VA_OLS
VB_OLS <- priors$VB_OLS
sigma_OLS <- priors$sigma_OLS
sigma_prmean <- sigma_OLS
sigma_prvar <- 4*diag(M)
S_prmean <- mm<-matrix(list(), M-1, 1)
S_prvar <- mm<-matrix(0, M-1, 1)

#######INITIALIZE MATRICES########
consQ <- 0.0001
consS <- 0.0001
consH <- 0.01
consW <- 0.0001


#  Set some hyperparameters 
k_Q <- 0.01
k_S <- 0.1
k_W <- 0.01

sizeW <- M # Size of matrix W
sizeS <- 1:M #Size of matrix S

B_0_prmean <- B_OLS
B_0_prvar <-  4*VB_OLS

# A_0 ~ N(A_OLS, 4Var(A_OLS))
A_0_prmean <- A_OLS
A_0_prvar <- 4*VA_OLS


# Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean <- ((k_Q)^2)*tau*VB_OLS
Q_prvar <- tau

# W ~ IG(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean <- ((k_W)^2)*matrix(1,nrow=M,ncol=1)
W_prvar = 2

ind <-1
for (ii in 2:M){
    #S is block diagonal 
S_prmean[[ii-1]] <-((k_S)^2)*(1 + sizeS[(ii-1)])*VA_OLS[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+
                                                            (ii-3)*(ii-2)/2):ind]
S_prvar[ii-1] = 1 + sizeS[ii-1]
    ind=ind+ii}

Ht <- kronecker(matrix(1,nrow=t,ncol=1),consH*diag(M)) 


Htchol <- kronecker(matrix(1,nrow=t,ncol=1),sqrt(consH)*diag(M))# Cholesky of Htdraw defined above

Qdraw <- consQ*diag(K) # Initialize Qdraw, a draw from the covariance matrix Q
Sdraw <- consS*diag(numa)# Initialize Sdraw, a draw from the covariance matrix S

Sblockdraw <- matrix(list(), M-1, 1)
ijc = 1
for (jj in 2:M){
    Sblockdraw[[jj-1]] <-Sdraw[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc]
    ijc <- ijc + jj
}

Wdraw <- consW*matrix(1,nrow=M,ncol=1) # Initialize Wdraw, a draw from the covariance matrix W
Btdraw <-matrix(0,nrow=K,ncol=t) # Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
Atdraw <-matrix(0,nrow=numa,ncol=t) # Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)
Sigtdraw <-matrix(0,nrow=t,ncol=M) # Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
# Matrix of the exponent of Sigtdraws (SIGMA(t))
sigt <- kronecker(matrix(0,nrow=t, ncol=1),0.01*diag(M)) 
statedraw <- 5*matrix(1,nrow=t,ncol=M) # initialize the draw of the indicator variable 
                               # (of 7-component mixture of Normals approximation)
Zs <-kronecker(matrix(1,nrow=t,ncol=1),diag(M))

#  Storage matrices for posteriors and stuff
Bt_postmean <-matrix(0,nrow=K,ncol=t)#  regression coefficients B(t)
At_postmean <-matrix(0,nrow=numa,ncol=t)#  lower triangular matrix A(t)
Sigt_postmean <-matrix(0,nrow=t,ncol=M) #  diagonal std matrix SIGMA(t)
Qmean <-matrix(0,nrow=K,ncol=K) #  covariance matrix Q of B(t)
Smean <-matrix(0,nrow=numa,ncol=numa)#  covariance matrix S of A(t)
Wmean <-matrix(0,nrow=M,ncol=1)#  covariance matrix W of SIGMA(t)

sigmean <-matrix(0,nrow=t,ncol=M)#  mean of the diagonal of the VAR covariance matrix
cormean <-matrix(0,nrow=t,ncol=numa) #  mean of the off-diagonal elements of the VAR cov matrix
sig2mo <-matrix(0,nrow=t,ncol=M)# squares of the diagonal of the VAR covariance matrix
cor2mo <-matrix(0,nrow=t,ncol=numa)#  squares of the off-diagonal elements of the VAR cov matrix


###### IMPULSE RESPONSES Storage matrices######

istore <-1
if (istore == 1) { 
    nhor <- 21; # Impulse response horizon
    bigj <- matrix(0, nrow=M, ncol=M*p)
    bigj[1:M,1:M] = diag(M)
    }
imprespt <- matrix(0,nrow=M*t,ncol=M*nhor)
#impres01 <- array(0, dim = c(M*t,M*nhor,batchsz))

################################# END OF PRELIMINARIES #################################
####################################START SAMPLING#####################################
# Progress Bar

total <- nrep + nburn
pb <- txtProgressBar(min = 0, max = total, style = 3)
for (irep  in 1:(nrep + nburn)){   
Sys.sleep(0.1)
################  STEP I: Sample B from p(B|y,A,Sigma,V)#################


# Accept draw
 Btdraw <- carter_kohn(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar)


#=====| Draw Q, the covariance of B(t) (from iWishart)
# Take the SSE in the state equation of B(t)

Btemp <- t(Btdraw[,2:t]) - t(Btdraw[,1:t-1])

sse_2 <- matrix(0,nrow=K,ncol=K)
for (i in 1:t-1){
    sse_2 <- sse_2+as.matrix(Btemp[1,])%*%t(Btemp[1,])
}

#  ...and subsequently draw Q, the covariance matrix of B(t)
Qinv <- solve(sse_2+Q_prmean,tol=1e-21)

Qinvdraw <- wish(Qinv,t-1+Q_prvar)
Qdraw <- solve(Qinvdraw) # this is a draw from Q

######################  STEP II: Draw A(t) from p(At|y,B,Sigma,V) ######################

# Drwas alpha

yhat <- matrix(0, nrow=M,ncol=t)
for (i in 1: t){
    yhat[,i] <- y[,i]-Z[(((i-1)*M+1):(i*M)),]%*%Btdraw[,i]
}

Zc <- -t(yhat)
sigma2temp <- exp(Sigtdraw)

 # Draw each block of A(t)
Atdraw <- carter_kohn(t(as.matrix(yhat[2,])),as.matrix(Zc[,1]),as.matrix(sigma2temp[,2]),
            Sblockdraw[[1]],sizeS[1],1,t,A_0_prmean[1,],A_0_prvar[1,1])
       

ind <- 3
for (ii in 3:M){
Atblockdraw <- carter_kohn(t(as.matrix(yhat[ii,])),Zc[,1:(ii-1)],as.matrix(sigma2temp[,ii]),
            Sblockdraw[[ii-1]],sizeS[ii-1],1,t,A_0_prmean[(((ii-1)+(ii-3)*(ii-2)/2):ind),],
                A_0_prvar[(((ii-1)+(ii-3)*(ii-2)/2):ind),(((ii-1)+(ii-3)*(ii-2)/2):ind)]);
        Atdraw <- rbind(Atdraw,Atblockdraw) # Atdraw is the final matrix of draws of A(t)
        ind = ind + ii;
    }



####### Draw S, the covariance of A(t) (from iWishart)
    #### Take the SSE in the state equation of A(t)

Attemp <- t(Atdraw[,2:t]) - t(Atdraw[,(1:(t-1))])
sse_2 <- matrix(0,nrow=numa,ncol=numa)

for (i in 1:(t-1)){
sse_2+(Attemp[i,])%*%t(Attemp[i,])
}

 
ijc <-1

for(jj in 2:M){
Sinv <- solve(sse_2[(((jj-1)+(jj-3)*(jj-2)/2):ijc),(((jj-1)+(jj-3)*(jj-2)/2):ijc)] + 
                   S_prmean[[jj-1]])

Sinvblockdraw <- wish(h=Sinv,n=t-1+S_prvar[jj-1])

Sblockdraw[[jj-1]] = solve(Sinvblockdraw)
ijc <- ijc+jj
}

########################### STEP III: Draw diagonal VAR covariance matrix log-SIGMA(t)#############

  

# First create capAt, the lower-triangular matrix A(t) with ones on the
# main diagonal. 

capAt <- diag(M)
aatemp <- as.matrix(Atdraw[,1])
    ic=1
    for (j in 2:M){
        capAt[j,(1:(j-1))] <- t(aatemp[(ic :(ic+j-2)),1])
        ic <- ic+j-1
    }

for (i in 2:t){
    capatemp <- diag(M)
    aatemp <- as.matrix(Atdraw[,i])
     ic=1
    for (j in 2:M){
        capatemp[j,(1:(j-1))] <- t(aatemp[(ic :(ic+j-2)),1])
        ic <- ic+j-1
    }
    capAt<-rbind(capAt,capatemp)
}

ytemps = capAt[(1:M),]%*%yhat[,1]
y2 <- (ytemps^2)


for (i in 2:t){
        ytemps <-  capAt[((i-1)*M+1):(i*M),]%*%yhat[,i]
        y2 <- cbind(y2, (ytemps^2))
  }

yss <-t(log(y2+ 10^(-6)))

for (j in 1:M){
Sigtdraw[,j] <- SVRW2(as.matrix(yss[,j]),Sigtdraw[,j],Wdraw[j,],sigma_prmean[j],
                      sigma_prvar[j,j],1)  
    }

sigt <- exp(0.5*Sigtdraw)
e2 <- Sigtdraw[-1,] - Sigtdraw[-(nrow(Sigtdraw)),]
W1 <- W_prvar + t - p - 1
W2 <- W_prmean + rowSums(t(e2^2))
Winvdraw <- matrix(rgamma(n=M,shape=W1/2,scale=2/W2))
Wdraw <- 1/Winvdraw

#  Create the VAR covariance matrix H(t). It holds that:
#  A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '

Ht <- matrix(0,nrow=M*t,M)
Htsd <- matrix(0,nrow=M*t,M)

for (i in 1:t){
    inva <- solve(capAt[(((i-1)*M+1):(i*M)),])
    stem <- diag(sigt[i,])
    Hsd <- inva%*%stem
    Hdraw <- Hsd%*%t(Hsd)
    Ht[(((i-1)*M+1):(i*M)),] <- Hdraw #H(t)
    Htsd[(((i-1)*M+1):(i*M)),] <- Hsd  #Cholesky of H(t)
}

##########SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES ##############
if (irep > nburn){               
Bt_postmean <- Bt_postmean + Btdraw # regression coefficients B(t)
At_postmean <- At_postmean + Atdraw # lower triangular matrix A(t)
Sigt_postmean <- Sigt_postmean + Sigtdraw  # diagonal std matrix SIGMA(t)
Qmean <- Qmean + Qdraw# covariance matrix Q of B(t)

ikc <-1
for(kk in 2:M){
    Sdraw[(((kk-1)+(kk-3)*(kk-2)/2):ikc),(((kk-1)+(kk-3)*(kk-2)/2):ikc)] <-Sblockdraw[[kk-1]]
    ikc <- ikc + kk
}
            

Smean <-  Smean + Sdraw #covariance matrix S of A(t)
Wmean <- Wmean + Wdraw #covariance matrix W of SIGMA(t)
#  Get time-varying correlations and variances

# for the first period
stemp6 <- matrix(0,nrow=M,ncol=1)
stemp8 <- corrvc(Ht[1:M,])
stemp7a <- matrix(0,nrow=M,ncol=0)
ic <-1
for(j in 1:M){
    if (j>1){
       stemp7a <- c(stemp7a,t(stemp8[j,(1:ic)])) 
       ic <- ic+1             
    } 
    stemp6[j,1] <- sqrt(Ht[j,j])
     
    }

stemp5 <- t(stemp6)
stemp7 <- t(stemp7a)                    
 

                  
for (i in 2:t) {
stemp6 <- matrix(0,nrow=M,ncol=1)
stemp8 <- corrvc(Ht[((i-1)*M+1):(i*M),])
stemp7a <- matrix(0,nrow=M,ncol=0)
ic <-1
for(j in 1:M){
    if (j>1){
       stemp7a <- c(stemp7a,t(stemp8[j,(1:ic)])) 
       ic <- ic+1             
    } 
    stemp6[j,1] <- sqrt(Ht[((i-1)*M+j),j])
     
    }
 stemp5 <- rbind(stemp5,t(stemp6))
 stemp7 <- rbind(stemp7,t(stemp7a))  
}                   


sigmean <- sigmean + stemp5 # diagonal of the VAR covariance matrix
cormean <- cormean + stemp7 # off-diagonal elements of the VAR cov matrix
sig2mo <- sig2mo + stemp5^2
cor2mo <- cor2mo + stemp7^2

# Impulse response analysis. 
if (istore==1){
            
           
biga <- matrix(0,nrow=M*p,ncol=M*p) # M= number of variables, P is nu

for(j in 1:(p-1)){
    biga[(j*M+1):(M*(j+1)), (M*(j-1)+1):(j*M)] <- diag(M)
}

for (i in 1:t){ #Get impulses recurssively for each time period
     bbtemp <- Btdraw[((M+1):K),i]
     splace <- 0
     for (ii in 1:p){
         for (iii in 1:M){
             biga[iii,(((ii-1)*M+1):(ii*M))] <- bbtemp[((splace+1):(splace+M))]
            splace <- splace + M 
         }
     }
#Identification code: 
# St dev matrix for structural VAR
Hsd <- Htsd[(((i-1)*M+1):(i*M)),(1:M)]  
Hsd
#  First shock is the Cholesky of the VAR covariance
diagonal = diag(diag(Hsd))
Hsd = solve(diagonal)%*%Hsd   # Unit initial shock

# Now get impulse responses for 1 through nhor future periods
impresp <- matrix(nrow=M,ncol=M*nhor)
impresp[1:M,1:M] <-Hsd

# First shock is the Cholesky of the VAR covariance
bigai <-biga

for (j in 1:(nhor-1)){
    impresp[,(j*M+1):((j+1)*M)] <- bigj%*%bigai%*%t(bigj)%*%Hsd
    bigai <- bigai%*%biga
}
if(i==1){ # please check
imprespt[1:M,] <-impresp
    } else{
    imprespt[((i-1)*M+1):(i*M),] <-impresp
}
 
#### the draws are saved as different files using large list package #####   
}
if((irep-nburn)<=bn){    
if((irep-nburn) == 1){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc1.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc1.llo", append = TRUE)
    }
} else if ((irep-nburn)>bn & (irep-nburn)<=(2*bn)){
    if((irep-nburn) == (bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc2.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc2.llo", append = TRUE)
    }
} else if ((irep-nburn)>(2*bn) & (irep-nburn)<=(3*bn)){
   if((irep-nburn) == (2*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc3.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc3.llo", append = TRUE)
    }  
} else if ((irep-nburn)>(3*bn) & (irep-nburn)<=(4*bn)){
   if((irep-nburn) == (3*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc4.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc4.llo", append = TRUE)
    }  
} else if ((irep-nburn)>(4*bn) & (irep-nburn)<=(5*bn)){
   if((irep-nburn) == (4*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc5.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc5.llo", append = TRUE)
    }  
} else if ((irep-nburn)>(5*bn) & (irep-nburn)<=(6*bn)){
   if((irep-nburn) == (5*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc6.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc6.llo", append = TRUE)
    }  
}else if ((irep-nburn)>(6*bn) & (irep-nburn)<=(7*bn)){
   if((irep-nburn) == (6*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc7.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc7.llo", append = TRUE)
    }  
}else if ((irep-nburn)>(7*bn) & (irep-nburn)<=(8*bn)){
   if((irep-nburn) == (7*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc8.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc8.llo", append = TRUE)
    }  
}else if ((irep-nburn)>(8*bn) & (irep-nburn)<=(9*bn)){
   if((irep-nburn) == (8*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc9.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc9.llo", append = TRUE)
    }  
}else {
   if((irep-nburn) == (9*bn+1)){  
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc10.llo", append = FALSE, compress = TRUE)
    } else { 
    name1 <- paste("S",irep-nburn,sep="")
    ar12 <- list(A=imprespt)
    names(ar12)[1] <- name1[1]
    saveList(object = ar12, file = "irfmcmc10.llo", append = TRUE)
    }  
}
    
    
} #END the impulse response calculation section 

} #END saving after burn-in results 

# update progress bar
   setTxtProgressBar(pb, irep)
} #END main Gibbs loop (for irep = 1:nrep+nburn)

close(pb) 
###############################End of Sampling###########################


################ Codes to create quantiles of IRF and to save it as excel file##### 
# The codes to arrange the simple IRF draws saved in .llo files in the previous section

impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc1.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd1.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc1.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd1.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc2.llo", index = i)[[1]]
impres00[i,] <- hhh[1,] }  
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd2.llo", append = FALSE, compress = TRUE)    

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc2.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd2.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc3.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]   }
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd3.llo", append = FALSE, compress = TRUE)    

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.01)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc3.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd3.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)



impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc4.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]   }
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd4.llo", append = FALSE, compress = TRUE)    

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc4.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd4.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)


impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc5.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd5.llo", append = FALSE, compress = TRUE)    

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc5.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd5.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)


impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc6.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]}   
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd6.llo", append = FALSE, compress = TRUE)    

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc6.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd6.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)


impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc7.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]   }
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd7.llo", append = FALSE, compress = TRUE)    

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc7.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd7.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)


impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc8.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]   }
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd8.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc8.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd8.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)


impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc9.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]   }
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd9.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc9.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd9.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)


impres00 <- matrix(0, nrow=bn, ncol=M*nhor)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc10.llo", index = i)[[1]]
impres00[i,] <- hhh[1,]   }
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd10.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hhh <- readList(file = "irfmcmc10.llo", index = i)[[1]]
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd10.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

rm(impres00,ar12)

imp50 <- matrix(0,nrow= t*M,ncol=M*nhor)
imp16 <- matrix(0,nrow= t*M,ncol=M*nhor)
imp84 <- matrix(0,nrow= t*M,ncol=M*nhor)

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (i in 1: (t*M)){
Sys.sleep(0.1)
imp001 <- readList("irfarrd1.llo",index=i)[[1]]
imp001 <- rbind(imp001,readList("irfarrd2.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd3.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd4.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd5.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd6.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd7.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd8.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd9.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd10.llo",index=i)[[1]])

# calculates 16th, median and 84th quantile IRFs
imp50[i,] <- apply(imp001, 2, quantile, probs=0.5)
imp16[i,] <- apply(imp001, 2, quantile, probs=0.16)
imp84[i,] <- apply(imp001, 2, quantile, probs=0.84)
setTxtProgressBar(pb, i)
}
close(pb)
rm(imp001)

imp16_n1 <- seq(1, ncol(imp16), by=4) # first equation
imp16_1 <- imp16[,imp16_n1]

# first to first (response of op to op)
imp16_n11 <- seq(1, nrow(imp16_1), by=4)
imp16_11 <- imp16_1[imp16_n11,]#####################

# second to first(response of xoi to op)
imp16_n12 <- seq(2, nrow(imp16_1), by=4)
imp16_12 <- imp16_1[imp16_n12,]#####################

# Third to firtst (response of epi to op)
imp16_n13 <- seq(3, nrow(imp16_1), by=4)
imp16_13 <- imp16_1[imp16_n13,]#####################



# Fourth to firtst (response of sigma to op)
imp16_n14 <- seq(4, nrow(imp16_1), by=4)
imp16_14 <- imp16_1[imp16_n14,]#####################

rm(imp16_n11,imp16_n12,imp16_n13,imp16_1,imp16_n1,imp16_n14)

imp16_n2 <- seq(2, ncol(imp16), by=4) # second equation
imp16_2 <- imp16[,imp16_n2]

# second to first (Response of op to xoi)
imp16_n21 <- seq(1, nrow(imp16_2), by=4)
imp16_21 <- imp16_2[imp16_n21,]#####################

# second to second (Response of xoi to xoi)
imp16_n22 <- seq(2, nrow(imp16_2), by=4)
imp16_22 <- imp16_2[imp16_n22,]#####################

# second to third(Response of epi to xoi)
imp16_n23 <- seq(3, nrow(imp16_2), by=4)
imp16_23 <- imp16_2[imp16_n23,]#####################

# second to (Response of sig to xoi)
imp16_n24 <- seq(4, nrow(imp16_2), by=4)
imp16_24 <- imp16_2[imp16_n24,]#####################

rm(imp16_2,imp16_n21,imp16_n22,imp16_n23,imp16_n2,imp16_n24)

imp16_n3 <- seq(3, ncol(imp16), by=4) # 3rd equation
imp16_3 <- imp16[,imp16_n3]

# third to first (Response of op to epi)
imp16_n31 <- seq(1, nrow(imp16_3), by=4)
imp16_31 <- imp16_3[imp16_n31,]#####################

# third to second  (Response of xoi to epi)
imp16_n32 <- seq(2, nrow(imp16_3), by=4)
imp16_32 <- imp16_3[imp16_n32,]#####################


# third to third  (Response of epi to epi)
imp16_n33 <- seq(3, nrow(imp16_3), by=4)
imp16_33 <- imp16_3[imp16_n33,]#####################

# third to third  (Response of sig to epi)
imp16_n34 <- seq(4, nrow(imp16_3), by=4)
imp16_34 <- imp16_3[imp16_n34,]#####################
rm(imp16_3,imp16_n31,imp16_n32,imp16_n33,imp16_n3,imp16_n34)

imp16_n4 <- seq(4, ncol(imp16), by=4) # 3rd equation
imp16_4 <- imp16[,imp16_n4]

# fort to first (Response of op to sig)
imp16_n41 <- seq(1, nrow(imp16_4), by=4)
imp16_41 <- imp16_4[imp16_n41,]#####################

# third to second  (Response of xoi to sig)
imp16_n42 <- seq(2, nrow(imp16_4), by=4)
imp16_42 <- imp16_4[imp16_n42,]#####################

# third to third  (Response of epi to sig)
imp16_n43 <- seq(3, nrow(imp16_4), by=4)
imp16_43 <- imp16_4[imp16_n43,]#####################

# third to third  (Response of sig to sig)
imp16_n44 <- seq(4, nrow(imp16_4), by=4)
imp16_44 <- imp16_4[imp16_n44,]#####################

rm(imp16_4,imp16_n41,imp16_n42,imp16_n43,imp16_n4,imp16_n44)


imp50_n1 <- seq(1, ncol(imp50), by=4) # first equation
imp50_1 <- imp50[,imp50_n1]

# first to first (response of op to op)
imp50_n11 <- seq(1, nrow(imp50_1), by=4)
imp50_11 <- imp50_1[imp50_n11,]#####################

# second to first(response of xoi to op)
imp50_n12 <- seq(2, nrow(imp50_1), by=4)
imp50_12 <- imp50_1[imp50_n12,]#####################

# Third to firtst (response of epi to op)
imp50_n13 <- seq(3, nrow(imp50_1), by=4)
imp50_13 <- imp50_1[imp50_n13,]#####################

# Fourth to firtst (response of sigma to op)
imp50_n14 <- seq(4, nrow(imp50_1), by=4)
imp50_14 <- imp50_1[imp50_n14,]#####################

rm(imp50_n11,imp50_n12,imp50_n13,imp50_1,imp50_n1,imp50_n14)

imp50_n2 <- seq(2, ncol(imp50), by=4) # second equation
imp50_2 <- imp50[,imp50_n2]

# second to first (Response of op to xoi)
imp50_n21 <- seq(1, nrow(imp50_2), by=4)
imp50_21 <- imp50_2[imp50_n21,]#####################

# second to second (Response of xoi to xoi)
imp50_n22 <- seq(2, nrow(imp50_2), by=4)
imp50_22 <- imp50_2[imp50_n22,]#####################

# second to third(Response of epi to xoi)
imp50_n23 <- seq(3, nrow(imp50_2), by=4)
imp50_23 <- imp50_2[imp50_n23,]#####################

# second to (Response of sig to xoi)
imp50_n24 <- seq(4, nrow(imp50_2), by=4)
imp50_24 <- imp50_2[imp50_n24,]#####################

rm(imp50_2,imp50_n21,imp50_n22,imp50_n23,imp50_n2,imp50_n24)

imp50_n3 <- seq(3, ncol(imp50), by=4) # 3rd equation
imp50_3 <- imp50[,imp50_n3]

# third to first (Response of op to epi)
imp50_n31 <- seq(1, nrow(imp50_3), by=4)
imp50_31 <- imp50_3[imp50_n31,]#####################

# third to second  (Response of xoi to epi)
imp50_n32 <- seq(2, nrow(imp50_3), by=4)
imp50_32 <- imp50_3[imp50_n32,]#####################

# third to third  (Response of epi to epi)
imp50_n33 <- seq(3, nrow(imp50_3), by=4)
imp50_33 <- imp50_3[imp50_n33,]#####################

# third to third  (Response of sig to epi)
imp50_n34 <- seq(4, nrow(imp50_3), by=4)
imp50_34 <- imp50_3[imp50_n34,]#####################

rm(imp50_n31,imp50_n32,imp50_n33,imp50_n3,imp50_n34)

imp50_n4 <- seq(4, ncol(imp50), by=4) # 3rd equation
imp50_4 <- imp50[,imp50_n4]

# fort to first (Response of op to sig)
imp50_n41 <- seq(1, nrow(imp50_4), by=4)
imp50_41 <- imp50_4[imp50_n41,]#####################

# third to second  (Response of xoi to sig)
imp50_n42 <- seq(2, nrow(imp50_4), by=4)
imp50_42 <- imp50_4[imp50_n42,]#####################

# third to third  (Response of epi to sig)
imp50_n43 <- seq(3, nrow(imp50_4), by=4)
imp50_43 <- imp50_4[imp50_n43,]#####################

# third to third  (Response of sig to sig)
imp50_n44 <- seq(4, nrow(imp50_3), by=4)
imp50_44 <- imp50_4[imp50_n44,]#####################

rm(imp50_4,imp50_n41,imp50_n42,imp50_n43,imp50_n4,imp50_n44)


imp84_n1 <- seq(1, ncol(imp84), by=4) # first equation
imp84_1 <- imp84[,imp84_n1]

# first to first (response of op to op)
imp84_n11 <- seq(1, nrow(imp84_1), by=4)
imp84_11 <- imp84_1[imp84_n11,]#####################

# second to first(response of xoi to op)
imp84_n12 <- seq(2, nrow(imp84_1), by=4)
imp84_12 <- imp84_1[imp84_n12,]#####################

# Third to firtst (response of epi to op)
imp84_n13 <- seq(3, nrow(imp84_1), by=4)
imp84_13 <- imp84_1[imp84_n13,]#####################

# Fourth to firtst (response of sigma to op)
imp84_n14 <- seq(4, nrow(imp84_1), by=4)
imp84_14 <- imp84_1[imp84_n14,]#####################

rm(imp84_n11,imp84_n12,imp84_n13,imp84_1,imp84_n1,imp84_n14)

imp84_n2 <- seq(2, ncol(imp84), by=4) # second equation
imp84_2 <- imp84[,imp84_n2]

# second to first (Response of op to xoi)
imp84_n21 <- seq(1, nrow(imp84_2), by=4)
imp84_21 <- imp84_2[imp84_n21,]#####################

# second to second (Response of xoi to xoi)
imp84_n22 <- seq(2, nrow(imp84_2), by=4)
imp84_22 <- imp84_2[imp84_n22,]#####################

# second to third(Response of epi to xoi)
imp84_n23 <- seq(3, nrow(imp84_2), by=4)
imp84_23 <- imp84_2[imp84_n23,]#####################

# second to (Response of sig to xoi)
imp84_n24 <- seq(4, nrow(imp84_2), by=4)
imp84_24 <- imp84_2[imp84_n24,]#####################

rm(imp84_2,imp84_n21,imp84_n22,imp84_n23,imp84_n2,imp84_n24)

imp84_n3 <- seq(3, ncol(imp84), by=4) # 3rd equation
imp84_3 <- imp84[,imp84_n3]

# third to first (Response of op to epi)
imp84_n31 <- seq(1, nrow(imp84_3), by=4)
imp84_31 <- imp84_3[imp84_n31,]#####################

# third to second  (Response of xoi to epi)
imp84_n32 <- seq(2, nrow(imp84_3), by=4)
imp84_32 <- imp84_3[imp84_n32,]#####################

# third to third  (Response of epi to epi)
imp84_n33 <- seq(3, nrow(imp84_3), by=4)
imp84_33 <- imp84_3[imp84_n33,]#####################

# third to third  (Response of sig to epi)
imp84_n34 <- seq(4, nrow(imp84_3), by=4)
imp84_34 <- imp84_3[imp84_n34,]#####################

rm(imp84_3,imp84_n31,imp84_n32,imp84_n33,imp84_n3,imp84_n34)

imp84_n4 <- seq(4, ncol(imp84), by=4) # 3rd equation
imp84_4 <- imp84[,imp84_n4]

# fort to first (Response of op to sig)
imp84_n41 <- seq(1, nrow(imp84_4), by=4)
imp84_41 <- imp84_4[imp84_n41,]#####################

# third to second  (Response of xoi to sig)
imp84_n42 <- seq(2, nrow(imp84_4), by=4)
imp84_42 <- imp84_4[imp84_n42,]#####################

# third to third  (Response of epi to sig)
imp84_n43 <- seq(3, nrow(imp84_4), by=4)
imp84_43 <- imp84_4[imp84_n43,]#####################

# third to third  (Response of sig to sig)
imp84_n44 <- seq(4, nrow(imp84_4), by=4)
imp84_44 <- imp84_4[imp84_n44,]#####################



rm(imp84_4,imp84_n41,imp84_n42,imp84_n43,imp84_n4,imp84_n44)

#Saves the 16th qunatile simple IRF as csv files###






write.csv(imp16_11 , file = " imp16_11.csv")
write.csv(imp16_12 , file = " imp16_12.csv")
write.csv(imp16_13 , file =" imp16_13.csv")
write.csv(imp16_14 , file = " imp16_14.csv")
write.csv(imp16_21 , file = " imp16_21.csv")
write.csv(imp16_22 , file = " imp16_22.csv")
write.csv(imp16_23 , file = " imp16_23.csv")
write.csv(imp16_24 , file = " imp16_24.csv")
write.csv(imp16_31, file = " imp16_31.csv")
write.csv(imp16_32 , file = " imp16_32.csv")
write.csv(imp16_33 , file = " imp16_33.csv")
write.csv(imp16_34 , file = " imp16_34.csv")
write.csv(imp16_41 , file = " imp16_41.csv")
write.csv(imp16_42 , file = " imp16_42.csv")
write.csv(imp16_43 , file = " imp16_43.csv")
write.csv(imp16_44 , file = " imp16_44.csv")

#Saves the 50 th qunatile simple IRF as csv files###
write.csv(imp50_11 , file = " imp50_11.csv")
write.csv(imp50_12 , file = " imp50_12.csv")
write.csv(imp50_13 , file =" imp50_13.csv")
write.csv(imp50_14 , file = " imp50_14.csv")
write.csv(imp50_21 , file = " imp50_21.csv")
write.csv(imp50_22 , file = " imp50_22.csv")
write.csv(imp50_23 , file = " imp50_23.csv")
write.csv(imp50_24 , file = " imp50_24.csv")
write.csv(imp50_31, file = " imp50_31.csv")
write.csv(imp50_32 , file = " imp50_32.csv")
write.csv(imp50_33 , file = " imp50_33.csv")
write.csv(imp50_34 , file = " imp50_34.csv")
write.csv(imp50_41 , file = " imp50_41.csv")
write.csv(imp50_42 , file = " imp50_42.csv")
write.csv(imp50_43 , file = " imp50_43.csv")
write.csv(imp50_44 , file = " imp50_44.csv")

#Saves the 50 th qunatile simple IRF as csv files###
write.csv(imp84_11 , file = " imp84_11.csv")
write.csv(imp84_12 , file = " imp84_12.csv")
write.csv(imp84_13 , file =" imp84_13.csv")
write.csv(imp84_14 , file = " imp84_14.csv")
write.csv(imp84_21 , file = " imp84_21.csv")
write.csv(imp84_22 , file = " imp84_22.csv")
write.csv(imp84_23 , file = " imp84_23.csv")
write.csv(imp84_24 , file = " imp84_24.csv")
write.csv(imp84_31, file = " imp84_31.csv")
write.csv(imp84_32 , file = " imp84_32.csv")
write.csv(imp84_33 , file = " imp84_33.csv")
write.csv(imp84_34 , file = " imp84_34.csv")
write.csv(imp84_41 , file = " imp84_41.csv")
write.csv(imp84_42 , file = " imp84_42.csv")
write.csv(imp84_43 , file = " imp84_43.csv")
write.csv(imp84_44 , file = " imp84_44.csv")

####################END of codes to create IRF ####################


################ Codes to create quantiles of cummulative IRF and to save it as .csv file##### 



impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc1.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd1c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc1.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd1c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)


impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc2.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd2c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc2.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd2c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc3.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd3c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc3.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd3c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc4.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd4c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc4.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd4c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc5.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd5c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc5.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd5c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc6.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd6c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc6.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd6c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc7.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd7c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc7.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd7c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc8.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd8c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc8.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd8c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc9.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd9c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc9.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd9c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

impres00 <- matrix(0, nrow=bn, ncol=20)
for (i in 1: bn){
hh <- readList(file = "irfmcmc10.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)
    
impres00[i,] <- hhh[1,]  } 
ar12 <- list(r1=impres00)
saveList(object = ar12, file = "irfarrd10c.llo", append = FALSE, compress = TRUE)    


pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (j in 2: (t*M)){
Sys.sleep(0.1)
for (i in 1: bn){
hh <- readList(file = "irfmcmc10.llo", index = i)[[1]]
h1 <- hh[,1:4]
h3 <- hh[,1:4]+hh[,5:8]+hh[,9:12];
h6 <- h3+hh[,13:16]+hh[,17:20]+hh[,21:24]
h9 <-h6+hh[,25:28]+hh[,29:32]+hh[,33:36]
h12 <-h9+hh[,37:40]+hh[,41:44]+hh[,45:48]

hhh <- cbind(h1,h3,h6,h9,h12)    
impres00[i,] <- hhh[j,]
    }
name1 <- paste("r",j,sep="")  
ar12 <- list(A=impres00)
names(ar12)[1] <- name1[1]
saveList(object = ar12, file = "irfarrd10c.llo", append = TRUE)
setTxtProgressBar(pb, j)
}
close(pb)

imp50 <- matrix(0,nrow= t*M,ncol=20)
imp16 <- matrix(0,nrow= t*M,ncol=20)
imp84 <- matrix(0,nrow= t*M,ncol=20)

pb <- txtProgressBar(min = 0, max = t*M, style = 3)
for (i in 1: (t*M)){
Sys.sleep(0.1)
imp001 <- readList("irfarrd1c.llo",index=i)[[1]]
imp001 <- rbind(imp001,readList("irfarrd2c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd3c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd4c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd5c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd6c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd7c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd8c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd9c.llo",index=i)[[1]])
imp001 <- rbind(imp001,readList("irfarrd10c.llo",index=i)[[1]])

imp50[i,] <- apply(imp001, 2, quantile, probs=0.5)
imp16[i,] <- apply(imp001, 2, quantile, probs=0.16)
imp84[i,] <- apply(imp001, 2, quantile, probs=0.84)
setTxtProgressBar(pb, i)
}
close(pb)

rm(imp001)

imp16_n1 <- seq(1, ncol(imp16), by=4) # first equation
imp16_1 <- imp16[,imp16_n1]

# first to first (response of op to op)
imp16_n11 <- seq(1, nrow(imp16_1), by=4)
imp16_11 <- imp16_1[imp16_n11,]#####################

# second to first(response of xoi to op)
imp16_n12 <- seq(2, nrow(imp16_1), by=4)
imp16_12 <- imp16_1[imp16_n12,]#####################

# Third to firtst (response of epi to op)
imp16_n13 <- seq(3, nrow(imp16_1), by=4)
imp16_13 <- imp16_1[imp16_n13,]#####################



# Fourth to firtst (response of sigma to op)
imp16_n14 <- seq(4, nrow(imp16_1), by=4)
imp16_14 <- imp16_1[imp16_n14,]#####################

rm(imp16_n11,imp16_n12,imp16_n13,imp16_1,imp16_n1,imp16_n14)

imp16_n2 <- seq(2, ncol(imp16), by=4) # second equation
imp16_2 <- imp16[,imp16_n2]

# second to first (Response of op to xoi)
imp16_n21 <- seq(1, nrow(imp16_2), by=4)
imp16_21 <- imp16_2[imp16_n21,]#####################

# second to second (Response of xoi to xoi)
imp16_n22 <- seq(2, nrow(imp16_2), by=4)
imp16_22 <- imp16_2[imp16_n22,]#####################

# second to third(Response of epi to xoi)
imp16_n23 <- seq(3, nrow(imp16_2), by=4)
imp16_23 <- imp16_2[imp16_n23,]#####################

# second to (Response of sig to xoi)
imp16_n24 <- seq(4, nrow(imp16_2), by=4)
imp16_24 <- imp16_2[imp16_n24,]#####################

rm(imp16_2,imp16_n21,imp16_n22,imp16_n23,imp16_n2,imp16_n24)

imp16_n3 <- seq(3, ncol(imp16), by=4) # 3rd equation
imp16_3 <- imp16[,imp16_n3]

# third to first (Response of op to epi)
imp16_n31 <- seq(1, nrow(imp16_3), by=4)
imp16_31 <- imp16_3[imp16_n31,]#####################

# third to second  (Response of xoi to epi)
imp16_n32 <- seq(2, nrow(imp16_3), by=4)
imp16_32 <- imp16_3[imp16_n32,]#####################


# third to third  (Response of epi to epi)
imp16_n33 <- seq(3, nrow(imp16_3), by=4)
imp16_33 <- imp16_3[imp16_n33,]#####################

# third to third  (Response of sig to epi)
imp16_n34 <- seq(4, nrow(imp16_3), by=4)
imp16_34 <- imp16_3[imp16_n34,]#####################
rm(imp16_3,imp16_n31,imp16_n32,imp16_n33,imp16_n3,imp16_n34)

imp16_n4 <- seq(4, ncol(imp16), by=4) # 3rd equation
imp16_4 <- imp16[,imp16_n4]

# fort to first (Response of op to sig)
imp16_n41 <- seq(1, nrow(imp16_4), by=4)
imp16_41 <- imp16_4[imp16_n41,]#####################

# third to second  (Response of xoi to sig)
imp16_n42 <- seq(2, nrow(imp16_4), by=4)
imp16_42 <- imp16_4[imp16_n42,]#####################

# third to third  (Response of epi to sig)
imp16_n43 <- seq(3, nrow(imp16_4), by=4)
imp16_43 <- imp16_4[imp16_n43,]#####################

# third to third  (Response of sig to sig)
imp16_n44 <- seq(4, nrow(imp16_4), by=4)
imp16_44 <- imp16_4[imp16_n44,]#####################

rm(imp16_4,imp16_n41,imp16_n42,imp16_n43,imp16_n4,imp16_n44)

imp50_n1 <- seq(1, ncol(imp50), by=4) # first equation
imp50_1 <- imp50[,imp50_n1]

# first to first (response of op to op)
imp50_n11 <- seq(1, nrow(imp50_1), by=4)
imp50_11 <- imp50_1[imp50_n11,]#####################

# second to first(response of xoi to op)
imp50_n12 <- seq(2, nrow(imp50_1), by=4)
imp50_12 <- imp50_1[imp50_n12,]#####################

# Third to firtst (response of epi to op)
imp50_n13 <- seq(3, nrow(imp50_1), by=4)
imp50_13 <- imp50_1[imp50_n13,]#####################

# Fourth to firtst (response of sigma to op)
imp50_n14 <- seq(4, nrow(imp50_1), by=4)
imp50_14 <- imp50_1[imp50_n14,]#####################

rm(imp50_n11,imp50_n12,imp50_n13,imp50_1,imp50_n1,imp50_n14)

imp50_n2 <- seq(2, ncol(imp50), by=4) # second equation
imp50_2 <- imp50[,imp50_n2]

# second to first (Response of op to xoi)
imp50_n21 <- seq(1, nrow(imp50_2), by=4)
imp50_21 <- imp50_2[imp50_n21,]#####################

# second to second (Response of xoi to xoi)
imp50_n22 <- seq(2, nrow(imp50_2), by=4)
imp50_22 <- imp50_2[imp50_n22,]#####################

# second to third(Response of epi to xoi)
imp50_n23 <- seq(3, nrow(imp50_2), by=4)
imp50_23 <- imp50_2[imp50_n23,]#####################

# second to (Response of sig to xoi)
imp50_n24 <- seq(4, nrow(imp50_2), by=4)
imp50_24 <- imp50_2[imp50_n24,]#####################

rm(imp50_2,imp50_n21,imp50_n22,imp50_n23,imp50_n2,imp50_n24)

imp50_n3 <- seq(3, ncol(imp50), by=4) # 3rd equation
imp50_3 <- imp50[,imp50_n3]

# third to first (Response of op to epi)
imp50_n31 <- seq(1, nrow(imp50_3), by=4)
imp50_31 <- imp50_3[imp50_n31,]#####################

# third to second  (Response of xoi to epi)
imp50_n32 <- seq(2, nrow(imp50_3), by=4)
imp50_32 <- imp50_3[imp50_n32,]#####################

# third to third  (Response of epi to epi)
imp50_n33 <- seq(3, nrow(imp50_3), by=4)
imp50_33 <- imp50_3[imp50_n33,]#####################

# third to third  (Response of sig to epi)
imp50_n34 <- seq(4, nrow(imp50_3), by=4)
imp50_34 <- imp50_3[imp50_n34,]#####################

rm(imp50_n31,imp50_n32,imp50_n33,imp50_n3,imp50_n34)

imp50_n4 <- seq(4, ncol(imp50), by=4) # 3rd equation
imp50_4 <- imp50[,imp50_n4]

# fort to first (Response of op to sig)
imp50_n41 <- seq(1, nrow(imp50_4), by=4)
imp50_41 <- imp50_4[imp50_n41,]#####################

# third to second  (Response of xoi to sig)
imp50_n42 <- seq(2, nrow(imp50_4), by=4)
imp50_42 <- imp50_4[imp50_n42,]#####################

# third to third  (Response of epi to sig)
imp50_n43 <- seq(3, nrow(imp50_4), by=4)
imp50_43 <- imp50_4[imp50_n43,]#####################

# third to third  (Response of sig to sig)
imp50_n44 <- seq(4, nrow(imp50_3), by=4)
imp50_44 <- imp50_4[imp50_n44,]#####################

rm(imp50_4,imp50_n41,imp50_n42,imp50_n43,imp50_n4,imp50_n44)


imp84_n1 <- seq(1, ncol(imp84), by=4) # first equation
imp84_1 <- imp84[,imp84_n1]

# first to first (response of op to op)
imp84_n11 <- seq(1, nrow(imp84_1), by=4)
imp84_11 <- imp84_1[imp84_n11,]#####################

# second to first(response of xoi to op)
imp84_n12 <- seq(2, nrow(imp84_1), by=4)
imp84_12 <- imp84_1[imp84_n12,]#####################

# Third to firtst (response of epi to op)
imp84_n13 <- seq(3, nrow(imp84_1), by=4)
imp84_13 <- imp84_1[imp84_n13,]#####################

# Fourth to firtst (response of sigma to op)
imp84_n14 <- seq(4, nrow(imp84_1), by=4)
imp84_14 <- imp84_1[imp84_n14,]#####################

rm(imp84_n11,imp84_n12,imp84_n13,imp84_1,imp84_n1,imp84_n14)

imp84_n2 <- seq(2, ncol(imp84), by=4) # second equation
imp84_2 <- imp84[,imp84_n2]

# second to first (Response of op to xoi)
imp84_n21 <- seq(1, nrow(imp84_2), by=4)
imp84_21 <- imp84_2[imp84_n21,]#####################

# second to second (Response of xoi to xoi)
imp84_n22 <- seq(2, nrow(imp84_2), by=4)
imp84_22 <- imp84_2[imp84_n22,]#####################

# second to third(Response of epi to xoi)
imp84_n23 <- seq(3, nrow(imp84_2), by=4)
imp84_23 <- imp84_2[imp84_n23,]#####################

# second to (Response of sig to xoi)
imp84_n24 <- seq(4, nrow(imp84_2), by=4)
imp84_24 <- imp84_2[imp84_n24,]#####################

rm(imp84_2,imp84_n21,imp84_n22,imp84_n23,imp84_n2,imp84_n24)

imp84_n3 <- seq(3, ncol(imp84), by=4) # 3rd equation
imp84_3 <- imp84[,imp84_n3]

# third to first (Response of op to epi)
imp84_n31 <- seq(1, nrow(imp84_3), by=4)
imp84_31 <- imp84_3[imp84_n31,]#####################

# third to second  (Response of xoi to epi)
imp84_n32 <- seq(2, nrow(imp84_3), by=4)
imp84_32 <- imp84_3[imp84_n32,]#####################

# third to third  (Response of epi to epi)
imp84_n33 <- seq(3, nrow(imp84_3), by=4)
imp84_33 <- imp84_3[imp84_n33,]#####################

# third to third  (Response of sig to epi)
imp84_n34 <- seq(4, nrow(imp84_3), by=4)
imp84_34 <- imp84_3[imp84_n34,]#####################

rm(imp84_3,imp84_n31,imp84_n32,imp84_n33,imp84_n3,imp84_n34)

imp84_n4 <- seq(4, ncol(imp84), by=4) # 3rd equation
imp84_4 <- imp84[,imp84_n4]

# fort to first (Response of op to sig)
imp84_n41 <- seq(1, nrow(imp84_4), by=4)
imp84_41 <- imp84_4[imp84_n41,]#####################

# third to second  (Response of xoi to sig)
imp84_n42 <- seq(2, nrow(imp84_4), by=4)
imp84_42 <- imp84_4[imp84_n42,]#####################

# third to third  (Response of epi to sig)
imp84_n43 <- seq(3, nrow(imp84_4), by=4)
imp84_43 <- imp84_4[imp84_n43,]#####################

# third to third  (Response of sig to sig)
imp84_n44 <- seq(4, nrow(imp84_4), by=4)
imp84_44 <- imp84_4[imp84_n44,]#####################



rm(imp84_4,imp84_n41,imp84_n42,imp84_n43,imp84_n4,imp84_n44)


write.csv(imp16_11 , file = " imp16_11c.csv")
write.csv(imp16_12 , file = " imp16_12c.csv")
write.csv(imp16_13 , file =" imp16_13c.csv")
write.csv(imp16_14 , file = " imp16_14c.csv")
write.csv(imp16_21 , file = " imp16_21c.csv")
write.csv(imp16_22 , file = " imp16_22c.csv")
write.csv(imp16_23 , file = " imp16_23c.csv")
write.csv(imp16_24 , file = " imp16_24c.csv")
write.csv(imp16_31, file = " imp16_31c.csv")
write.csv(imp16_32 , file = " imp16_32c.csv")
write.csv(imp16_33 , file = " imp16_33c.csv")
write.csv(imp16_34 , file = " imp16_34c.csv")
write.csv(imp16_41 , file = " imp16_41c.csv")
write.csv(imp16_42 , file = " imp16_42c.csv")
write.csv(imp16_43 , file = " imp16_43c.csv")
write.csv(imp16_44 , file = " imp16_44c.csv")


write.csv(imp50_11 , file = " imp50_11c.csv")
write.csv(imp50_12 , file = " imp50_12c.csv")
write.csv(imp50_13 , file =" imp50_13c.csv")
write.csv(imp50_14 , file = " imp50_14c.csv")
write.csv(imp50_21 , file = " imp50_21c.csv")
write.csv(imp50_22 , file = " imp50_22c.csv")
write.csv(imp50_23 , file = " imp50_23c.csv")
write.csv(imp50_24 , file = " imp50_24c.csv")
write.csv(imp50_31, file = " imp50_31c.csv")
write.csv(imp50_32 , file = " imp50_32c.csv")
write.csv(imp50_33 , file = " imp50_33c.csv")
write.csv(imp50_34 , file = " imp50_34c.csv")
write.csv(imp50_41 , file = " imp50_41c.csv")
write.csv(imp50_42 , file = " imp50_42c.csv")
write.csv(imp50_43 , file = " imp50_43c.csv")
write.csv(imp50_44 , file = " imp50_44c.csv")


write.csv(imp84_11 , file = " imp84_11c.csv")
write.csv(imp84_12 , file = " imp84_12c.csv")
write.csv(imp84_13 , file =" imp84_13c.csv")
write.csv(imp84_14 , file = " imp84_14c.csv")
write.csv(imp84_21 , file = " imp84_21c.csv")
write.csv(imp84_22 , file = " imp84_22c.csv")
write.csv(imp84_23 , file = " imp84_23c.csv")
write.csv(imp84_24 , file = " imp84_24c.csv")
write.csv(imp84_31, file = " imp84_31c.csv")
write.csv(imp84_32 , file = " imp84_32c.csv")
write.csv(imp84_33 , file = " imp84_33c.csv")
write.csv(imp84_34 , file = " imp84_34c.csv")
write.csv(imp84_41 , file = " imp84_41c.csv")
write.csv(imp84_42 , file = " imp84_42c.csv")
write.csv(imp84_43 , file = " imp84_43c.csv")
write.csv(imp84_44 , file = " imp84_44c.csv")
################## END of Codes################################


