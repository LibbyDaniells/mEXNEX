library(rjags)
library(textmineR)
library(matrixStats)

#This code calibrates the cut-off Delta, for accepting/rejecting the null at the conclusion of the trial, to control error rates
#Target type I error rate = 10%
#Each method is calibrated seperately under the null


#Independent Model
Ind_cutoff <- function(p,n,run,pw){
  no.successes <- rep(0,length(p))
  cut <- matrix(,nrow=run,ncol=length(p))
  true <- rep(0,length(p))
  for(l in 1:length(p)){ 
    if(p[l]<=0.15){ 
      true[l] <- 0
    }else{
      true[l] <- 1
    }
  }
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate Data
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    N <- length(y)
    mu <- log(pw/(1-pw))
    jags.data <- list('n'=n,'y'=y,'N'=N,'mu'=mu)
    jags.fit <- jags.model(file='Ind.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesEXNEX <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE) #Fit the model
    samplesEXNEX <- as.data.frame(samplesEXNEX[[1]])
    pmat <- as.matrix(samplesEXNEX[,1:length(p)])
    cut[j,] <- apply(pmat,2,Fun)
    print(j)
  }
  cut_off <- colQuantiles(cut,probs=0.9) #Calibrating to ensure 10% type I error under the null
  return(cut_off)
}

#BHM
BHM_cutoff <- function(p,n,run){
  no.successes <- rep(0,length(p))
  cut <- matrix(,nrow=run,ncol=length(p))
  true <- rep(0,length(p))
  for(l in 1:length(p)){
    if(p[l]<=0.15){
      true[l] <- 0
    }else{
      true[l] <- 1
    }
  }
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])}
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    N <- length(y)
    jags.data <- list('n'=n,'y'=y,'N'=N)
    jags.fit <- jags.model(file='BHM.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesBHM <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE)
    samplesBHM <- as.data.frame(samplesBHM[[1]])
    pmat <- as.matrix(samplesBHM[,1:length(p)])
    cut[j,] <- apply(pmat,2,Fun)
    print(j)
  }
  cut_off <- colQuantiles(cut,probs=0.9) #Calibrating to ensure 10% type I error under the null
  return(cut_off)
}

#CBHM
CBHM_cutoff <- function(p,n,run,a,b){
  no.successes <- rep(0,length(p))
  cut <- matrix(,nrow=run,ncol=length(p))
  true <- rep(0,length(p))
  for(l in 1:length(p)){
    if(p[l]<=0.15){
      true[l] <- 0
    }else{
      true[l] <- 1
    }
  }
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])}
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    N <- length(p)
    phat <- sum(y)/sum(n)
    obs <- cbind(y,n-y)
    E <- cbind(n*phat,n*(1-phat))
    Test <- sum((abs(obs-E))^2/E)
    if(is.nan(Test)|(Test<1)){Test <- 1}
    sigma <- exp(a+b*log(Test))
      jags.data <- list('n'=n,'y'=y,'N'=N,'sigma2'=sigma)
      jags.fit <- jags.model(file='CBHM.txt',data=jags.data,n.adapt=1000,n.chains=1)
      samplesCBHM <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE)
      samplesCBHM <- as.data.frame(samplesCBHM[[1]])
      pmat <- as.matrix(samplesCBHM[,1:length(p)])}
    cut[j,] <- apply(pmat,2,Fun)
    print(j)
  }
  cut_off <- colQuantiles(cut,probs=0.9) #Calibrating to ensure 10% type I error under the null
  return(cut_off)
}

#EXNEX
EXNEX_cutoff <- function(p,n,run,pw){
  no.successes <- rep(0,length(p))
  cut <- matrix(,nrow=run,ncol=length(p))
  true <- rep(0,length(p))
  for(l in 1:length(p)){
    if(p[l]<=0.15){
      true[l] <- 0
    }else{
      true[l] <- 1
    }
  }
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])}
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    nexmu <- rep(log(pw/(1-pw)),length(p)) #NEX mu parameter
    nexsigma <- rep((1/pw)+(1/(1-pw)),length(p)) #NEX sigma parameter
    prob <- c(0.5,0.5) #Fixed weights
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N,'nexmu'=nexmu,'nexsigma'=nexsigma,'prob'=prob)
    jags.fit <- jags.model(file='EXNEX.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesEXNEX <- coda.samples(jags.fit,variable.names = c('p','weight'),n.iter=M,silent=TRUE)
    samplesEXNEX <- as.data.frame(samplesEXNEX[[1]])
    pmat <- as.matrix(samplesEXNEX[,1:length(p)])
    cut[j,] <- apply(pmat,2,Fun)
    print(j)
  }
  cut_off <- colQuantiles(cut,probs=0.9) #Calibrating to ensure 10% type I error under the null
  return(cut_off)
}

#mEXNEX
diff <- function(vec,n){  #Compute pairwise differences in response rates
  len <- 1:length(vec)
  diff_mat <- matrix(,nrow=length(vec),ncol=length(vec))
  for(i in 1:length(vec)){
    for(j in len[-i]){
      diff_mat[i,j] <- abs(vec[i]-vec[j])
    }
  }
  return(diff_mat)
}

pi_fun_H <- function(vec,n,c){
  data <- data.frame(Responses=vec,No.Patients=n)
  y <- data$Responses
  x <- seq(0,1,by=0.0001)
  p <- matrix(,nrow=length(vec),ncol=length(x)) #Get posterior densities
  for(i in 1:length(vec)){
    p[i,] <- dbeta(x,y[i]+1,n[i]-y[i]+1)
  }
  x <- p
  est <- vec/n
  diff_mat <- diff(est,n)
  rem <- c()
  keep <- c()
  pi <- rep('NA',length(vec))
  for(i in 1:length(vec)){
    min <- as.numeric(min(diff_mat[i,],na.rm=TRUE))
    if(min>(c+0.00001)){ #If min pairwise difference is greater than c_e then treat as independent
      rem[i] <- i
      pi[i] <- 0
    }else{
      keep[i] <- i
    }
  }
  keep <- keep[!is.na(keep)]
  rem <- rem[!is.na(rem)]
  if(length(rem)==length(vec)){ #If all 'extreme' treat all independent
    pi <- rep(0,length(vec))
  }else if(length(rem)==0){ #If all not 'extreme' compute all pairwise Hellinger distances and take average
    mat <- 1-CalcHellingerDist(x)
    for(i in 1:length(vec)){
      pi[i] <- mean(mat[i,-i])
    }
  }else{ #Compute pairwise Hellinger distances between those not deemed 'extreme' and treat the rest independent
    y <- x[-rem,]
    mat2 <- 1-CalcHellingerDist(y)
    for(i in 1:length(keep)){
      pi[keep[i]] <- mean(mat2[i,-i])
    }
  }
  pi <- as.numeric(pi) #Probability vector for the EXNEX model
  return(pi)
}

mEXNEX_cutoff <- function(p,n,run,pw,c){
  no.successes <- rep(0,length(p))
  cut <- matrix(,nrow=run,ncol=length(p))
  true <- rep(0,length(p))
  for(l in 1:length(p)){
    if(p[l]<=0.15){
      true[l] <- 0
    }else{
      true[l] <- 1
    }
  }
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])}
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    nexmu <- rep(log(pw/(1-pw)),length(p)) #NEX mu parameter
    nexsigma <- rep((1/pw)+(1/(1-pw)),length(p)) #NEX sigma parameter
    prob <- pi_fun_H(y,n,c) #Compute the probability vector using the Hellinger distances
    prob <- cbind(prob,1-prob)
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N,'nexmu'=nexmu,'nexsigma'=nexsigma,'prob'=prob)
    jags.fit <- jags.model(file='mEXNEX.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesmEXNEX <- coda.samples(jags.fit,variable.names = c('p','weight'),n.iter=M,silent=TRUE)
    samplesmEXNEX <- as.data.frame(samplesmEXNEX[[1]])
    pmat <- as.matrix(samplesmEXNEX[,1:length(p)])
    cut[j,] <- apply(pmat,2,Fun)
    print(j)
  }
  cut_off <- colQuantiles(cut,probs=0.9) #Calibrating to ensure 10% type I error under the null
  return(cut_off)
}

#mEXNEXmin
pi_fun_H_min <- function(vec,n){
  data <- data.frame(Responses=vec,No.Patients=n)
  N <- dim(data)[1]
  K <- data$No.Patients
  y <- data$Responses
  M <- 100000
  x <- seq(0,1,by=0.0001)
  p <- matrix(,nrow=length(vec),ncol=length(x))
  for(i in 1:length(vec)){
    p[i,] <- dbeta(x,y[i]+1,n[i]-y[i]+1) #Get posterior densities
  }
  x <- p
  Hell <- 1-CalcHellingerDist(x) 
  pi <- c()
  for(j in 1:length(vec)){
    pi[j] <- max(Hell[j,-j]) #Get minimum Hellinger distance
  }
  return(pi)
}

mEXNEXmin_cutoff <- function(p,n,run,pw){
  no.successes <- rep(0,length(p))
  cut <- matrix(,nrow=run,ncol=length(p))
  true <- rep(0,length(p))
  for(l in 1:length(p)){
    if(p[l]<=0.15){
      true[l] <- 0
    }else{
      true[l] <- 1
    }
  }
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])}
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    nexmu <- rep(log(pw/(1-pw)),length(p)) #NEX mu parameter
    nexsigma <- rep((1/pw)+(1/(1-pw)),length(p)) #NEX sigma parameter
    prob <- pi_fun_H_min(y,n) #Probability vector using minimum Hellinger distances
    prob <- cbind(prob,1-prob)
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N,'nexmu'=nexmu,'nexsigma'=nexsigma,'prob'=prob)
    jags.fit <- jags.model(file='mEXNEX.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesmEXNEXmin <- coda.samples(jags.fit,variable.names = c('p','weight'),n.iter=M,silent=TRUE)
    samplesmEXNEXmin <- as.data.frame(samplesmEXNEXmin[[1]])
    pmat <- as.matrix(samplesmEXNEXmin[,1:length(p)])
    cut[j,] <- apply(pmat,2,Fun)
    print(j)
  }
  cut_off <- colQuantiles(cut,probs=0.9)  #Calibrating to ensure 10% type I error under the null
  return(cut_off)
}

#BMA
model_averaging <- function(y,n,piA,cut){
  a0 <- piA
  b0 <- 1-a0
  perm <- expand.grid(rep(list(0:1),length(y)))
  rem <- which(apply(perm,1,sum)==1)
  perm <- perm[-rem,]
  perm <- as.matrix(perm) #Possible models
  dimM <- dim(perm)[1]
  omega <- matrix('NA',nrow=dimM,ncol=length(y)) 
  for(i in 1:dimM){ #Split of baskets into distinct response rates
    vec <- perm[i,]
    ind <- which(vec==0)
    if(sum(vec)==0){
      vec <- c(1:length(vec))
    }else{
      for(j in 1:length(ind)){
        vec[ind[j]] <- j+1
      }
    }
    omega[i,] <- vec
  }
  P <- apply(omega,1,function(x) sum(!duplicated(x))) #No. distinct response rates in each model
  post <- as.list(numeric(dimM*length(y))) #Posterior distributions
  dim(post) <- c(dimM,length(y))
  for(i in 1:dimM){
    for(j in 1:P[i]){
      id <- which(omega[i,]==j)
      a <- a0+sum(y[id])
      b <- b0+sum(n[id]-y[id])
      post[[i,j]] <- c(a,b)
    }
    for(k in 1:length(y)){
      if(sum(post[[i,k]]==0)==1){post[[i,k]]<- 'NA'}
    }
  }
  marg <- c() #Marginal likelihoods
  for(i in 1:dimM){
    const <- prod(choose(n,y))
    betaf <- c()
    for(j in 1:P[i]){
      ajp <- post[[i,j]][1]
      bjp <- post[[i,j]][2]
      betaf[j] <- beta(ajp,bjp)
    }
    marg[i] <- prod(betaf/beta(a0,b0))*const
  }
  PostM <- c() #Posterior probability of the models
  priorM <- P^2
  tot <- sum(marg*priorM)
  for(i in 1:dimM){
    PostM[i] <- (marg[i]*priorM[i])/tot
  }
  DecisionProb <- c()  #Decision probabilities
  for(i in 1:length(y)){
    Om <- omega[,i]
    Om <- as.numeric(Om)
    cdf <- c()
    for(j in 1:dimM){
      a <- post[[j,Om[j]]][1]
      b <- post[[j,Om[j]]][2]
      cdf[j] <- 1-pbeta(cut,a,b)
    }
    DecisionProb[i] <- sum(cdf*PostM)
  }
  mean <- c()
  var <- c()
  for(i in 1:length(y)){
    Om <- omega[,i]  #Which group is basket 1 in 
    Om <- as.numeric(Om)
    meanM <- c()
    varM <- c()
    for(j in 1:dim(perm)[1]){
      a <- post[[j,Om[j]]][1]  #alpha value for corresponding group
      b <- post[[j,Om[j]]][2]  #beta value for corresponding group
      meanM[j] <- a/(a+b) #mean 
      varM[j] <- (a*b)/((a+b+1)*(a+b)^2) 
    }
    mean[i] <- sum(meanM*PostM)
    var[i] <- sum((varM+meanM^2)*PostM)-mean[i]^2
  }
  sd <- sqrt(var)
  my_list <- list('Models'=perm,'Omega'=omega,'No. Distinct Responses'=P,'Posterior Distribution Parameters'=post,
                  'Marginal Likelihood'=marg,'Posterior probabiliity of the models'=PostM,'Decision Probabilities'=DecisionProb,'Mean'=mean,'Sd'=sd)
  return(my_list)
}

BMA_cutoff <- function(p,n,run,piA){
  cut <- matrix(,nrow=run,ncol=length(p))
  no.successes <- rep(0,length(p))
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])}
    mod <- model_averaging(no.successes,n,piA,0.15)
    cut[j,] <- mod$'Decision Probabilities'
  }
  cut <- colQuantiles(cut,probs=0.9)  #Calibrating to ensure 10% type I error under the null
  return(cut)
}


#Simulation for planned sample size of n_k=13
#Data Input
p <- rep(0.15,5)
n <- rep(13,5)
pw <- rep(0.35,5)
piA <- 0.45
run <- 10000
a <- -7.248794
b <- 5.857633

IndCut <- Ind_cutoff(p,n,run,pw=0.35)
BHMCut <- BHM_cutoff(p,n,run)
CBHMCut <- CBHM_cutoff(p,n,run,a,b)
EXNEXCut <- EXNEX_cutoff(p,n,run,pw)
mEXNEXCut <- mEXNEX_cutoff(p,n,run,pw,1/13)
mEXNEXcCut <- mEXNEX_cutoff(p,n,run,pw,0)
mEXNEXminCut <- mEXNEXmin_cutoff(p,n,run,pw)
BMACut <- BMA_cutoff(p,n,run,piA)

