library(rjags)
library(textmineR)
library(matrixStats)

#Simulation study for the planned sample size of 13 patients per basket.
#Metrics considered: type I error rate, power, % all correct inference, FWER.
#Cut-off values for the mEXNEX model chosen using the c_e-analysis.R code.
#Delta values calibrated in the Delta-Calibration.R code.

#Independent Model
Ind_OC <- function(p,n,cut_ind,run){
  no.successes <- rep(0,length(p))
  hypo <- matrix(,nrow=run,ncol=length(p))
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
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate data
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N)
    jags.fit <- jags.model(file='Ind.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesInd <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE) #Fit the model
    samplesInd <- as.data.frame(samplesInd[[1]])
    pmat <- as.matrix(samplesInd[,1:length(p)])
    hypo[j,] <- as.integer(apply(pmat,2,Fun)>cut_ind) #Reject/accept the null hypothesis
    print(j)
  }
  perfect <- 0
  for(k in 1:run){
    if(all(hypo[k,]==true)){ #Correct inference across all baskets
      perfect <- perfect+1}
  }
  fwer_true <- which(true==0) #Compute FWER
  if(sum(true)==length(p)){
    fwer <- rep('NA',length(p))
  }else{
    fwer <- 0
    for(a in 1:run){
      error <- rep(0,length(fwer_true))
      for(b in 1:length(fwer_true)){
        if(hypo[a,fwer_true[b]]==1){
          error[b] <- 1
        }else{
          error[b] <- 0
        }
      }
      if(sum(error)!=0){ #Make at least one type error
        fwer <- fwer+1
      }
    }
    fwer <- fwer/run
  }
  my_list <- list('Error Rates'=colMeans(hypo),'Perfect'=perfect/run,'FWER'=fwer)
  return(my_list)
}

#BHM
BHM_OC <- function(p,n,cut_bhm,run){
  no.successes <- rep(0,length(p))
  hypo <- matrix(,nrow=run,ncol=length(p))
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
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate data
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N)
    jags.fit <- jags.model(file='BHM.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesBHM <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE) #Fit the model
    samplesBHM <- as.data.frame(samplesBHM[[1]])
    pmat <- as.matrix(samplesBHM[,1:length(p)])
    hypo[j,] <- as.integer(apply(pmat,2,Fun)>cut_bhm) #Accept/reject the null hypothesis
    print(j)
  }
  perfect <- 0
  for(k in 1:run){
    if(all(hypo[k,]==true)){ #Correct inference across all baskets
      perfect <- perfect+1}
  }
  fwer_true <- which(true==0) #Compute FWER
  if(sum(true)==length(p)){
    fwer <- rep('NA',length(p))
  }else{
    fwer <- 0
    for(a in 1:run){
      error <- rep(0,length(fwer_true))
      for(b in 1:length(fwer_true)){
        if(hypo[a,fwer_true[b]]==1){
          error[b] <- 1
        }else{
          error[b] <- 0
        }
      }
      if(sum(error)!=0){ #Make at least one type I error
        fwer <- fwer+1
      }
    }
    fwer <- fwer/run
  }
  my_list <- list('Error Rates'=colMeans(hypo),'Perfect'=perfect/run,'FWER'=fwer)
  return(my_list)
}

#CBHM
CBHM_OC <- function(p,n,cut_cbhm,run,a,b){
  no.successes <- rep(0,length(p))
  hypo <- matrix(,nrow=run,ncol=length(p))
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
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate data
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    N <- length(p)
    O0 <- n-y #Observed failures
    O1 <- y #Observed successes
    E0 <- c() #Expected failures
    E1 <- c() #Expected successes
    for(k in 1:length(n)){
    E0[k] <- n[k]*((sum(n)-sum(y))/sum(n))
    E1[k] <- n[k]*((sum(y))/sum(n))}
    if(sum(E1==0)){
      Test <- 0
    }else{
      first <- c()
      last <- c()
      for(o in 1:length(p)){
        first[o] <- (O0[o]-E0[o])^2/E0[o]
        last[o] <- (O1[o]-E1[o])^2/E1[o]
      }
      Test <- sum(first)+sum(last) #Chi squared test statistic for homogeneity 
    }
    if(Test!=0){
      sigma <- exp(a+b*log(Test))
      jags.data <- list('n'=n,'y'=y,'N'=N,'sigma2'=sigma)
      jags.fit <- jags.model(file='CBHM.txt',data=jags.data,n.adapt=1000,n.chains=1)
      samplesCBHM <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE)
      samplesCBHM <- as.data.frame(samplesCBHM[[1]])
      pmat <- as.matrix(samplesCBHM[,1:length(p)])}
    else{
      jags.data <- list('n'=n,'y'=y,'N'=N)
      jags.fit <- jags.model(file='CBHMNoSigma.txt',data=jags.data,n.adapt=1000,n.chains=1)
      samplesCBHM <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE)
      samplesCBHM <- as.data.frame(samplesCBHM[[1]])
      pmat <- as.matrix(samplesCBHM[,1:length(p)])}
    hypo[j,] <- as.integer(apply(pmat,2,Fun)>cut_cbhm) #Accept/reject the null hypothesis
    print(j)
  }
  perfect <- 0
  for(k in 1:run){
    if(all(hypo[k,]==true)){ #Correct inference across all baskets
      perfect <- perfect+1}
  }
  fwer_true <- which(true==0) #Calculate FWER
  if(sum(true)==length(p)){
    fwer <- rep('NA',length(p))
  }else{
    fwer <- 0
    for(a in 1:run){
      error <- rep(0,length(fwer_true))
      for(b in 1:length(fwer_true)){
        if(hypo[a,fwer_true[b]]==1){
          error[b] <- 1
        }else{
          error[b] <- 0
        }
      }
      if(sum(error)!=0){ #Make at least one type I error
        fwer <- fwer+1
      }
    }
    fwer <- fwer/run
  }
  my_list <- list('Error Rates'=colMeans(hypo),'Perfect'=perfect/run,'FWER'=fwer)
  return(my_list)
}

#EXNEX
EXNEX_OC <- function(p,n,cut_exnex,run,pw){
  no.successes <- rep(0,length(p))
  hypo <- matrix(,nrow=run,ncol=length(p))
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
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate data
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    nexmu <- rep(log(pw/(1-pw)),length(p)) #NEX mu parameter
    nexsigma <- rep((1/pw)+(1/(1-pw)),length(p)) #NEX sigma parameter
    prob <- c(0.5,0.5) #Fixed mixture weights
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N,'nexmu'=nexmu,'nexsigma'=nexsigma,'prob'=prob)
    jags.fit <- jags.model(file='EXNEX.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesEXNEX <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE)
    samplesEXNEX <- as.data.frame(samplesEXNEX[[1]])
    pmat <- as.matrix(samplesEXNEX[,1:length(p)])
    hypo[j,] <- as.integer(apply(pmat,2,Fun)>cut_exnex) #Accept/reject the null hypothesis
    print(j)
  }
  perfect <- 0
  for(k in 1:run){
    if(all(hypo[k,]==true)){ #Correct inference across all baskets
      perfect <- perfect+1}
  }
  fwer_true <- which(true==0) #Calculate FWER
  if(sum(true)==length(p)){ 
    fwer <- rep('NA',length(p))
  }else{
    fwer <- 0
    for(a in 1:run){
      error <- rep(0,length(fwer_true))
      for(b in 1:length(fwer_true)){
        if(hypo[a,fwer_true[b]]==1){
          error[b] <- 1
        }else{
          error[b] <- 0
        }
      }
      if(sum(error)!=0){ #Make at least one type I error
        fwer <- fwer+1
      }
    }
    fwer <- fwer/run
  }
  my_list <- list('Error Rates'=colMeans(hypo),'Perfect'=perfect/run,'FWER'=fwer)
  return(my_list)
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

mEXNEX_OC <- function(p,n,cut_mexnex,run,pw,c){
  no.successes <- rep(0,length(p))
  hypo <- matrix(,nrow=run,ncol=length(p))
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
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate data
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    nexmu <- rep(log(pw/(1-pw)),length(p)) #NEX mu parameter
    nexsigma <- rep((1/pw)+(1/(1-pw)),length(p)) #NEX sigma parameter
    prob <- pi_fun_H(y,n,c) #Compute EX/NEX probabilities using averaged Hellinger distances
    prob <- cbind(prob,1-prob)
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N,'nexmu'=nexmu,'nexsigma'=nexsigma,'prob'=prob)
    jags.fit <- jags.model(file='mEXNEX.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesmEXNEX <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE)
    samplesmEXNEX <- as.data.frame(samplesmEXNEX[[1]])
    pmat <- as.matrix(samplesmEXNEX[,1:length(p)])
    hypo[j,] <- as.integer(apply(pmat,2,Fun)>cut_mexnex) #Accept/reject the null hypothesis
    print(j)
  }
  perfect <- 0
  for(k in 1:run){
    if(all(hypo[k,]==true)){ #Correct inference across all baskets
      perfect <- perfect+1}
  }
  fwer_true <- which(true==0) #Calculated FWER
  if(sum(true)==length(p)){
    fwer <- rep('NA',length(p))
  }else{
    fwer <- 0
    for(a in 1:run){
      error <- rep(0,length(fwer_true))
      for(b in 1:length(fwer_true)){
        if(hypo[a,fwer_true[b]]==1){
          error[b] <- 1
        }else{
          error[b] <- 0
        }
      }
      if(sum(error)!=0){ #Make at least one type I error
        fwer <- fwer+1
      }
    }
    fwer <- fwer/run
  }
  my_list <- list('Error Rates'=colMeans(hypo),'Perfect'=perfect/run,'FWER'=fwer)
  return(my_list)
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

mEXNEXmin_OC <- function(p,n,cut_mexnexmin,run,pw){
  no.successes <- rep(0,length(p))
  hypo <- matrix(,nrow=run,ncol=length(p))
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
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate data
    M <- 10000
    Fun <- function(x){
      fun <- sum(x>0.15)/M
      return(fun)
    }
    y <- no.successes
    nexmu <- rep(log(pw/(1-pw)),length(p)) #NEX mu parameter
    nexsigma <- rep((1/pw)+(1/(1-pw)),length(p)) #NEX sigma parameter
    prob <- pi_fun_H_min(y,n) #Compute EX/NEX probabilities using minimum Hellinger distances
    prob <- cbind(prob,1-prob)
    N <- length(p)
    jags.data <- list('n'=n,'y'=y,'N'=N,'nexmu'=nexmu,'nexsigma'=nexsigma,'prob'=prob)
    jags.fit <- jags.model(file='mEXNEX.txt',data=jags.data,n.adapt=1000,n.chains=1)
    samplesmEXNEXmin <- coda.samples(jags.fit,variable.names = c('p'),n.iter=M,silent=TRUE)
    samplesmEXNEXmin <- as.data.frame(samplesmEXNEXmin[[1]])
    pmat <- as.matrix(samplesmEXNEXmin[,1:length(p)])
    hypo[j,] <- as.integer(apply(pmat,2,Fun)>cut_mexnexmin) #Accept/reject the null hypothesis
    print(j)
  }
  perfect <- 0
  for(k in 1:run){
    if(all(hypo[k,]==true)){ #Correct inference across all baskets
      perfect <- perfect+1}
  }
  fwer_true <- which(true==0) #Calculate FWER
  if(sum(true)==length(p)){
    fwer <- rep('NA',length(p))
  }else{
    fwer <- 0
    for(a in 1:run){
      error <- rep(0,length(fwer_true))
      for(b in 1:length(fwer_true)){
        if(hypo[a,fwer_true[b]]==1){
          error[b] <- 1
        }else{
          error[b] <- 0
        }
      }
      if(sum(error)!=0){ #Make at least one type I error
        fwer <- fwer+1
      }
    }
    fwer <- fwer/run
  }
  my_list <- list('Error Rates'=colMeans(hypo),'Perfect'=perfect/run,'FWER'=fwer)
  return(my_list)
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

BMA_OC <- function(p,n,run,piA,cut){
  no.successes <- rep(0,length(p))
  hypo <- matrix(,nrow=run,ncol=length(p))
  true <- rep(0,length(p))
  for(l in 1:length(p)){
    if(p[l]<=0.2){
      true[l] <- 0
    }else{
      true[l] <- 1
    }
  }
  for(j in 1:run){
    for(i in 1:length(p)){
      no.successes[i] <- rbinom(1,n[i],p[i])} #Generate data
    mod <- model_averaging(no.successes,n,piA,0.15) #Conduct Bayesian model averaging
    hypo[j,] <- as.numeric(mod$`Decision Probabilities`>cut_bma) #Accept/reject the null hypothesis
  }
  perfect <- 0
  for(k in 1:run){
    if(all(hypo[k,]==true)){ #Correct inference across all baskets
      perfect <- perfect+1}
  }
  fwer_true <- which(true==0) #Calculate FWER
  if(sum(true)==length(p)){
    fwer <- rep('Na',length(p))
  }else{
    fwer <- 0
    for(a in 1:run){
      error <- rep(0,length(fwer_true))
      for(b in 1:length(fwer_true)){
        if(hypo[a,fwer_true[b]]==1){
          error[b] <- 1
        }else{
          error[b] <- 0
        }
      }
      if(sum(error)!=0){ #Make at least one type I error
        fwer <- fwer+1
      }
    }
    fwer <- fwer/run
  }
  my_list <- list('Error Rates'=colMeans(hypo),'Perfect'=perfect/run,'FWER'=fwer)
  return(my_list)
}

#Inputs
p1 <- c(0.15,0.15,0.15,0.15,0.15)
p2 <- c(0.45,0.15,0.15,0.15,0.15)
p3 <- c(0.45,0.45,0.15,0.15,0.15)
p4 <- c(0.45,0.45,0.45,0.15,0.15)
p5 <- c(0.45,0.45,0.45,0.45,0.15)
p6 <- c(0.45,0.45,0.45,0.45,0.45)
n <- rep(13,5)
pw <- rep(0.35,5)
piA <- 0.45
run <- 10000
a <- -7.248794
b <- 5.857633

#Cut-Off Values 
#Obtained using the Delta-Calibration.R code
cut_ind <- 0.951162
cut_bhm <- 0.80945
cut_cbhm <- 0.894508
cut_exnex <- 0.877432
cut_mexnexc <- 0.919946
cut_mexnex <- 0.881284
cut_mexnexmin <- 0.839586
cut_bma <- 0.8709033

#Scenario 1
IndSc1 <- Ind_OC(p1,n,cut_ind,run)
BHMSc1 <- BHM_OC(p1,n,cut_bhm,run)
CBHMSc1 <- CBHM_OC(p1,n,cut_cbhm,run,a,b)
EXNEXSc1 <- EXNEX_OC(p1,n,cut_exnex,run,pw)
mEXNEXcSc1 <- mEXNEX_OC(p1,n,cut_mexnexc,run,pw,0.05)
mEXNEXSc1 <- mEXNEX_OC(p1,n,cut_mexnex,run,pw,0.1)
mEXNEXminSc1 <- mEXNEXmin_OC(p1,n,cut_mexnexmin,run,pw)
BMASc1 <- BMA_OC(p1,n,run,piA,cut_bma)

#Scenario 2
IndSc2 <- Ind_OC(p2,n,cut_ind,run)
BHMSc2 <- BHM_OC(p2,n,cut_bhm,run)
CBHMSc2 <- CBHM_OC(p2,n,cut_cbhm,run,a,b)
EXNEXSc2 <- EXNEX_OC(p2,n,cut_exnex,run,pw)
mEXNEXcSc2 <- mEXNEX_OC(p2,n,cut_mexnexc,run,pw,0.05)
mEXNEXSc2 <- mEXNEX_OC(p2,n,cut_mexnex,run,pw,0.1)
mEXNEXminSc2 <- mEXNEXmin_OC(p2,n,cut_mexnexmin,run,pw)
BMASc2 <- BMA_OC(p2,n,run,piA,cut_bma)

#Scenario 3
IndSc3 <- Ind_OC(p3,n,cut_ind,run)
BHMSc3 <- BHM_OC(p3,n,cut_bhm,run)
CBHMSc3 <- CBHM_OC(p3,n,cut_cbhm,run,a,b)
EXNEXSc3 <- EXNEX_OC(p3,n,cut_exnex,run,pw)
mEXNEXcSc3 <- mEXNEX_OC(p3,n,cut_mexnexc,run,pw,0.05)
mEXNEXSc3 <- mEXNEX_OC(p3,n,cut_mexnex,run,pw,0.1)
mEXNEXminSc3 <- mEXNEXmin_OC(p3,n,cut_mexnexmin,run,pw)
BMASc3 <- BMA_OC(p3,n,run,piA,cut_bma)

#Scenario 4
IndSc4 <- Ind_OC(p4,n,cut_ind,run)
BHMSc4 <- BHM_OC(p4,n,cut_bhm,run)
CBHMSc4 <- CBHM_OC(p4,n,cut_cbhm,run,a,b)
EXNEXSc4 <- EXNEX_OC(p4,n,cut_exnex,run,pw)
mEXNEXcSc4 <- mEXNEX_OC(p4,n,cut_mexnexc,run,pw,0.05)
mEXNEXSc4 <- mEXNEX_OC(p4,n,cut_mexnex,run,pw,0.1)
mEXNEXminSc4 <- mEXNEXmin_OC(p4,n,cut_mexnexmin,run,pw)
BMASc4 <- BMA_OC(p4,n,run,piA,cut_bma)

#Scenario 5
IndSc5 <- Ind_OC(p5,n,cut_ind,run)
BHMSc5 <- BHM_OC(p5,n,cut_bhm,run)
CBHMSc5 <- CBHM_OC(p5,n,cut_cbhm,run,a,b)
EXNEXSc5 <- EXNEX_OC(p5,n,cut_exnex,run,pw)
mEXNEXcSc5 <- mEXNEX_OC(p5,n,cut_mexnexc,run,pw,0.05)
mEXNEXSc5 <- mEXNEX_OC(p5,n,cut_mexnex,run,pw,0.1)
mEXNEXminSc5 <- mEXNEXmin_OC(p5,n,cut_mexnexmin,run,pw)
BMASc5 <- BMA_OC(p5,n,run,piA,cut_bma)

#Scenario 6
IndSc6 <- Ind_OC(p6,n,cut_ind,run)
BHMSc6 <- BHM_OC(p6,n,cut_bhm,run)
CBHMSc6 <- CBHM_OC(p6,n,cut_cbhm,run,a,b)
EXNEXSc6 <- EXNEX_OC(p6,n,cut_exnex,run,pw)
mEXNEXcSc6 <- mEXNEX_OC(p6,n,cut_mexnexc,run,pw,0.05)
mEXNEXSc6 <- mEXNEX_OC(p6,n,cut_mexnex,run,pw,0.1)
mEXNEXminSc6 <- mEXNEXmin_OC(p6,n,cut_mexnexmin,run,pw)
BMASc6 <- BMA_OC(p6,n,run,piA,cut_bma)
