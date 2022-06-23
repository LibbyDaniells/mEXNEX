library(rjags)
library(textmineR)
library(matrixStats)

#This code computes operating characteristics for the mEXNEX model under different cut-off values c_e
#Metrics computed: type I error rate, power, % all correct inference, FWER
#Computations done under planned sample size with null and target response rates

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

c <- c(0,1/13,2/13,3/13,4/13)
p <- rep(0.15,5)
n <- rep(13,5)
pw <- rep(0.35,5)
run <- 10000
cut <- matrix(,nrow=length(c),ncol=length(p))
for(i in 1:length(c)){
  cut[i,] <- mEXNEX_cutoff(p,n,run,pw,c[i])
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

p1 <- rep(0.15,5)
p2 <- c(0.45,0.15,0.15,0.15,0.15)
p3 <- c(0.45,0.45,0.15,0.15,0.15)
p4 <- c(0.45,0.45,0.45,0.15,0.15)
p5 <- c(0.45,0.45,0.45,0.45,0.15)
p6 <- c(0.45,0.45,0.45,0.45,0.45)
n <- rep(13,5)
pw <- rep(0.35,5)
run <- 10000

#Scenario 1
Reject <- matrix(,nrow=length(c),ncol=length(p1))
Perfect <- c()
FWER <- c()
for(i in 1:length(c)){
  OpCar <- mEXNEX_OC(p1,n,cut[i,],run,pw,c[i])
  Reject[i,] <- OpCar$`Error Rates`
  Perfect[i] <- OpCar$Perfect
  FWER[i] <- OpCar$FWER
}
Sc1 <- data.frame(Cuts=cut,Reject=Reject,Perfect=Perfect,FWER=FWER)

#Scenario 2
Reject <- matrix(,nrow=length(c),ncol=length(p1))
Perfect <- c()
FWER <- c()
for(i in 1:length(c)){
  OpCar <- mEXNEX_OC(p2,n,cut[i,],run,pw,c[i])
  Reject[i,] <- OpCar$`Error Rates`
  Perfect[i] <- OpCar$Perfect
  FWER[i] <- OpCar$FWER
}
Sc2 <- data.frame(Cuts=cut,Reject=Reject,Perfect=Perfect,FWER=FWER)

#Scenario 3
Reject <- matrix(,nrow=length(c),ncol=length(p1))
Perfect <- c()
FWER <- c()
for(i in 1:length(c)){
  OpCar <- mEXNEX_OC(p3,n,cut[i,],run,pw,c[i])
  Reject[i,] <- OpCar$`Error Rates`
  Perfect[i] <- OpCar$Perfect
  FWER[i] <- OpCar$FWER
}
Sc3 <- data.frame(Cuts=cut,Reject=Reject,Perfect=Perfect,FWER=FWER)

#Scenario 4
Reject <- matrix(,nrow=length(c),ncol=length(p1))
Perfect <- c()
FWER <- c()
for(i in 1:length(c)){
  OpCar <- mEXNEX_OC(p4,n,cut[i,],run,pw,c[i])
  Reject[i,] <- OpCar$`Error Rates`
  Perfect[i] <- OpCar$Perfect
  FWER[i] <- OpCar$FWER
}
Sc4 <- data.frame(Cuts=cut,Reject=Reject,Perfect=Perfect,FWER=FWER)

#Scenario 5
Reject <- matrix(,nrow=length(c),ncol=length(p1))
Perfect <- c()
FWER <- c()
for(i in 1:length(c)){
  OpCar <- mEXNEX_OC(p5,n,cut[i,],run,pw,c[i])
  Reject[i,] <- OpCar$`Error Rates`
  Perfect[i] <- OpCar$Perfect
  FWER[i] <- OpCar$FWER
}
Sc5 <- data.frame(Cuts=cut,Reject=Reject,Perfect=Perfect,FWER=FWER)

#Scenario 6
Reject <- matrix(,nrow=length(c),ncol=length(p1))
Perfect <- c()
for(i in 1:length(c)){
  OpCar <- mEXNEX_OC(p6,n,cut[i,],run,pw,c[i])
  Reject[i,] <- OpCar$`Error Rates`
  Perfect[i] <- OpCar$Perfect
}
Sc6 <- data.frame(Cuts=cut,Reject=Reject,Perfect=Perfect)

