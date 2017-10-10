
p2model <- function (net, sender = NULL, receiver = NULL , density = NULL, reciprocity = NULL) 
{
  y <- net[row(net)!=col(net)]
  nact <- ncol(net)
  n <- nact*(nact-1)
  if(!is.null(sender)){
    ns <- dim(model.matrix(sender))[2]-1
    XS <- matrix(model.matrix(sender)[,2:(ns+1)])
    X1 <- matrix(rep(NA, n*ns), ncol=ns)
    for (i in 1:ns){
      S <- XS[,i]%*%t(rep(1,nact))
      X1[,i] <- S[row(S)!=col(S)]
    }
  }  else {
    ns <- 0
    X1 <- NULL
  }
  if(!is.null(receiver)){
    nre <- dim(model.matrix(receiver))[2]-1
    XRE <- matrix(model.matrix(receiver)[,2:(nre+1)])
    X2 <- matrix(rep(NA, n*nre), ncol=nre)
    for (i in 1:nre){
      RE <- rep(1,nact)%*%t(XRE[,i])
      X2[,i] <- RE[row(RE)!=col(RE)]
    }
  } else {
    nre <- 0
    X2 <- NULL
  }  
  if(!is.null(density)){
    nd <- 1+ ((dim(model.matrix(density))[2]-1)/nact)
    X3 <- matrix(rep(NA, n*nd), ncol=nd)
    X3[,1] <- 1
    for (i in 2:nd){
      Dtmp <- model.matrix(density)[,(2+(i-2)*nact): (1+(i-1)*nact)]
      X3[,i] <- Dtmp[row(Dtmp)!=col(Dtmp)]
    }
  } else {
    nd <- 1
    X3 <- matrix(rep(1, n*nd), ncol=nd)
  }  
  XC1 <- matrix(rep(NA, n*nact), ncol=nact)
  XC2 <- matrix(rep(NA, n*nact), ncol=nact)
  for (i in 1:nact){
    SC <- matrix(rep(0, nact*nact), ncol=nact)
    SC[i,] <-1
    XC1[,i] <- SC[row(SC)!=col(SC)]
    REC <- matrix(rep(0, nact*nact), ncol=nact)
    REC[,i] <-1
    XC2[,i] <- REC[row(REC)!=col(REC)]
  }
  X <- cbind(X1, X2, X3, XC1, XC2)
  
  if(!is.null(reciprocity)){
    nr <- 1+ ((dim(model.matrix(reciprocity))[2]-1)/nact)
    X4 <- matrix(rep(NA, n*nr), ncol=nr)
    X4[,1] <- 1
    for (i in 2:nr){
      Rtmp <- model.matrix(reciprocity)[,(2+(i-2)*nact): (1+(i-1)*nact)]
      X4[,i] <- Rtmp[row(Rtmp)!=col(Rtmp)]
    }
  } else {
    nr <- 1
    X4 <- matrix(rep(1, n*nr), ncol=nr)
  }  
  return(list(y=y, X=X, X1=X1, X2=X2, X3=X3, X4=X4, nact=nact, ns=ns, nre=nre, nd=nd, nr=nr))
}

llp2 <- function (y, X, X4, beta, g4, M, My, R, rInd){
  eXb <- exp(X%*%beta)
  M[row(M)!=col(M)] <- eXb
  R[row(R)!=col(R)] <- exp(X4%*%g4)
  Mt <- t(M)
  K <- 1 + M + Mt + M*(Mt*R)
  My[row(My)!=col(My)] <- (1-y) + y*eXb
  Pl <-((My*t(My))*((1-rInd)+ (rInd*R)))/K
  ll <- sum(log(Pl[lower.tri(Pl)])) 
  return(ll)
}

effectiveEst <- function(vec)
{
  l <- length(vec)
  lb <- l/25
  vecM <- as.vector(rep(NA, l/lb))
  varvec <- var(vec)
  for (i in 1:(l/lb)){
    vecM[i] <- mean(vec[(((i-1)*lb)+1): (i*lb)])
  }
  return((l*varvec)/(var(vecM)*lb))
}