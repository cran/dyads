j2model <- function (net, sender = NULL, receiver = NULL , density = NULL, reciprocity = NULL, center = NULL) 
{
  y <- net[row(net)!=col(net)]
  nact <- ncol(net)
  n <- nact*(nact-1)
  if(!is.null(sender)){
    ns <- dim(model.matrix(sender))[2]-1
    XS <- matrix(model.matrix(sender)[,2:(ns+1)])
    X1 <- matrix(rep(NA, n*ns), ncol=ns)
    for (i in 1:ns){
      XS[,i] <- scale(XS[,i], center=center, scale=FALSE)
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
      XRE[,i] <- scale(XRE[,i], center=center, scale=FALSE)
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
      X3[,i] <- scale(Dtmp[row(Dtmp)!=col(Dtmp)], center=center, scale=FALSE)
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
      X4[,i] <- scale(Rtmp[row(Rtmp)!=col(Rtmp)], center=center, scale=FALSE)
    }
  } else {
    nr <- 1
    X4 <- matrix(rep(1, n*nr), ncol=nr)
  }  
  return(list(y=y, X=X, X1=X1, X2=X2, X3=X3, X4=X4, nact=nact, ns=ns, nre=nre, nd=nd, nr=nr))
}



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

b2modelB2ML1l <- function (nets, actor = NULL , density = NULL, center = NULL, separate = NULL, densVar = NULL) 
{
  current.na.action <- options('na.action')
  options(na.action='na.pass')
  nnets <- length(nets)
  net <- blockMatrixDiagonal(nets)
  if (nnets==1 | densVar==FALSE) {nrandd <- 0; separate <- FALSE} else {nrandd <- nnets}
  if (separate == FALSE){overalld <- 1} else {overalld <- 0}
  y <- unlist(lapply(1:nnets, function(x){nets[[x]][row(nets[[x]])!=col(nets[[x]])]}))
  yl <- list()
  for (i in 1:nnets){yl[[i]] <- nets[[i]][row(nets[[i]])!=col(nets[[i]])]}
  nact <- unlist(lapply(1:nnets, function(x){ncol(nets[[x]])}))
  nactcum <- cumsum(nact)
  nactsum <- sum(nact)
  maxact <- max(nact)
  n <- unlist(lapply(1:nnets, function(x){nact[x]*(nact[x]-1)}))
  ncum <- cumsum(n)
  nsum <- sum(n)
  if(!is.null(actor)){
    ns <- dim(model.matrix(actor))[2]-1
    XS <- matrix(model.matrix(actor)[,2:(ns+1)], ncol = ns)
    X1 <- matrix(rep(NA, nsum*ns), ncol=ns)
    for (i in 1:ns){
      XS[,i] <- scale(XS[,i], center=center, scale=FALSE)
      S <- XS[,i]%*%t(rep(1,nactsum))
      X1[,i] <- S[row(S)!=col(S) & (!is.na(net))]
    }
  }  else {
    ns <- 0
    X1 <- NULL
    XS <- NULL
  }
  if(!is.null(actor)){
    nre <- dim(model.matrix(actor))[2]-1
    XRE <- matrix(model.matrix(actor)[,2:(nre+1)], ncol = nre)
    X2 <- matrix(rep(NA, nsum*nre), ncol=nre)
    for (i in 1:nre){
      XRE[,i] <- scale(XRE[,i], center=center, scale=FALSE)
      RE <- rep(1,nactsum)%*%t(XRE[,i])
      X2[,i] <- RE[(row(RE)!=col(RE)) & (!is.na(net))]
    }
  } else {
    nre <- 0
    X2 <- NULL
    XRE <- NULL
  }
  startrow <- ncum - n +1
  colnum <- nactcum -nact
  if((!is.null(density)) & (density != ~1)){
    if (separate == FALSE){
      nd <- 1 + nrandd + ((dim(model.matrix(density))[2]-1)/maxact)
      X3 <- matrix(rep(0, nsum*nd), ncol=nd)
      X3[,1] <- 1
    } else {      
      nd <- nrandd + ((dim(model.matrix(density))[2]-1)/maxact)
      X3 <- matrix(rep(0, nsum*nd), ncol=nd)
    }
    if (nrandd > 0){
      for (k in 1:nrandd){
        X3[startrow[k]:ncum[k],k+overalld] <- 1 
        X3[,k+overalld] <- X3[,k+overalld]
      }
    }
    for (i in 2:(nd+(1-overalld)-nrandd)){
      for (j in 1:nnets){
        Dtmp <- model.matrix(density)[((nactcum-nact)[j]+1):nactcum[j],(2+(i-2)*maxact): (1+(i-2)*maxact+nact[j])]
        X3[startrow[j]:ncum[j],i-(1-overalld)+nrandd] <- scale(Dtmp[row(Dtmp)!=col(Dtmp)], center=center, scale=FALSE)
      }
    }
  } else {
    if (separate == FALSE){
      nd <- nrandd +1
      X3 <- matrix(rep(0, nsum*nd), ncol=nd)
      X3[,1] <- 1
    } else {
      nd <- nrandd
      X3 <- matrix(rep(0, nsum*nd), ncol=nd)
    }
    if (nrandd > 0){
      for (k in 1:nrandd){
        X3[startrow[k]:ncum[k],k+overalld] <- 1 
        X3[,k+overalld] <-X3[,k+overalld]
      }  
    }
  }
  X3l <- list()
  for (k in 1:nnets){
    X3l[[k]] <- cbind(X3[startrow[k]:ncum[k], ])
  }  
  XC1 <- matrix(rep(0, nsum*nactsum), ncol=nactsum)
  XC2 <- matrix(rep(0, nsum*nactsum), ncol=nactsum)
  for (k in 1:nnets){
    for (i in 1:nact[k]){
      SC <- matrix(rep(0, nact[k]*nact[k]), ncol=nact[k])
      SC[i,] <-1
      XC1[startrow[k]:ncum[k],i+ colnum[k]] <- SC[row(SC)!=col(SC)]
      REC <- matrix(rep(0, nact[k]*nact[k]), ncol=nact[k])
      REC[,i] <-1
      XC2[startrow[k]:ncum[k],i+ colnum[k]] <-REC[row(REC)!=col(REC)]
    }
  }
  X <- cbind(X1, X2, X3, XC1, XC2)
  Xl <- list()
  for (k in 1:nnets){
    Xl[[k]] <- cbind(X1[startrow[k]:ncum[k],], X2[startrow[k]:ncum[k],], X3[startrow[k]:ncum[k], ], XC1[startrow[k]:ncum[k], (colnum[k]+1):nactcum[k]], XC2[startrow[k]:ncum[k],(colnum[k]+1):nactcum[k]])
  }  
  nr <- 0
  nrandr <- 0
  X4l <- list()
  netnums <- matrix(rep(NA, nactsum), ncol=1)
  for (k in 1:nnets){
    netnums[(nactcum[k]-nact[k]+1):nactcum[k],1] <- k
  }  
  return(list(y=y, X=X, X1=X1, X2=X2, X3=X3, yl=yl, Xl=Xl, X3l=X3l, X4l=X4l, XS=XS, XRE=XRE, nrows=n, nact=nactsum, netnums=netnums, nactnets=nact, ns=ns, nre=nre, nd=nd, nr=nr, nrandd= nrandd, nrandr=nrandr))
  options(na.action='current.na.action')
}

llj2 <- function (y, X, X4, beta, g4, M, My, R, cSign){
  eXb <- exp(X%*%beta)
  M[row(M)!=col(M)] <- eXb
  R[row(R)!=col(R)] <- exp(X4%*%g4)
  Mt <- t(M)
  MM <- M*Mt
  K <- 1 + M + Mt + MM
  My[row(My)!=col(My)] <- (1-y) + y*eXb
  Pla <- (My*t(My))
  G <- 1+ MM + M*R + Mt*R
  H <- R -1
  C <- (-1*sqrt(G*G-(4*H*H*MM))+G)/(2*H)
  C[is.nan(C)] <- 0
  Pl <- (Pla + (cSign*C))/K
  ll <- sum(log(Pl[lower.tri(Pl)]))  
  return(ll)
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

llb2MLC <- function (y, Y, X, beta, M, My){
  eXb <- exp(X%*%beta)
  M[(row(M)!=col(M)) & (!is.na(Y))] <- eXb
  K <- 1 + M
  My <- (1-Y)/K + Y*M/K
  return(log(My))
}

lpC <- function(var, list, grp, vec) {
  for (i in 1:max(grp)){
    vec[i] <-  sum(Rfast::dmvnorm(var[grp==i,], mu= c(0,0), sigma= list[[i]], logged = TRUE))
  }    
  return(vec)
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

effectiveEst3 <- function(vec)
{
  vec <- vec[!is.na(vec)]
  l <- length(vec)
  lb <- trunc(sqrt(l))
  vecM <- as.vector(rep(NA, l/lb))
  varvec <- var(vec)
  for (i in 1:(l/lb)){
    vecM[i] <- mean(vec[(((i-1)*lb)+1): (i*lb)])
  }
  return((l*varvec)/(var(vecM)*lb))
}

blockMatrixDiagonal<-function(...){  
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(NA,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}

Ccenter <- function(var,grp) {
  for (i in 1:max(grp)){
    var[grp==i,] <-  scale(var[grp==i,], center = TRUE, scale = FALSE)
  }    
  return(as.vector(var))
}

summary.b2ML <- function(z){
  MCMCsims <- z$MCMCsims[ , c(1, if (z$drandd>0) 1+z$drandd, 
                              if (z$ns>0) {(1+z$drandd+1):(1+z$drandd+z$ns)},
                              if (z$separate == FALSE) {1+z$drandd+z$ns+z$nre+1} else {(1+z$drandd+z$ns+z$nre+1):(1+z$drandd+z$ns+z$nre+z$nrandd)}, 
                              if (z$nd>0) {if (z$separate == FALSE) {(1+z$drandd+z$ns+z$nre+z$nrandd+2):(1+z$drandd+z$ns+z$nre+z$nrandd+1+z$nd)} else {(1+z$drandd+z$ns+z$nre+z$nrandd+1):(1+z$drandd+z$ns+z$nre+z$nrandd+z$nd)}})]
  output.matrix <- matrix(NA, dim(MCMCsims)[2], 10)
  colnames(output.matrix) <- c("Estimate", "SE", "Q.05", "Q2.5", "Q25", "Q50", "Q75", "Q97.5", "Q99.5", "Neff")
  rownames(output.matrix) <- colnames(MCMCsims)
  output.matrix[,1] <- colMeans(MCMCsims, na.rm = T) 
  output.matrix[,2] <- sqrt(diag(cov(MCMCsims, use = "pairwise.complete.obs")))
  output.matrix[,3:9] <- t(apply(MCMCsims, 2, quantile, probs = c(0.5, 2.5, 25, 50, 75, 97.5, 99.5)/100, na.rm=T))
  output.matrix[,10] <- t(apply(MCMCsims, 2, effectiveEst3))
  return(round(output.matrix, digits=3))  
}

