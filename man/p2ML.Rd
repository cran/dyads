\name{p2ML}
\alias{p2ML}

\title{
  MCMC estimates for the (multilevel) p2 model
}
\description{
  Estimates the (multilevel) p2 model parameters,as described in Zijlstra, Van Duijn and Snijders (2006) <doi: 10.1027/1614-2241.2.1.42>.
}
\usage{
p2ML(nets, sender = NULL, receiver = NULL, density =~ 1, reciprocity =~ 1, 
adapt = NULL, burnin = NULL, center = NULL, separate= NULL, seed = NULL) 

}

\arguments{
  \item{nets}{List with n dichotomous dependent directed networks.}
  \item{sender}{Optional matrix with a stacked actor-level sender covariate, corresponding to the n networks. Multiple sender covariates can be added as a formula object, see example below}
  \item{receiver}{Optional matrix with a stacked actor-level receiver covariate, corresponding to the n networks. Multiple receiver covariates can be added as a formula object}
  \item{density}{Optional stacked matrix with a density covariate, with dimensions similar to the n dependent networks. Multiple density covariates can be added as a formula object, see example below}
  \item{reciprocity}{Optional stacked matrix with a symmetric reciprocity covariate, with dimensions similar to the n dependent networks. Multiple reciprocity covariates can be added as a formula object}
  \item{adapt}{Optional number of adaptive sequenses (default is 125).}
  \item{burnin}{Optional specification of number of burn-in iterations (default is 2500).}
  \item{center}{Optional argument for centering predictors (default is TRUE).}
  \item{separate}{Optional argument for estimating separate coefficients for the n dependent networks (default is FALSE).}
  \item{seed}{Optonal specification of random seed (delfault is 1).}
}

\value{
Returns a fitted model of class 2ML, to be opened with the function summary().
}
\references{
Zijlstra, B. J., Van Duijn, M. A., & Snijders, T. A. (2006). The Multilevel p2 Model A random effects model for the analysis of multiple social networks. \emph{Methodology: European Journal of Research Methods for the Behavioral and Social Sciences}, 2(1), 42.
}
\author{
Bonne J.H. Zijlstra \email{b.j.h.zijlstra@uva.nl}
}

\examples{

# create two very small networks with covariates for illustrative purposes
Y1 <- matrix(c(0,1,0,1,0,1,0,1,0,0,
               0,0,1,1,0,1,0,1,0,1,
               1,1,0,0,1,0,0,0,0,0,
               1,1,1,0,1,0,0,0,0,1,
               1,0,1,0,0,1,1,0,1,1,
               0,0,0,0,0,0,1,1,1,1,
               0,0,0,0,0,1,0,1,0,1,
               1,0,0,0,0,1,1,0,1,1,
               0,1,0,1,0,1,0,1,0,0,
               0,0,1,1,1,0,0,0,0,0), ncol=10)
Y2 <- matrix(c(0,0,1,0,1,0,0,1,0,0,
               0,0,0,0,0,0,0,1,1,0,
               0,0,0,1,0,1,0,1,0,1,
               0,0,1,0,0,0,1,1,0,0,
               1,0,0,1,0,0,1,0,0,1,
               0,0,1,0,0,0,1,1,0,0,
               0,1,0,0,1,0,0,0,0,0,
               1,0,1,0,1,1,1,0,0,1,
               0,1,0,1,0,0,0,0,0,0,
               0,1,0,1,0,0,0,1,0,0), ncol=10) 
Y <- list(Y1, Y2)                
Sa1 <- c(1,0,1,0,1,1,0,1,0,1)
Sa2 <- c(1,0,0,1,0,0,1,1,0,1)
Sa <- list(Sa1, Sa2)
Sat <- Rat <-  do.call(plyr::rbind.fill.matrix, Sa)
Sb1 <- c(0,1,1,0,1,0,1,0,1,0)
Sb2 <- c(1,0,1,0,0,1,0,1,0,1)
Sb <- list(Sb1, Sb2)
Sbt <- do.call(plyr::rbind.fill.matrix, Sb)
Da1 <- abs(matrix(rep(Sa1,10), byrow = FALSE, ncol= 10) -
              matrix(rep(Sa1,10), byrow = TRUE, ncol= 10))
Da2 <- abs(matrix(rep(Sa2,10), byrow = FALSE, ncol= 10) -
              matrix(rep(Sa2,10), byrow = TRUE, ncol= 10))
Da <- list(Da1, Da2)
Dat <- do.call(plyr::rbind.fill.matrix, Da)

# estimate p2 model for two networks
M1 <- p2ML(Y, sender= ~ Sat + Sbt,  receiver= ~ Rat, density = ~ Dat, adapt = 10, burnin = 100)
summary(M1)
# Notice: burn-in, and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: 
\dontrun{
M1 <- p2ML(Y,sender= ~ Sat + Sbt,  receiver= ~ Rat, density = ~ Dat)
summary(M1)
}

# estimate p2 model for a single network
M2 <- p2ML(list(Y[[1]]),sender= ~ Sat[1:10,] + Sbt[1:10,],  receiver= ~ Rat[1:10,],
density = ~ Dat[1:10,], adapt = 10, burnin = 100)
summary(M2)
# Notice: burn-in, and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: 
\dontrun{
M2 <- p2ML(list(Y[[1]]),sender= ~ Sat[1:10,] + Sbt[1:10,],  receiver= ~ Rat[1:10,], 
density = ~ Dat[1:10,])
summary(M2)
}
}
\keyword{p2 model}

