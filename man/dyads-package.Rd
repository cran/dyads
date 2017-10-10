\name{dyads-package}
\alias{dyads-package}
\alias{dyads}
\docType{package}
\title{
dyads
}
\description{
Package for Dyadic Network Analysis.
}
\details{


\packageDESCRIPTION{dyads}
\packageIndices{dyads}
Package for Dyadic Network Analysis.
}
\author{
Bonne J.H. Zijlstra
Maintainer: Bonne J.H. Zijlstra <B.J.H.Zijlstra@uva.nl>
}
\references{
Zijlstra, B.J.H., Duijn, M.A.J. van, and  Snijders, T.A.B. (2009). MCMC estimation for the $p_2$ network regression model with crossed random effects. \emph{British Journal of Mathematical and Statistical Psychology}, 62, 143-166.
}


\keyword{p2 model}

\examples{

# create a very small network with covariates for illustrative purposes
S <- c(1,0,1,0,1,1,0,1,0,1)
REC <- (S*-1)+1
D1 <- matrix(c(0,1,0,1,0,1,0,1,0,1,
              0,0,0,1,0,1,0,1,0,1,
              1,1,0,0,1,0,0,0,0,0,
              1,1,1,0,1,0,0,0,0,1,
              1,0,1,0,0,1,1,0,1,0,
              0,0,0,0,0,0,1,1,1,1,
              0,0,0,0,0,1,0,1,0,1,
              1,0,0,0,0,1,1,0,1,1,
              0,1,0,1,0,1,0,1,0,0,
              1,0,1,1,1,0,0,0,0,0), ncol=10)
D2 <- abs(matrix(rep(S,10), byrow = FALSE, ncol= 10) -
            matrix(rep(REC,10), byrow = TRUE, ncol= 10))
R <- D1*t(D1)
Y <- matrix(c(0,1,1,1,1,1,0,0,1,1,
              0,0,0,1,1,1,0,0,1,0,
              1,1,0,1,1,1,0,0,1,1,
              1,1,1,0,1,1,0,1,1,0,
              1,1,1,1,0,1,1,0,1,1,
              0,1,1,1,1,0,1,1,1,0,
              1,0,1,0,1,1,0,1,0,1,
              0,1,1,1,0,1,1,0,1,1,
              1,0,1,0,1,0,1,1,0,1,
              1,1,1,0,0,1,1,1,1,0), ncol=10) 

# estimate p2 model
p2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D1 + D2, reciprocity= ~ R,
   burnin = 1000, sample = 4000, adapt = 20)
# Notice: burn-in, sample size and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: p2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D, reciprocity= ~ R)
}
