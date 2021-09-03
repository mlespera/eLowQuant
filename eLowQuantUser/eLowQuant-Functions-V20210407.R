# MLE functions for Poisson Calibration model (not shown)

# a*log(b) returns zero if both a and b are zero
alogb <- function(a,b){
  a * (a!=0) * log(b+(a==0 & b==0))
}

#Log Likelihood for model not through origin
# betas=c(intercept,slope) starting values; nd=number detected; 
# S=Starting Quantity; nn=number of replicates
# nn=number of technical replicates

Calib.LLik <- function(betas, nd, S, nn) {
  sum(alogb(nd, 1 - exp(-(betas[1] + betas[2]*S))) - (nn - nd)*(betas[1] + betas[2]*S))
}

Calib.LLikbak <- function(betas, nd, S, nn) {
  sum(alogb(nd, 1 - exp(-(betas[1] + betas[2]*S))) - (nn - nd)*(betas[1] + betas[2]*S))
}

Calib.LLik.bak1 <- function(betas, nd, S, nn) {
  sum(nd*log(1 - exp(-(betas[1] + betas[2]*S))) - (nn - nd)*(betas[1] + betas[2]*S))
}

#Log Likelihood for model through origin - use for NTC detects very near zero
# betas=c(slope) starting values; nd=number detected; 
# S=Starting Quantity; nn=number of replicates
# For S_j = 0, if(nd_j = 0) L_j = 1 (LL_j=0) else L_j = 0 (LL_j=-INF) 
CalibOr.LLik <- function(betas, nd, S, nn) {
  bool <- S != 0
  nd <- nd[bool]; nn <- nn[bool]; S <- S[bool]
  sum(alogb(nd, 1 - exp(-(betas*S))) - (nn - nd)*(betas*S))
}


#Log Relative Likelihood interval of level Plevel; solve for roots
# Plevel in (0,1)
# betas=(intercept, slope, Shat) obtained from optimizing CalibS0.LLik not through origin
Calib.LRL <- function(S, betas, nd0, nn0, Plevel){
  Calib.LLik(betas[1:2], nd0, S, nn0) - Calib.LLik(betas[1:2], nd0, betas[3], nn0) - log(Plevel)
}

#Log Relative Likelihood interval of level Plevel; solve for roots
# Plevel in (0,1)
# betas=(slope, Shat) obtained from optimizing CalibOrS0.LLik through origin
CalibOr.LRL <- function(S, betas, nd0, nn0, Plevel){
  CalibOr.LLik(betas[1], nd0, S, nn0) - CalibOr.LLik(betas[1], nd0, betas[2], nn0) - log(Plevel)
}

#Derivative of Log likelihood for model not through origin
Calib.dLLik <- function(betas, nd, S, nn){
  pj <- 1 - exp(-(betas[1]+betas[2]*S)); pj1 <- 1 - pj
  g <- sum(nd * pj1 / pj - (nn - nd))
  g <- c(g, sum(nd * S * pj1 / pj - (nn - nd)*S ))
  return(g)
}

#Derivative of Log likelihood for model through origin
# For S_j = 0, if(nd_j = 0) L_j = 1 else L_j = 0
CalibOr.dLLik <- function(betas, nd, S, nn){
  bool <- S != 0
  nd <- nd[bool]; nn <- nn[bool]; S <- S[bool]
  pj <- 1 - exp(-(betas*S)); pj1 <- 1 - pj
  #  g <- sum(nd * pj1 / pj - (nn - nd))
  g <- c(sum(nd * S * pj1 / pj - (nn - nd)*S ))
  return(g)
}

#Second Derivatives/Hessian of Log likelihood for model not through origin
Calib.ddLLik <- function(betas, nd, S, nn){
  H <- matrix(0, nrow=2, ncol=2)
  pj <- 1 - exp(-(betas[1]+betas[2]*S)); pj1 <- 1 - pj
  H[1,1] <- -sum(nd * pj1 / (pj^2))
  H[2,2] <- -sum(nd * (S^2) * pj1 / (pj^2))
  H[1,2] <- H[2,1] <- -sum(nd * S * pj1 / (pj^2))
  return(H)
}


#Second Derivatives/Hessian of Log likelihood for model through origin
# For S_j = 0, if(nd_j = 0) L_j = 1 else L_j = 0
CalibOr.ddLLik <- function(betas, nd, S, nn){
  bool <- S != 0
  nd <- nd[bool]; nn <- nn[bool]; S <- S[bool]
  pj <- 1 - exp(-(betas*S)); pj1 <- 1 - pj
  H <- -sum(nd * (S^2) * pj1 / (pj^2))
  return(H)
}



#Binomial Log Likelihood
Bin.LLik <- function(nd, nn){ 
  #  sum(nd*log(nd/nn) + (nn - nd)*log((nn-nd)/nn))
  sum(alogb(nd, nd/nn) + alogb((nn - nd), (nn-nd)/nn))
}

#Log Likelihood for model not through origin, plus new observation, (nd0, nn0)
CalibS0.LLik <-  function(beta.S0, nd, S, nn, nd0, nn0){
  betas <- beta.S0[1:2]; S0 <- beta.S0[3]
  Calib.LLik(betas, nd, S, nn) + Calib.LLik(betas, nd0, S0, nn0)
}

#Log Likelihood for model  through origin, plus new observation, (nd0, nn0, S0)
# For S_j = 0, if(nd_j = 0) L_j = 1 else L_j = 0
CalibS0Or.LLik <-  function(beta.S0, nd, S, nn, nd0, nn0){
  bool <- S != 0
  nd <- nd[bool]; nn <- nn[bool]; S <- S[bool]
  betas <- beta.S0[1]; S0 <- beta.S0[2]
  CalibOr.LLik(betas, nd, S, nn) + CalibOr.LLik(betas, nd0, S0, nn0)
}

#Derivative of Log likelihood for model not through origin, plus new observation, (nd0, nn0)
CalibS0.dLLik <- function(beta.S0, nd, S, nn, nd0, nn0){
  g <- vector("numeric", 3)
  betas <- beta.S0[1:2]; S0 <- beta.S0[3]
  pj0 <- 1 - exp(-(betas[1]+betas[2]*S0)); pj01 <- 1 - pj0
  g[1:2] <- Calib.dLLik(betas, nd, S, nn)
  g[1] <- g[1] + nd0 * pj01 / pj0 - (nn0 - nd0)
  g[2] <- g[2] + nd0 * S0 * pj01 / pj0 - (nn0 - nd0) * S0
  g[3] <- nd0 * betas[2] * pj01 / pj0 - (nn0 - nd0) * betas[2]
  return(g)
}

#Derivative of Log likelihood for model through origin, plus new observation, (nd0, nn0, S0)
# For S_j = 0, if(nd_j = 0) L_j = 1 else L_j = 0
CalibS0Or.dLLik <- function(beta.S0, nd, S, nn, nd0, nn0){
  bool <- S != 0
  nd <- nd[bool]; nn <- nn[bool]; S <- S[bool]
  g <- vector("numeric", 2)
  betas <- beta.S0[1]; S0 <- beta.S0[2]
  pj0 <- 1 - exp(-(betas*S0)); pj01 <- 1 - pj0
  g[1] <- CalibOr.dLLik(betas, nd, S, nn)
  g[1] <- g[1] + nd0 * S0 * pj01 / pj0 - (nn0 - nd0) * S0
  g[2] <- nd0 * betas * pj01 / pj0 - (nn0 - nd0) * betas
  return(g)
}


#Second Derivatives/Hessian of Log likelihood for model not through origin, plus new observation, (nd0, nn0)
CalibS0.ddLLik <- function(beta.S0, nd, S, nn, nd0, nn0){
  H <- matrix(0, nrow=3, ncol=3)
  betas <- beta.S0[1:2]; S0 <- beta.S0[3]
  pj0 <- 1 - exp(-(betas[1]+betas[2]*S0)); pj01 <- 1 - pj0
  H[1:2, 1:2] <- Calib.ddLLik(betas, nd, S, nn)
  H[1,1] <- H[1,1] - nd0 * pj01 / (pj0^2)
  H[2,2] <- H[2,2] - nd0 * pj01 * (S0^2) / (pj0^2)
  H[1,2] <- H[1,2] - nd0 * pj01 * S0 / (pj0^2)    
  H[2,1] <- H[1,2]
  H[1,3] <- - nd0 * pj01 * betas[2] / (pj0^2)
  H[3,1] <- H[1,3]
  H[2,3] <- (nd0/pj0 - nn0) - (nd0 * S0 * betas[2] * pj01) / (pj0^2) 
  H[3,2] <- H[2,3]
  H[3,3] <- H[1,3]* betas[2]
  return(H)
}

#Second Derivatives/Hessian of Log likelihood for model through origin, plus new observation, (nd0, nn0)
CalibS0Or.ddLLik <- function(beta.S0, nd, S, nn, nd0, nn0){
  H <- matrix(0, nrow=2, ncol=2)
  betas <- beta.S0[1]; S0 <- beta.S0[2]
  pj0 <- 1 - exp(-(betas*S0)); pj01 <- 1 - pj0
  H[1,1] <- CalibOr.ddLLik(betas, nd, S, nn)
  H[1,1] <- H[1,1] - nd0 * pj01 * (S0^2) / (pj0^2)
  H[1,2] <- (nd0/pj0 - nn0) - (nd0 * S0 * betas * pj01) / (pj0^2) 
  H[2,1] <- H[1,2]
  H[2,2] <- - nd0 * pj01 * (betas^2) / (pj0^2)
  return(H)
}

