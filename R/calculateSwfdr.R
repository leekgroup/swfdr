#' Calculate the science-wise FDR (swfdr)
#' 
#' @param pValues Numerical vector of p-values
#' @param truncated Vector of 0s and 1s with indices corresponding to those in pValues; 1 indicates that the p-values is truncated, 0 that it is not truncated
#' @param rounded Vector of 0s and 1s with indices corresponding to those in pValues; 1 indicates that the p-values is rounded, 0 that it is not rounded
#' @param pi0 Initial prior probability that a hypothesis is null (default is 0.5)
#' @param alpha Initial value of parameter alpha from Beta(alpha, beta) true positive distribution (default is 1)
#' @param beta Initial value of parameter beta from Beta(alpha, beta) true positive distribution (default is 50)
#' @param numEmIterations The number of EM iterations (default is 100)
#' 
#' @return pi0 Final value of prior probability - estimated from EM - that a hypothesis is null, i.e. estimated swfdr
#' @return alpha Final value of parameter alpha - estimated from EM - from Beta(alpha, beta) true positive distribution
#' @return beta Final value of parameter beta - estimated from EM - from Beta(alpha, beta) true positive distribution
#' @return z Indicator vector of 0s and 1s - estimated from EM - with 1s indicating p-values from the null distribution
#' @return n0 Expected number of null p-values - estimated from EM - between certain cutpoints (0.005, 0.015, 0.025, 0.035, 0.045, 0.051)
#' @return n Number of p-values between certain cutpoints (0.005, 0.015, 0.025, 0.035, 0.045, 0.051)
#' 
#' @import stats4
#' 
#' @export
calculateSwfdr = function(pValues,truncated,rounded,pi0 = 0.5,alpha=1,beta=50,numEmIterations=100){
  pp = pValues
  tt = truncated
  rr = rounded
  
  
  ll = function(a,b){
    tmp1 = rep(0,length(pp))
    tmp1[tt==0 & rr==0] = log(dbeta(pp[tt==0 & rr==0],a,b)/pbeta(0.05,a,b))
    tmp1[tt > 0 & rr==0] = log(pbeta(pp[tt > 0 & rr==0],a,b)/pbeta(0.05,a,b))
    tmp1 = -sum((1-z)*tmp1)
    
    probvec = (pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),a,b) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),a,b))/pbeta(0.05,a,b)
    probvec = rev(probvec)
    tmp2 = sum(-n1*log(probvec))
    
    
    return(tmp1 + tmp2)
  }
  
  n = table(cut(pp[rr > 0],c(-0.01,0.005,0.015,0.025,0.035,0.045,0.051)))
  
  
  for(i in 1:numEmIterations){
    
    ## E-step
    
    probvec1 = (pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),alpha,beta) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),alpha,beta))/pbeta(0.05,alpha,beta)
    probvec1 = rev(probvec1)
    probvec0 = c(0.005,0.01,0.01,0.01,0.01,0.005)*20
    
    pij0 = pi0*probvec0/(probvec0*pi0 + probvec1*(1-pi0))
    n0 = n*pij0
    n1 = n - n0
    
    z = rep(0,length(pp))
    z[tt == 0 & rr ==0] <- pi0*20/(pi0*20 + (1-pi0)*dbeta(pp[tt==0 & rr == 0],alpha,beta)/pbeta(0.05,alpha,beta))
    z[tt > 0 & rr ==0] <- pi0*20*pp[tt > 0 & rr ==0]/(pi0*20*pp[tt > 0 & rr==0] + (1-pi0)*pbeta(pp[tt > 0 & rr==0],alpha,beta)/pbeta(0.05,alpha,beta))
    
    ## M-step
    
    pi0 = (sum(n0) + sum(z))/(sum(n) + sum(rr == 0))
    tmp = mle(ll,start=list(a=0.05,b=100),lower=c(0.001,1),upper=c(1,500),method="L-BFGS-B")
    alpha = coef(tmp)[1]
    beta = coef(tmp)[2]
  }
  return(list(pi0 = pi0, alpha=alpha, beta = beta, z=z,n0=n0,n=n))
}




